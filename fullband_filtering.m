clc
clear all
close all

%% preparation
% load all packages
addpath(genpath('additional-packages'));
target_fs = 16e3;

startTwoEars  % optional

% load impulse response
irFilename = 'impulse_responses/surrey_cortex_rooms/UniS_Room_D_BRIR_16k.sofa';
irName = split(irFilename,'/');

irFilenameDir = char(strcat('\', irName(1),'\', irName(2),'\', irName(3)));
irDir = char(strcat(pwd, '\additional-packages\TwoEars\BinauralSimulator\tmp',irFilenameDir));

irName = char(irName(end));
irName = irName(1:end-5);

if ~exist(irDir, "file")
    irFilename = db.downloadFile(irFilename);
    ir = SOFAload(irFilename);
else
    ir = SOFAload(strcat(irName,'.sofa'));
end

ir_left = squeeze(ir.Data.IR(1, 1, :));
ir_right = squeeze(ir.Data.IR(1, 2, :));

% import input audio
[x, fs] = audioread("sample_mono_clean_speech.wav");
startSpeech = find(x, 1, "first");
endSpeech = find(x, 1, "last");

x = x(startSpeech:endSpeech);

% create reverberant audio
y = filter(ir_right, 1, x);

%% extract LP-residual
% split audio into several frames
lenFrames = 0.032*fs;
lenOverlap = 0.5*lenFrames;
nFrames = floor((length(y) - lenOverlap)/lenOverlap) + 1;

estl = zeros(1, length(y));
yres = zeros(1, length(y));
for n = 1:nFrames
    % define start and end index
    idxStart = 1 + (n-1) * (lenFrames - lenOverlap);
    idxEnd = idxStart + lenFrames -1;

    if idxEnd > length(y)
        idxEnd = length(y);
    end

    l_frame = y(idxStart:idxEnd);

    % apply LPC
    la = lpc(l_frame, 12);
    estl_frame = filter([0 -la(2:end)], 1, l_frame);
    resl_frame = l_frame - estl_frame;

    estl(idxStart:idxEnd) = estl_frame;
    yres(idxStart:idxEnd) = resl_frame;

end

%% apply inverse filtering
filterOrder = 250;
stepSize = 0.2*filterOrder;
blockSize = stepSize + filterOrder - 1;

lgh = [1 zeros(1, length((1:filterOrder))-1)]';
% lgh = ones(1, length(1:filterOrder))';
lgh = lgh./sqrt(sum(abs(lgh).^2));
LGh = fft([lgh; zeros(stepSize-1, 1)]);

nIter = 150;
mu = 3e-6;

lzKurt = zeros(nIter, 1);
lzr2 = zeros(blockSize, 1);
lzr3 = zeros(blockSize, 1);
lzr4 = zeros(blockSize, 1);
lzr5 = zeros(blockSize, 1);
lzr6 = zeros(blockSize, 1);

hos = 3;

tic
for m = 1:nIter
    lenh_fine = zeros(blockSize, 1);
    for n = 1:stepSize:length(yres)
        %idxStart = 1 + (n-1) * (lenFrames - lenOverlap);
        %idxEnd = idxStart + lenOverlap - 1;

        if n+stepSize-1 > length(yres)
            li = [yres(n:length(yres)) zeros(1, 13)];
        else
            li = yres(n:n+stepSize-1);
        end

        % processing left signal
        lenh_fine(1:filterOrder-1) = lenh_fine(end-filterOrder+2:end);
        lenh_fine(filterOrder:end) = li;

        LENH_FINE = fft(lenh_fine);
        cLENH_FINE = conj(LENH_FINE);
        lzrn = ifft(LGh.*LENH_FINE);

        lzrn(1:filterOrder-1) = 0;
        lzr2(filterOrder:end) = lzrn(filterOrder:end).^2;
        lzr3(filterOrder:end) = lzrn(filterOrder:end).^3;
        lzr4(filterOrder:end) = lzrn(filterOrder:end).^4;
        lzr5(filterOrder:end) = lzrn(filterOrder:end).^5;
        lzr6(filterOrder:end) = lzrn(filterOrder:end).^6;
        
        if hos == 3
            LZ2 = mean(lzr2(filterOrder:end));
            LZ3 = mean(lzr3(filterOrder:end));

            LZ2t = fft(lzr2).*cLENH_FINE;
            LZRN = fft(lzrn).*cLENH_FINE;

            lzKurt(m) = max(lzKurt(m), LZ3/(LZ2^(3/2) + eps));

            LgJ = 3 * (LZ2 * LZ2t - LZ3*LZRN) / (LZ2^3 + 1e-20);
        end

        if hos == 4
            LZ2 = mean(lzr2(filterOrder:end));
            LZ4 = mean(lzr4(filterOrder:end));

            LZ3 = fft(lzr3).*cLENH_FINE;
            LZRN = fft(lzrn).*cLENH_FINE;

            lzKurt(m) = max(lzKurt(m), LZ4/(LZ2^2 + eps));

            LgJ = 4 * (LZ2 * LZ3 - LZ4*LZRN) / (LZ2^3 + 1e-20);
        end

        if hos == 5
            LZ2 = mean(lzr2(filterOrder:end));
            LZ5 = mean(lzr5(filterOrder:end));

            LZ4 = fft(lzr4).*cLENH_FINE;
            LZRN = fft(lzrn).*cLENH_FINE;

            lzKurt(m) = max(lzKurt(m), LZ5/(LZ2^(2.5) + eps));

            LgJ = 5 * (LZ2 * LZ4 - LZ5*LZRN) / (LZ2^(3.5) + 1e-20);
        end

        if hos == 6
            LZ2 = mean(lzr2(filterOrder:end));
            LZ6 = mean(lzr6(filterOrder:end));

            LZ5 = fft(lzr5).*cLENH_FINE;
            LZRN = fft(lzrn).*cLENH_FINE;

            lzKurt(m) = max(lzKurt(m), LZ6/(LZ2^3 + eps));

            LgJ = 6 * (LZ2 * LZ5 - LZ6*LZRN) / (LZ2^4 + 1e-20);
        end
        
        % LGh = LGh + mu*LgJ / stepSize;
        LGh = LGh + (1/(eps+norm(LgJ)^2))*LgJ / stepSize;

        LGh=LGh./sqrt(sum(abs(LGh).^2)/blockSize);

    end

end
toc

lgh = ifft(LGh);
ln = filter(lgh, 1, y);

imp = filter(lgh, 1, ir_right);

figure
plot(ir_right)
hold on
plot(imp)

figure
plot(lzKurt)

figure
plot(lgh)

figure
subplot(311)
plot(x)
subplot(312)
plot(y)
subplot(313)
plot(ln)

%% Spectral Substraction
wlen = 64e-3*fs;
hop = 8e-3*fs;
nfft = 1024;
[Sz, f_stft, t_stft] = stft(ln, fs, Window=hamming(wlen), OverlapLength=hop, FFTLength=nfft);
Pz=abs(Sz).^2;
phz=angle(Sz);

ro_w=7;
a_w=5;
i_w=-ro_w:15;
wS=(i_w+a_w)/(a_w^2).*exp(-.5*(i_w/a_w+1).^2);
wS(i_w<-a_w)=0;
figure;plot(i_w,wS)

gS=.32;
epS=1e-3;
Pl=gS*filter(wS,1,Pz.').';

P_att=(Pz-Pl)./Pz;
P_att(P_att<epS)=epS;

Pxh=Pz.*P_att;
Sxh=sqrt(Pxh).*exp(1i*phz);

nu1=.05;
nu2=4;
Ez=sum(abs(Sz).^2, 1)/size(Sz,1);
Exh=sum(abs(Sxh).^2, 1)/size(Sxh,1);

P_att=ones(size(Sz,2),1);
P_att(Ez<nu1 & Ez./Exh>nu2)=1e-3;
Sxh=Sxh*spdiags(P_att,0,length(P_att),length(P_att));

[xh, t_istft] = istft(Sxh, fs, Window=hamming(wlen), OverlapLength=hop, FFTLength=nfft);
xh=real(xh)/max(abs(xh));

%% Plot


figure; 
subplot(221); spectrogram(x(1:fs*2),hamming(400),200,'yaxis');colorbar off
title('Clean speech, x(t)')
subplot(223); spectrogram(y(1:fs*2),hamming(400),200,'yaxis');colorbar off
title('Reverbarant speech, y(t)')
subplot(222); spectrogram(ln(1:fs*2),hamming(400),200,'yaxis');colorbar off
title('Inverse filtered speech, z(t)')
subplot(224); spectrogram(xh(1:fs*2),hamming(400),200,'yaxis');colorbar off
title('Final processed speech, xh(t)')