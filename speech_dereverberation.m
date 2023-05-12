clc
clear all
close all
addpath(genpath('additional-packages'));
eps = 1e-18;

% load input audio
[x, fs] = audioread("full-speech-sample.wav");
[ref, fs] = audioread("sample_anechoic_min90.wav");
tsig = 0:(1/fs):7;
debug = 1;

%ref(ref==0) = eps;

% load impulse response
irFilename = 'impulse_responses/surrey_cortex_rooms/UniS_Room_D_BRIR_16k.sofa';
irName = split(irFilename,'/');

irFilenameDir = char(strcat('\', irName(1),'\', irName(2),'\', irName(3)));
irDir = char(strcat(pwd, '\additional-packages\two-ears\BinauralSimulator\tmp',irFilenameDir));

irName = char(irName(end));
irName = irName(1:end-5);

if ~exist(irDir, "file")
    irFilename = db.downloadFile(irFilename);
    ir = SOFAload(irFilename);
else
    ir = SOFAload(strcat(irName,'.sofa'));
end

ir_left = resample(squeeze(ir.Data.IR(1, 1, :)), fs, ir.Data.SamplingRate);
ir_right = resample(squeeze(ir.Data.IR(1, 2, :)), fs, ir.Data.SamplingRate);

y = [conv(squeeze(ir_left), x) ...
    conv(squeeze(ir_right), x)];

%y(y==0) = eps;

% LP residue analysis of anechoic signal
%refLeft = [ref(:,1); zeros(size(y,1)-size(ref,1),1)];
%refRight = [ref(:,2); zeros(size(y,1)-size(ref,1),1)];

refLeft = [x; zeros(size(y,1)-size(ref,1),1)];
refRight = [x; zeros(size(y,1)-size(ref,1),1)];

[refLeftResidue, refLeft] = lpres(refLeft);
[refRightResidue, refRight] = lpres(refRight);

refLeftResidue(isnan(refLeftResidue)) = 0;
refRightResidue(isnan(refRightResidue)) = 0;

% LP residue analysis of reverberation signal
yLeft = y(:,1);
yRight = y(:,2);

[yLeftResidue, yLeft] = lpres(yLeft);
[yRightResidue, yRight] = lpres(yRight);

yLeftResidue(isnan(yLeftResidue)) = 0;
yRightResidue(isnan(yRightResidue)) = 0;

if debug
    figure
    subplot(421);plot(tsig(1:length(refLeft)),refLeft)
    xlabel('time(s)');title('Left clean speech, xL(t)')
    ylim([-0.2 0.2])
    xlim([0 3])
    subplot(422);plot(tsig(1:length(refLeft)),refRight)
    xlabel('time(s)');title('Right clean speech, xR(t)')
    ylim([-0.2 0.2])
    xlim([0 3])

    subplot(423);plot(tsig(1:length(yLeft)),yLeft)
    xlabel('time(s)');title('Left reverberant speech, yL(t)')
    ylim([-0.2 0.2])
    xlim([0 3])
    subplot(424);plot(tsig(1:length(yRight)),yRight)
    xlabel('time(s)');title('Right reverberant speech, yR(t)')
    ylim([-0.2 0.2])
    xlim([0 3])

    subplot(425);plot(tsig(1:length(refLeftResidue)),refLeftResidue)
    xlabel('time(s)');title('LP residual of left clean speech, x_r(t)')
    ylim([-0.2 0.2])
    xlim([0 3])
    subplot(426);plot(tsig(1:length(refRightResidue)),refRightResidue)
    xlabel('time(s)');title('LP residual of right clean speech, x_r(t)')
    ylim([-0.2 0.2])
    xlim([0 3])

    subplot(427);plot(tsig(1:length(yLeftResidue)),yLeftResidue)
    xlabel('time(s)');title('LP residual of left reverberant speech, y_r(t)')
    ylim([-0.2 0.2])
    xlim([0 3])
    subplot(428);plot(tsig(1:length(yRightResidue)),yRightResidue)
    xlabel('time(s)');title('LP residual of right reverberant speech, y_r(t)')
    ylim([-0.2 0.2])
    xlim([0 3])
end

lseg=400;
theta=1e-6;

refLeftKurt = kurt(refLeftResidue,lseg,theta);
refRightKurt = kurt(refRightResidue,lseg,theta);
yLeftKurt = kurt(yLeftResidue,lseg,theta);
yRightKurt = kurt(yRightResidue,lseg,theta);

if debug
    figure
    subplot(421);plot(tsig(1:length(refLeftResidue)),refLeftResidue)
    xlabel('time(s)');title('LP residual of left clean speech, x_r(t)')
    subplot(422);plot(tsig(1:length(refRightResidue)),refRightResidue)
    xlabel('time(s)');title('LP residual of right clean speech, x_r(t)')

    subplot(423);plot(tsig(1:length(yLeftResidue)),yLeftResidue)
    xlabel('time(s)');title('LP residual of left reverberant speech, y_r(t)')
    subplot(424);plot(tsig(1:length(yRightResidue)),yRightResidue)
    xlabel('time(s)');title('LP residual of right reverberant speech, y_r(t)')


    subplot(425);plot(tsig(1:length(refLeftResidue)),refLeftKurt)
    xlabel('time(s)');title('Kurtosis of LP residual of left clean speech, x_r(t)')
    subplot(426);plot(tsig(1:length(refRightResidue)),refRightKurt)
    xlabel('time(s)');title('Kurtosis of LP residual of right clean speech, x_r(t)')


    subplot(427);plot(tsig(1:length(yLeftResidue)),yLeftKurt)
    xlabel('time(s)');title('Kurtosis of LP residual of left reverberant speech, y_r(t)')
    subplot(428);plot(tsig(1:length(yRightResidue)),yRightKurt)
    xlabel('time(s)');title('Kurtosis of LP residual of right reverberant speech, y_r(t)')


end

% processing left channel
Lf=250;
niter=150;
mu=3e-6;
p=Lf;                       %filter order
ss=Lf;                      %stepsize
bs=ss+p-1;                  %blocksize
% gh=[1;zeros(p-1,1)];
gh=(1:p).';
gh=gh./sqrt(sum(abs(gh).^2));
Gh=fft([gh; zeros(ss-1,1)]);
zkurt=zeros(niter,1);

zr2=zeros(bs,1);
zr3=zeros(bs,1);
zr4=zeros(bs,1);

tic
for m=1:niter
    yrn=zeros(bs,1);
    for k=1:ss:size(ref,1)
        yrn(1:p-1)=yrn(end-p+2:end);
        yrn(p:end)=yLeftResidue(k:k+ss-1);
        
        Yrn=fft(yrn);
        cYrn=conj(Yrn);
        zrn=ifft(Gh.*Yrn);
        
        zrn(1:p-1)=0;
        zr2(p:end)=zrn(p:end).^2;
        zr3(p:end)=zrn(p:end).^3;
        zr4(p:end)=zrn(p:end).^4;
        
        Z2=sum(zr2(p:end));
        Z4=sum(zr4(p:end));
        
        zkurt(m)=max(zkurt(m),Z4/(Z2^2+1e-15)*ss);
        
%         z3y=fft(zr3).*cYrn;
%         zy=fft(zrn).*cYrn;
% 
%         gJ=4*(Z2*z3y-Z4*zy)/(Z2^3+1e-20)*ss;
%         Gh=Gh+mu*gJ;
        
                
        z3y=ifft(fft(zr3).*cYrn);
        z3y(p+1:end)=0;
        zy=ifft(fft(zrn).*cYrn);
        zy(p+1:end)=0;

        gJ=4*(Z2*z3y-Z4*zy)/(Z2^3+1e-20)*ss;
        Gh=Gh+mu*fft(gJ);
        
        
        Gh=Gh./sqrt(sum(abs(Gh).^2)/bs);
    end
end
toc

gh=ifft(Gh);
gh=gh(1:p);
save ifilt.mat gh
figure
plot(zkurt);
xlabel('Iterations'); ylabel('Maximum kurtosis');

load ifilt.mat gh
figure;plot((1:length(gh))/fs,gh);
xlabel('time(s)');title('Inverse filter, \bf{g}');

figure
plot((1:length(squeeze(ir_left)))/fs,squeeze(ir_left)); hold on
plot((1:length(squeeze(ir_left)))/fs,filter(gh,1,squeeze(ir_left)), 'LineWidth',1)
xlabel('time(s)');title('Original vs Inverse filtered RIR');
legend('Original RIR', 'Inverse Filtered RIR');

figure
plot((1:length(squeeze(ir_left)))/fs,filter(gh,1,squeeze(ir_left)), 'LineWidth', 1)
xlabel('time(s)');title('Inverse filtered RIR');

z=filter(gh, 1, yLeft);
% z=z/max(abs(z));
zr=lpres(z);
zkurt=kurt(z,lseg,theta);

figure; 
subplot(311); spectrogram(refLeft,hamming(400),200,'yaxis');colorbar off
title('Clean speech, x(t)')
subplot(312); spectrogram(yLeft,hamming(400),200,'yaxis');colorbar off
title('Reverbarant speech, y(t)')
subplot(313); spectrogram(z,hamming(400),200,'yaxis');colorbar off
title('Inverse filtered speech, z(t)')
