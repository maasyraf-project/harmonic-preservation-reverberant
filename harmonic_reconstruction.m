clc
clear all 
close all

%% import and auralize full audio
[h, fs] = audioread("sample_rev_high_min90.wav");
hl = h(:,1);
hr = h(:,2);
t = (0:size(h,1)-1)/fs;

[l, fs] = audioread("sample_rev_low_min90.wav");
l = [l' zeros(2, length(h)-length(l))]';
ll = l(:,1);
lr = l(:,2);

[x, fs] = audioread("sample_anechoic_min90.wav");
x = [x' zeros(2, length(h)-length(x))]';
xl = x(:,1);
xr = x(:,2);

% conduct VAD on clean speech
Nw = 0.03*fs;
Nsh = 0.015*fs;
yl = G729(xl, fs, Nw, Nsh);
yr = G729(xr, fs, Nw, Nsh);

% segmentting and windowing signal
xl = buffer(xl, Nw, Nsh, 'nodelay');
xr = buffer(xr, Nw, Nsh, 'nodelay');

win = hanning(0.03*fs);

VADPar = InitVADPar(fs,Nw,Nsh);

nfft = 4096;
f = (fs/2) * linspace(0,1,nfft/2+1);

syl = zeros(1, size(x,1) + nfft);
for i = 1:size(xl,2)
    xl(:,i) = xl(:,i) .* win;

    % conduct fft
    SXL = fft(xl(:,i), nfft);
    XL = abs(SXL(1:numel(f)));

    % harmonic reconstruction

    % plot
    if i == 100
        figure()
        plot(f, XL);
    end

    % conduct ifft
    sxl = real(ifft(SXL, nfft));
    outIndex = i*Nsh+1:i*Nsh+nfft;

    syl(outIndex) = syl(outIndex) + sxl';

end
