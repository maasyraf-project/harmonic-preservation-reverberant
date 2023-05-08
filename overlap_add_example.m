Fs = 4000; % sampling rate in Hz
Nsig = 128; % signal length in samples
L = 32; % filter length in taps
fc = 600; % cutoff frequency in Hz
M = 32; % window length
Nfft = 2^(ceil(log2( M + L - 1))); % FFT Length
R = M/2; % Hop Size
Nframes = floor( Nsig/R ); % No. of frames
% Generate a test signal containing five impulses
sig = zeros(1,Nsig); period = round(Nsig/5);
sig(1:period:Nsig) = ones(size(1:period:Nsig));

% *** design a lowpass filter using the window method
epsilon = .0001; % avoids 0 / 0
nfilt = (-(L-1)/2:(L-1)/2) + epsilon;

hideal = sin( 2*pi*fc*nfilt/Fs ) ./ (pi*nfilt); % sinc fn
h = hamming( L )' .* hideal; % windowed sinc fn
hzp = [h zeros(1,Nfft-L)]; % Zero-pad h to FFT size
H = fft(hzp); % Filter frequency response
y = zeros(1,Nsig + Nfft); % Allocate the output+’ringing’ vector

figure(1);clf
subplot(211); stem(sig); axis( [0 length(sig) -0.2 1.2]);
subplot(212); stem(h); grid; axis( [0 length(h) -0.2 0.5]);


% *** Overlap add loop, with plots
figure(2);
subplot(411); stem(sig); axis( [0 length(sig) -0.2 1.2]);
title('signal'); ylabel('amplitude');
window = hanning(M)'; % Signal window (could also be rect)

for m = 0:(Nframes-1)
index = m*R+1:min(m*R+M,Nsig); % index for the mth frame
xm = sig(index) .* window(index-m*R); % windowed mth frame
xmzp = [ xm zeros(1,Nfft-length(xm))]; % zero pad the signal
Xm = fft( xmzp );
Xfilt = Xm .* H; % freq domain multiplication
xmfilt = real(ifft(Xfilt)); % inverse transform
outindex = m*R+1:(m*R+Nfft); %
y(outindex) = y(outindex) + xmfilt; % overlap add


xmplot = [ zeros(1,m*R) xm zeros(1,Nsig-(m*R+M)) ];
winplot = [ zeros(1,m*R) window zeros(1,Nsig-(m*R+M)) ];
xmfiltplot = [ zeros(1,m*R) xmfilt zeros(1,Nsig-(m*R+Nfft)) ];
subplot(412);stem(xmplot); axis( [0 length(sig) -0.2 1.2]);
hold on; plot(winplot,'--'); hold off
title('windowed frame');
subplot(413);plot(xmfiltplot); grid; axis( [0 length(sig) -0.1 0.4]);
title('windowed frame - filtered');
subplot(414);plot(y); grid; axis( [0 length(sig) -0.1 0.4]);
title('overlapp add buffer'); xlabel('samples');
cmd = sprintf('print -deps ola%d',m); disp(cmd); eval(cmd);
end
