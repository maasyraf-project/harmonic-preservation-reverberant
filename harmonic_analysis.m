clc
clear all
close all

%% import and auralize full audio
[x, fs] = audioread("full-speech-sample.wav");



%% speech analysis
% conduct VAD on clean speech
Nw = 512;
Nsh = 256;
y = G729(x, fs, Nw, Nsh);

figure(1)
plot(x)
hold on
plot(y)

% select a voiced segment from clean speech
x_voiced = x(y==1);

shortFrame = 30e-3 * fs;
longFrame = 80e-3 * fs;

idx = 8;

idx_start_short = 1 + idx*shortFrame;
idx_end_short = idx_start_short + shortFrame - 1;
x_voiced_shortFrame = x_voiced(idx_start_short:idx_end_short);

idx_start_long = 1 + idx*longFrame;
idx_end_long = idx_start_long + longFrame - 1;
x_voiced_longFrame = x_voiced(idx_start_long:idx_end_long);

% make fft observation
F_short = fs/2 * linspace(0,1,shortFrame/2+1);
X_short = fft(x_voiced_shortFrame, shortFrame);
X_short = abs(X_short(1:shortFrame/2+1));

F_long = fs/2 * linspace(0,1,longFrame/2+1);
X_long = fft(x_voiced_longFrame, longFrame);
X_long = abs(X_long(1:longFrame/2+1));

figure(2)
subplot(5,1,1)
plot(F_short, X_short)
xlim([0 3000])

figure(3)
subplot(5,1,1)
plot(F_long, X_long)
xlim([0 3000])

% estimate F0 of voiced audio sample
