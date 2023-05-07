clc
clear all
close all
addpath(genpath('additional-packages'));

%% import and auralize full audio
[x, fs] = audioread("sample_anechoic_min90.wav");
tx = (0:size(x,1)-1)/fs;

[l, fs] = audioread("sample_rev_low_min90.wav");
tl = (0:size(l,1)-1)/fs;

[h, fs] = audioread("sample_rev_high_min90.wav");
th = (0:size(h,1)-1)/fs;

figure(1)
subplot(3,2,1)
plot(tx, x(:,1))
xlim([0 3.2])
ylim([-0.15 0.15])
hold on

subplot(3,2,2)
plot(tx, x(:,2))
xlim([0 3.2])
ylim([-0.15 0.15])
hold on

subplot(3,2,3)
plot(tl, l(:,1))
xlim([0 3.2])
ylim([-0.15 0.15])

subplot(3,2,4)
plot(tl, l(:,2))
xlim([0 3.2])
ylim([-0.15 0.15])

subplot(3,2,5)
plot(th, h(:,1))
xlim([0 3.2])
ylim([-0.15 0.15])

subplot(3,2,6)
plot(th, h(:,2))
xlim([0 3.2])
ylim([-0.15 0.15])

%% speech analysis
% conduct VAD on clean speech
Nw = 512;
Nsh = 256;
yl = 0.05*G729(x(:,1), fs, Nw, Nsh);
yr = 0.05*G729(x(:,2), fs, Nw, Nsh);

figure(1)
subplot(3,2,1)
plot(tx, yl)
subplot(3,2,2)
plot(tx, yr)

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
