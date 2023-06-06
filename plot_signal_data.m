clc
clear all
close all

%% load all necessary packages and files
addpath(genpath('additional-packages'))

% load anechoic signal
[an, fs] = audioread("sample_anechoic_min90.wav");

% load low and high reverberant signal
[lrev, fs] = audioread("sample_rev_low_min90.wav");
[hrev, fs] = audioread("sample_rev_high_min90.wav");

% concatenate all loaded data
raw_data = {an, lrev, hrev};

% calculate start- and end-point on the audio
idx_start = [];
idx_end = [];
for i = 1:length(raw_data)
    l_start = find(raw_data{i}(:,1), 1, "first");
    r_start = find(raw_data{i}(:,2), 1, "first");
    
    l_end = find(raw_data{i}(:,1), 1, "last");
    r_end = find(raw_data{i}(:,2), 1, "last");
    
    idx_start = [idx_start; [l_start, r_start]]; 
    idx_end = [idx_end; [l_end, r_end]]; 
end

idx_start = min(idx_start(:));
idx_end = min(idx_end(:));
for i = 1:length(raw_data)
    l_sig = raw_data{i}(idx_start:idx_end,1);
    r_sig = raw_data{i}(idx_start:idx_end,2);
    raw_data{i} = [l_sig, r_sig]; 
end

len_data = length(idx_start:idx_end);
t = 0:(1/fs):(len_data-1)/fs;

%% plot the data
figure(1)
subplot(3,2,1)
plot(t, raw_data{1}(:,1))
ylim([-0.2 0.2])
ylabel("Amplitude")
title("Waveform on Left Channel")
grid on
grid minor

subplot(3,2,2)
plot(t, raw_data{1}(:,2))
ylim([-0.2 0.2])
ylabel("Amplitude")
title("Waveform on Right Channel")
grid on
grid minor

subplot(3,2,3)
plot(t, raw_data{2}(:,1))
ylim([-0.2 0.2])
ylabel("Amplitude")
grid on
grid minor

subplot(3,2,4)
plot(t, raw_data{2}(:,2))
ylim([-0.2 0.2])
ylabel("Amplitude")
grid on
grid minor

subplot(3,2,5)
plot(t, raw_data{3}(:,1))
ylim([-0.2 0.2])
xlabel("Time (s)")
ylabel("Amplitude")
grid on
grid minor

subplot(3,2,6)
plot(t, raw_data{3}(:,2))
ylim([-0.2 0.2])
xlabel("Time (s)")
ylabel("Amplitude")
grid on
grid minor

%% plot the spectrogarm
figure(2)
subplot(3,2,1)
spectrogram(raw_data{1}(:,1), 512, 256, 512, 'yaxis', fs)
title("Spectrogram on Right Channel")

subplot(3,2,2)
spectrogram(raw_data{1}(:,2), 512, 256, 512, 'yaxis', fs)
title("Spectrogram on Right Channel")

subplot(3,2,3)
spectrogram(raw_data{2}(:,1), 512, 256, 512, 'yaxis', fs)
title("Spectrogram on Right Channel")

subplot(3,2,4)
spectrogram(raw_data{2}(:,2), 512, 256, 512, 'yaxis', fs)
title("Spectrogram on Right Channel")

subplot(3,2,5)
spectrogram(raw_data{3}(:,1), 512, 256, 512, 'yaxis', fs)
title("Spectrogram on Right Channel")

subplot(3,2,6)
spectrogram(raw_data{3}(:,2), 512, 256, 512, 'yaxis', fs)
title("Spectrogram on Right Channel")