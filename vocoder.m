clc
clear all
close all

%% load all necessary packages and files
addpath(genpath('additional-packages'))

filename = "its_ind_fena_e101_0001_min90.wav";
% load anechoic signal
[an, fs] = audioread(strcat("stimuli\unvocoded\UniS_Anechoic_BRIR_16k\", filename));

% load low and high reverberant signal
[revA, fs] = audioread(strcat("stimuli\unvocoded\UniS_Room_A_BRIR_16k\", filename));
[revB, fs] = audioread(strcat("stimuli\unvocoded\UniS_Room_B_BRIR_16k\", filename));
[revC, fs] = audioread(strcat("stimuli\unvocoded\UniS_Room_C_BRIR_16k\", filename));
[revD, fs] = audioread(strcat("stimuli\unvocoded\UniS_Room_D_BRIR_16k\", filename));

raw_data = {an, revA, revB, revC, revD};

%% calculate start- and end-point on the audio
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

%% apply orignial vocoder
[an_bi, ~, ~] = BiCI_EE(raw_data{1}, fs, 'Med-El', 'none', 0);
[revA_bi, ~, ~] = BiCI_EE(raw_data{2}, fs, 'Med-El', 'none', 0);
[revB_bi, ~, ~] = BiCI_EE(raw_data{3}, fs, 'Med-El', 'none', 0);
[revC_bi, ~, ~] = BiCI_EE(raw_data{4}, fs, 'Med-El', 'none', 0);
[revD_bi, ~, ~] = BiCI_EE(raw_data{5}, fs, 'Med-El', 'none', 0);

%% apply fullband enhancement vocoder
[revA_bi_full3, ~, ~] = BiCI_EE(raw_data{2}, fs, 'Med-El', 'fullband', 3);
[revB_bi_full3, ~, ~] = BiCI_EE(raw_data{3}, fs, 'Med-El', 'fullband', 3);
[revC_bi_full3, ~, ~] = BiCI_EE(raw_data{4}, fs, 'Med-El', 'fullband', 3);
[revD_bi_full3, ~, ~] = BiCI_EE(raw_data{5}, fs, 'Med-El', 'fullband', 3);

[revA_bi_full4, ~, ~] = BiCI_EE(raw_data{2}, fs, 'Med-El', 'fullband', 4);
[revB_bi_full4, ~, ~] = BiCI_EE(raw_data{3}, fs, 'Med-El', 'fullband', 4);
[revC_bi_full4, ~, ~] = BiCI_EE(raw_data{4}, fs, 'Med-El', 'fullband', 4);
[revD_bi_full4, ~, ~] = BiCI_EE(raw_data{5}, fs, 'Med-El', 'fullband', 4);

[revA_bi_full5, ~, ~] = BiCI_EE(raw_data{2}, fs, 'Med-El', 'fullband', 5);
[revB_bi_full5, ~, ~] = BiCI_EE(raw_data{3}, fs, 'Med-El', 'fullband', 5);
[revC_bi_full5, ~, ~] = BiCI_EE(raw_data{4}, fs, 'Med-El', 'fullband', 5);
[revD_bi_full5, ~, ~] = BiCI_EE(raw_data{5}, fs, 'Med-El', 'fullband', 5);

[revA_bi_full6, ~, ~] = BiCI_EE(raw_data{2}, fs, 'Med-El', 'fullband', 6);
[revB_bi_full6, ~, ~] = BiCI_EE(raw_data{3}, fs, 'Med-El', 'fullband', 6);
[revC_bi_full6, ~, ~] = BiCI_EE(raw_data{4}, fs, 'Med-El', 'fullband', 6);
[revD_bi_full6, ~, ~] = BiCI_EE(raw_data{5}, fs, 'Med-El', 'fullband', 6);


%% apply subband signal enhancement
[revA_bi_FS3, ~, ~] = BiCI_EE(raw_data{2}, fs, 'Med-El', 'FS', 3);
[revB_bi_FS3, ~, ~] = BiCI_EE(raw_data{3}, fs, 'Med-El', 'FS', 3);
[revC_bi_FS3, ~, ~] = BiCI_EE(raw_data{4}, fs, 'Med-El', 'FS', 3);
[revD_bi_FS3, ~, ~] = BiCI_EE(raw_data{5}, fs, 'Med-El', 'FS', 3);

[revA_bi_FS4, ~, ~] = BiCI_EE(raw_data{2}, fs, 'Med-El', 'FS', 4);
[revB_bi_FS4, ~, ~] = BiCI_EE(raw_data{3}, fs, 'Med-El', 'FS', 4);
[revC_bi_FS4, ~, ~] = BiCI_EE(raw_data{4}, fs, 'Med-El', 'FS', 4);
[revD_bi_FS4, ~, ~] = BiCI_EE(raw_data{5}, fs, 'Med-El', 'FS', 4);

[revA_bi_FS5, ~, ~] = BiCI_EE(raw_data{2}, fs, 'Med-El', 'FS', 5);
[revB_bi_FS5, ~, ~] = BiCI_EE(raw_data{3}, fs, 'Med-El', 'FS', 5);
[revC_bi_FS5, ~, ~] = BiCI_EE(raw_data{4}, fs, 'Med-El', 'FS', 5);
[revD_bi_FS5, ~, ~] = BiCI_EE(raw_data{5}, fs, 'Med-El', 'FS', 5);

[revA_bi_FS6, ~, ~] = BiCI_EE(raw_data{2}, fs, 'Med-El', 'FS', 6);
[revB_bi_FS6, ~, ~] = BiCI_EE(raw_data{3}, fs, 'Med-El', 'FS', 6);
[revC_bi_FS6, ~, ~] = BiCI_EE(raw_data{4}, fs, 'Med-El', 'FS', 6);
[revD_bi_FS6, ~, ~] = BiCI_EE(raw_data{5}, fs, 'Med-El', 'FS', 6);

%% apply subband envelope enhancement
[revA_bi_EE3, ~, ~] = BiCI_EE(raw_data{2}, fs, 'Med-El', 'EE', 3);
[revB_bi_EE3, ~, ~] = BiCI_EE(raw_data{3}, fs, 'Med-El', 'EE', 3);
[revC_bi_EE3, ~, ~] = BiCI_EE(raw_data{4}, fs, 'Med-El', 'EE', 3);
[revD_bi_EE3, ~, ~] = BiCI_EE(raw_data{5}, fs, 'Med-El', 'EE', 3);

[revA_bi_EE4, ~, ~] = BiCI_EE(raw_data{2}, fs, 'Med-El', 'EE', 4);
[revB_bi_EE4, ~, ~] = BiCI_EE(raw_data{3}, fs, 'Med-El', 'EE', 4);
[revC_bi_EE4, ~, ~] = BiCI_EE(raw_data{4}, fs, 'Med-El', 'EE', 4);
[revD_bi_EE4, ~, ~] = BiCI_EE(raw_data{5}, fs, 'Med-El', 'EE', 4);

[revA_bi_EE5, ~, ~] = BiCI_EE(raw_data{2}, fs, 'Med-El', 'EE', 5);
[revB_bi_EE5, ~, ~] = BiCI_EE(raw_data{3}, fs, 'Med-El', 'EE', 5);
[revC_bi_EE5, ~, ~] = BiCI_EE(raw_data{4}, fs, 'Med-El', 'EE', 5);
[revD_bi_EE5, ~, ~] = BiCI_EE(raw_data{5}, fs, 'Med-El', 'EE', 5);

[revA_bi_EE6, ~, ~] = BiCI_EE(raw_data{2}, fs, 'Med-El', 'EE', 6);
[revB_bi_EE6, ~, ~] = BiCI_EE(raw_data{3}, fs, 'Med-El', 'EE', 6);
[revC_bi_EE6, ~, ~] = BiCI_EE(raw_data{4}, fs, 'Med-El', 'EE', 6);
[revD_bi_EE6, ~, ~] = BiCI_EE(raw_data{5}, fs, 'Med-El', 'EE', 6);

%% calculate metric intelligibility
% baseline
an_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), an_bi(:,1), an_bi(:,2), fs);
revA_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revA_bi(:,1), revA_bi(:,2), fs);
revB_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revB_bi(:,1), revB_bi(:,2), fs);
revC_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revC_bi(:,1), revC_bi(:,2), fs);
revD_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revD_bi(:,1), revD_bi(:,2), fs);

% full enhancement
revA_full3_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revA_bi_full3(:,1), revA_bi_full3(:,2), fs);
revB_full3_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revB_bi_full3(:,1), revB_bi_full3(:,2), fs);
revC_full3_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revC_bi_full3(:,1), revC_bi_full3(:,2), fs);
revD_full3_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revD_bi_full3(:,1), revD_bi_full3(:,2), fs);

revA_full4_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revA_bi_full4(:,1), revA_bi_full4(:,2), fs);
revB_full4_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revB_bi_full4(:,1), revB_bi_full4(:,2), fs);
revC_full4_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revC_bi_full4(:,1), revC_bi_full4(:,2), fs);
revD_full4_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revD_bi_full4(:,1), revD_bi_full4(:,2), fs);

revA_full5_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revA_bi_full5(:,1), revA_bi_full5(:,2), fs);
revB_full5_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revB_bi_full5(:,1), revB_bi_full5(:,2), fs);
revC_full5_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revC_bi_full5(:,1), revC_bi_full5(:,2), fs);
revD_full5_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revD_bi_full5(:,1), revD_bi_full5(:,2), fs);

revA_full6_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revA_bi_full6(:,1), revA_bi_full6(:,2), fs);
revB_full6_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revB_bi_full6(:,1), revB_bi_full6(:,2), fs);
revC_full6_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revC_bi_full6(:,1), revC_bi_full6(:,2), fs);
revD_full6_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revD_bi_full6(:,1), revD_bi_full6(:,2), fs);

% fine subband enhancement
revA_FS3_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revA_bi_FS3(:,1), revA_bi_FS3(:,2), fs);
revB_FS3_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revB_bi_FS3(:,1), revB_bi_FS3(:,2), fs);
revC_FS3_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revC_bi_FS3(:,1), revC_bi_FS3(:,2), fs);
revD_FS3_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revD_bi_FS3(:,1), revD_bi_FS3(:,2), fs);

revA_FS4_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revA_bi_FS4(:,1), revA_bi_FS4(:,2), fs);
revB_FS4_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revB_bi_FS4(:,1), revB_bi_FS4(:,2), fs);
revC_FS4_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revC_bi_FS4(:,1), revC_bi_FS4(:,2), fs);
revD_FS4_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revD_bi_FS4(:,1), revD_bi_FS4(:,2), fs);

revA_FS5_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revA_bi_FS5(:,1), revA_bi_FS5(:,2), fs);
revB_FS5_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revB_bi_FS5(:,1), revB_bi_FS5(:,2), fs);
revC_FS5_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revC_bi_FS5(:,1), revC_bi_FS5(:,2), fs);
revD_FS5_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revD_bi_FS5(:,1), revD_bi_FS5(:,2), fs);

revA_FS6_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revA_bi_FS6(:,1), revA_bi_FS6(:,2), fs);
revB_FS6_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revB_bi_FS6(:,1), revB_bi_FS6(:,2), fs);
revC_FS6_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revC_bi_FS6(:,1), revC_bi_FS6(:,2), fs);
revD_FS6_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revD_bi_FS6(:,1), revD_bi_FS6(:,2), fs);

% envelope enhancement
revA_EE3_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revA_bi_EE3(:,1), revA_bi_EE3(:,2), fs);
revB_EE3_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revB_bi_EE3(:,1), revB_bi_EE3(:,2), fs);
revC_EE3_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revC_bi_EE3(:,1), revC_bi_EE3(:,2), fs);
revD_EE3_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revD_bi_EE3(:,1), revD_bi_EE3(:,2), fs);

revA_EE4_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revA_bi_EE4(:,1), revA_bi_EE4(:,2), fs);
revB_EE4_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revB_bi_EE4(:,1), revB_bi_EE4(:,2), fs);
revC_EE4_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revC_bi_EE4(:,1), revC_bi_EE4(:,2), fs);
revD_EE4_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revD_bi_EE4(:,1), revD_bi_EE4(:,2), fs);

revA_EE5_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revA_bi_EE5(:,1), revA_bi_EE5(:,2), fs);
revB_EE5_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revB_bi_EE5(:,1), revB_bi_EE5(:,2), fs);
revC_EE5_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revC_bi_EE5(:,1), revC_bi_EE5(:,2), fs);
revD_EE5_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revD_bi_EE5(:,1), revD_bi_EE5(:,2), fs);

revA_EE6_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revA_bi_EE6(:,1), revA_bi_EE6(:,2), fs);
revB_EE6_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revB_bi_EE6(:,1), revB_bi_EE6(:,2), fs);
revC_EE6_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revC_bi_EE6(:,1), revC_bi_EE6(:,2), fs);
revD_EE6_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), revD_bi_EE6(:,1), revD_bi_EE6(:,2), fs);

%% calculate interaural level difference
% baseline
an_ild = 20*log10(rms(an_bi(:,2))/rms(an_bi(:,1)));
revA_ild = 20*log10(rms(revA_bi(:,2))/rms(revA_bi(:,1)));
revB_ild = 20*log10(rms(revB_bi(:,2))/rms(revB_bi(:,1)));
revC_ild = 20*log10(rms(revC_bi(:,2))/rms(revC_bi(:,1)));
revD_ild = 20*log10(rms(revD_bi(:,2))/rms(revD_bi(:,1)));

% full enhancement
revA_full3_ild = 20*log10(rms(revA_bi_full3(:,2))/rms(revA_bi_full3(:,1)));
revB_full3_ild = 20*log10(rms(revB_bi_full3(:,2))/rms(revB_bi_full3(:,1)));
revC_full3_ild = 20*log10(rms(revC_bi_full3(:,2))/rms(revC_bi_full3(:,1)));
revD_full3_ild = 20*log10(rms(revD_bi_full3(:,2))/rms(revD_bi_full3(:,1)));

revA_full4_ild = 20*log10(rms(revA_bi_full4(:,2))/rms(revA_bi_full4(:,1)));
revB_full4_ild = 20*log10(rms(revB_bi_full4(:,2))/rms(revB_bi_full4(:,1)));
revC_full4_ild = 20*log10(rms(revC_bi_full4(:,2))/rms(revC_bi_full4(:,1)));
revD_full4_ild = 20*log10(rms(revD_bi_full4(:,2))/rms(revD_bi_full4(:,1)));

revA_full5_ild = 20*log10(rms(revA_bi_full5(:,2))/rms(revA_bi_full5(:,1)));
revB_full5_ild = 20*log10(rms(revB_bi_full5(:,2))/rms(revB_bi_full5(:,1)));
revC_full5_ild = 20*log10(rms(revC_bi_full5(:,2))/rms(revC_bi_full5(:,1)));
revD_full5_ild = 20*log10(rms(revD_bi_full5(:,2))/rms(revD_bi_full5(:,1)));

revA_full6_ild = 20*log10(rms(revA_bi_full6(:,2))/rms(revA_bi_full6(:,1)));
revB_full6_ild = 20*log10(rms(revB_bi_full6(:,2))/rms(revB_bi_full6(:,1)));
revC_full6_ild = 20*log10(rms(revC_bi_full6(:,2))/rms(revC_bi_full6(:,1)));
revD_full6_ild = 20*log10(rms(revD_bi_full6(:,2))/rms(revD_bi_full6(:,1)));

% subband enhancement
revA_EE3_ild = 20*log10(rms(revA_bi_EE3(:,2))/rms(revA_bi_EE3(:,1)));
revB_EE3_ild = 20*log10(rms(revB_bi_EE3(:,2))/rms(revB_bi_EE3(:,1)));
revC_EE3_ild = 20*log10(rms(revC_bi_EE3(:,2))/rms(revC_bi_EE3(:,1)));
revD_EE3_ild = 20*log10(rms(revD_bi_EE3(:,2))/rms(revD_bi_EE3(:,1)));

revA_EE4_ild = 20*log10(rms(revA_bi_EE4(:,2))/rms(revA_bi_EE4(:,1)));
revB_EE4_ild = 20*log10(rms(revB_bi_EE4(:,2))/rms(revB_bi_EE4(:,1)));
revC_EE4_ild = 20*log10(rms(revC_bi_EE4(:,2))/rms(revC_bi_EE4(:,1)));
revD_EE4_ild = 20*log10(rms(revD_bi_EE4(:,2))/rms(revD_bi_EE4(:,1)));

revA_EE5_ild = 20*log10(rms(revA_bi_EE5(:,2))/rms(revA_bi_EE5(:,1)));
revB_EE5_ild = 20*log10(rms(revB_bi_EE5(:,2))/rms(revB_bi_EE5(:,1)));
revC_EE5_ild = 20*log10(rms(revC_bi_EE5(:,2))/rms(revC_bi_EE5(:,1)));
revD_EE5_ild = 20*log10(rms(revD_bi_EE5(:,2))/rms(revD_bi_EE5(:,1)));

revA_EE6_ild = 20*log10(rms(revA_bi_EE6(:,2))/rms(revA_bi_EE6(:,1)));
revB_EE6_ild = 20*log10(rms(revB_bi_EE6(:,2))/rms(revB_bi_EE6(:,1)));
revC_EE6_ild = 20*log10(rms(revC_bi_EE6(:,2))/rms(revC_bi_EE6(:,1)));
revD_EE6_ild = 20*log10(rms(revD_bi_EE6(:,2))/rms(revD_bi_EE6(:,1)));

% fine enhancement
revA_FS3_ild = 20*log10(rms(revA_bi_FS3(:,2))/rms(revA_bi_FS3(:,1)));
revB_FS3_ild = 20*log10(rms(revB_bi_FS3(:,2))/rms(revB_bi_FS3(:,1)));
revC_FS3_ild = 20*log10(rms(revC_bi_FS3(:,2))/rms(revC_bi_FS3(:,1)));
revD_FS3_ild = 20*log10(rms(revD_bi_FS3(:,2))/rms(revD_bi_FS3(:,1)));

revA_FS4_ild = 20*log10(rms(revA_bi_FS4(:,2))/rms(revA_bi_FS4(:,1)));
revB_FS4_ild = 20*log10(rms(revB_bi_FS4(:,2))/rms(revB_bi_FS4(:,1)));
revC_FS4_ild = 20*log10(rms(revC_bi_FS4(:,2))/rms(revC_bi_FS4(:,1)));
revD_FS4_ild = 20*log10(rms(revD_bi_FS4(:,2))/rms(revD_bi_FS4(:,1)));

revA_FS5_ild = 20*log10(rms(revA_bi_FS5(:,2))/rms(revA_bi_FS5(:,1)));
revB_FS5_ild = 20*log10(rms(revB_bi_FS5(:,2))/rms(revB_bi_FS5(:,1)));
revC_FS5_ild = 20*log10(rms(revC_bi_FS5(:,2))/rms(revC_bi_FS5(:,1)));
revD_FS5_ild = 20*log10(rms(revD_bi_FS5(:,2))/rms(revD_bi_FS5(:,1)));

revA_FS6_ild = 20*log10(rms(revA_bi_FS6(:,2))/rms(revA_bi_FS6(:,1)));
revB_FS6_ild = 20*log10(rms(revB_bi_FS6(:,2))/rms(revB_bi_FS6(:,1)));
revC_FS6_ild = 20*log10(rms(revC_bi_FS6(:,2))/rms(revC_bi_FS6(:,1)));
revD_FS6_ild = 20*log10(rms(revD_bi_FS6(:,2))/rms(revD_bi_FS6(:,1)));


%% calcualte interaural time difference
%  baseline
an_itd = estimate_ITD_Broadband(an_bi, fs);
revA_itd = estimate_ITD_Broadband(revA_bi, fs);
revB_itd = estimate_ITD_Broadband(revB_bi, fs);
revC_itd = estimate_ITD_Broadband(revC_bi, fs);
revD_itd = estimate_ITD_Broadband(revD_bi, fs);

% full enhancement
revA_full3_itd = estimate_ITD_Broadband(revA_bi_full3, fs);
revB_full3_itd = estimate_ITD_Broadband(revB_bi_full3, fs);
revC_full3_itd = estimate_ITD_Broadband(revC_bi_full3, fs);
revD_full3_itd = estimate_ITD_Broadband(revD_bi_full3, fs);

revA_full4_itd = estimate_ITD_Broadband(revA_bi_full4, fs);
revB_full4_itd = estimate_ITD_Broadband(revB_bi_full4, fs);
revC_full4_itd = estimate_ITD_Broadband(revC_bi_full4, fs);
revD_full4_itd = estimate_ITD_Broadband(revD_bi_full4, fs);

revA_full5_itd = estimate_ITD_Broadband(revA_bi_full5, fs);
revB_full5_itd = estimate_ITD_Broadband(revB_bi_full5, fs);
revC_full5_itd = estimate_ITD_Broadband(revC_bi_full5, fs);
revD_full5_itd = estimate_ITD_Broadband(revD_bi_full5, fs);

revA_full6_itd = estimate_ITD_Broadband(revA_bi_full6, fs);
revB_full6_itd = estimate_ITD_Broadband(revB_bi_full6, fs);
revC_full6_itd = estimate_ITD_Broadband(revC_bi_full6, fs);
revD_full6_itd = estimate_ITD_Broadband(revD_bi_full6, fs);

% subband enhancement
revA_FS3_itd = estimate_ITD_Broadband(revA_bi_FS3, fs);
revB_FS3_itd = estimate_ITD_Broadband(revB_bi_FS3, fs);
revC_FS3_itd = estimate_ITD_Broadband(revC_bi_FS3, fs);
revD_FS3_itd = estimate_ITD_Broadband(revD_bi_FS3, fs);

revA_FS4_itd = estimate_ITD_Broadband(revA_bi_FS4, fs);
revB_FS4_itd = estimate_ITD_Broadband(revB_bi_FS4, fs);
revC_FS4_itd = estimate_ITD_Broadband(revC_bi_FS4, fs);
revD_FS4_itd = estimate_ITD_Broadband(revD_bi_FS4, fs);

revA_FS5_itd = estimate_ITD_Broadband(revA_bi_FS5, fs);
revB_FS5_itd = estimate_ITD_Broadband(revB_bi_FS5, fs);
revC_FS5_itd = estimate_ITD_Broadband(revC_bi_FS5, fs);
revD_FS5_itd = estimate_ITD_Broadband(revD_bi_FS5, fs);

revA_FS6_itd = estimate_ITD_Broadband(revA_bi_FS6, fs);
revB_FS6_itd = estimate_ITD_Broadband(revB_bi_FS6, fs);
revC_FS6_itd = estimate_ITD_Broadband(revC_bi_FS6, fs);
revD_FS6_itd = estimate_ITD_Broadband(revD_bi_FS6, fs);

% envelope enhancement
revA_EE3_itd = estimate_ITD_Broadband(revA_bi_EE3, fs);
revB_EE3_itd = estimate_ITD_Broadband(revB_bi_EE3, fs);
revC_EE3_itd = estimate_ITD_Broadband(revC_bi_EE3, fs);
revD_EE3_itd = estimate_ITD_Broadband(revD_bi_EE3, fs);

revA_EE4_itd = estimate_ITD_Broadband(revA_bi_EE4, fs);
revB_EE4_itd = estimate_ITD_Broadband(revB_bi_EE4, fs);
revC_EE4_itd = estimate_ITD_Broadband(revC_bi_EE4, fs);
revD_EE4_itd = estimate_ITD_Broadband(revD_bi_EE4, fs);

revA_EE5_itd = estimate_ITD_Broadband(revA_bi_EE5, fs);
revB_EE5_itd = estimate_ITD_Broadband(revB_bi_EE5, fs);
revC_EE5_itd = estimate_ITD_Broadband(revC_bi_EE5, fs);
revD_EE5_itd = estimate_ITD_Broadband(revD_bi_EE5, fs);

revA_EE6_itd = estimate_ITD_Broadband(revA_bi_EE6, fs);
revB_EE6_itd = estimate_ITD_Broadband(revB_bi_EE6, fs);
revC_EE6_itd = estimate_ITD_Broadband(revC_bi_EE6, fs);
revD_EE6_itd = estimate_ITD_Broadband(revD_bi_EE6, fs);

