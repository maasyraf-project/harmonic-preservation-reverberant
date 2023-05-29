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

%% apply orignial vocoder
[an_bi, ~, ~] = BiCI_Sim_EAF(an, fs, 'Med-El', 'none');
[lrev_bi, ~, ~] = BiCI_Sim_EAF(lrev, fs, 'Med-El', 'none');
[hrev_bi, ~, ~] = BiCI_Sim_EAF(hrev, fs, 'Med-El', 'none');

%% apply fullband speech enhancement
[an_bi_full, ~, ~] = BiCI_Sim_EAF(an, fs, 'Med-El', 'fullband');
[lrev_bi_full, ~, ~] = BiCI_Sim_EAF(lrev, fs, 'Med-El', 'fullband');
[hrev_bi_full, ~, ~] = BiCI_Sim_EAF(hrev, fs, 'Med-El', 'fullband');

%% apply subband speech enhancement
[an_bi_subenv, ~, ~] = BiCI_Sim_EAF(an, fs, 'Med-El', 'subband-envelope');
[lrev_bi_subenv, ~, ~] = BiCI_Sim_EAF(lrev, fs, 'Med-El', 'subband-envelope');
[hrev_bi_subenv, ~, ~] = BiCI_Sim_EAF(hrev, fs, 'Med-El', 'subband-envelope');

[an_bi_subsig, ~, ~] = BiCI_Sim_EAF(an, fs, 'Med-El', 'subband-signal');
[lrev_bi_subsig, ~, ~] = BiCI_Sim_EAF(lrev, fs, 'Med-El', 'subband-signal');
[hrev_bi_subsig, ~, ~] = BiCI_Sim_EAF(hrev, fs, 'Med-El', 'subband-signal');

%% calculate metric
% intelligibility
an_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), an_bi(:,1), an_bi(:,2), fs);
lrev_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), lrev_bi(:,1), lrev_bi(:,2), fs);
hrev_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), hrev_bi(:,1), hrev_bi(:,2), fs);

lrev_full_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), lrev_bi_full(:,1), lrev_bi_full(:,2), fs);
hrev_full_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), hrev_bi_full(:,1), hrev_bi_full(:,2), fs);

lrev_subenv_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), lrev_bi_subenv(:,1), lrev_bi_subenv(:,2), fs);
hrev_subenv_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), hrev_bi_subenv(:,1), hrev_bi_subenv(:,2), fs);

lrev_subsig_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), lrev_bi_subsig(:,1), lrev_bi_subsig(:,2), fs);
hrev_subsig_mbstoi = mbstoi(an_bi(:,1), an_bi(:,2), hrev_bi_subsig(:,1), hrev_bi_subsig(:,2), fs);

% spatial cues (ILD)
an_ild = 20*log10(rms(an_bi(:,2))/rms(an_bi(:,1)));
lrev_ild = 20*log10(rms(lrev_bi(:,2))/rms(lrev_bi(:,1)));
hrev_ild = 20*log10(rms(hrev_bi(:,2))/rms(hrev_bi(:,1)));

lrev_full_ild = 20*log10(rms(lrev_bi_full(:,2))/rms(lrev_bi_full(:,1)));
hrev_full_ild = 20*log10(rms(hrev_bi_full(:,2))/rms(hrev_bi_full(:,1)));

lrev_subenv_ild = 20*log10(rms(lrev_bi_subenv(:,2))/rms(lrev_bi_subenv(:,1)));
hrev_subenv_ild = 20*log10(rms(hrev_bi_subenv(:,2))/rms(hrev_bi_subenv(:,1)));

lrev_subsig_ild = 20*log10(rms(lrev_bi_subsig(:,2))/rms(lrev_bi_subsig(:,1)));
hrev_subsig_ild = 20*log10(rms(hrev_bi_subsig(:,2))/rms(hrev_bi_subsig(:,1)));

% spatial cues (ITD)
an_itd = estimate_ITD_Broadband(an_bi, fs);
lrev_itd = estimate_ITD_Broadband(lrev_bi, fs);
hrev_itd = estimate_ITD_Broadband(hrev_bi, fs);

lrev_full_itd = estimate_ITD_Broadband(lrev_bi_full, fs);
hrev_full_itd = estimate_ITD_Broadband(hrev_bi_full, fs);

lrev_subenv_itd = estimate_ITD_Broadband(lrev_bi_subenv, fs);
hrev_subenv_itd = estimate_ITD_Broadband(hrev_bi_subenv, fs);

lrev_subsig_itd = estimate_ITD_Broadband(lrev_bi_subsig, fs);
hrev_subsig_itd = estimate_ITD_Broadband(hrev_bi_subsig, fs);