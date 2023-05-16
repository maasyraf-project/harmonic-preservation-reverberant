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

% spatial cues