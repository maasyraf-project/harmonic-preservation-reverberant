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