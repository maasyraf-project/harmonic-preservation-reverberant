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

%% vocoder process
% parameter
par.voc_sampling_frequency_hz = 16e3;
par.gamma_order_stimulation = 3;
par.center_frequencies_hz_stimulation = [120 235 384 579 836 1175 1624 2222 3019 4084 5507 7410];
par.bandwidth_factor = [1 1 1 1 1 1 1 1 1 1 1 1].*3;
par.weights = [0.98 0.98 0.98 0.68 0.68 0.45 0.45 0.2 0.2 0.15 0.15 0.15]';

% resampling
if par.voc_sampling_frequency_hz ~= fs
    for i = 1:length(pe_data)
        pe_data{i} = resample(pe_data{i}, par.voc_sampling_frequency_hz, fs);
    end
end

% pre-emphasis filter
pe_data = cell(1, length(raw_data));
for i = 1:length(raw_data)
    w     = 2*1200/fs;
    [b,a] = butter(1,w,'high');
    l_sig = filter(b, a, raw_data{i}(:,1));
    r_sig = filter(b, a, raw_data{i}(:,2));
    pe_data{i} = [l_sig, r_sig]; 
end

% Create analyse -Filterbank
analyzer_stim = Gfb_Analyzer_new(par.voc_sampling_frequency_hz,...
    par.gamma_order_stimulation, ...
    par.center_frequencies_hz_stimulation,...
    par.bandwidth_factor);

an_data = cell(1, length(pe_data));
for i = 1:length(pe_data)
    [l_sig, ~] = Gfb_Analyzer_process(analyzer_stim, pe_data{i}(:,1)');
    [r_sig, ~] = Gfb_Analyzer_process(analyzer_stim, pe_data{i}(:,2)');
    an_data{i} = {l_sig, r_sig}; 
end

env_data = cell(1, length(an_data));
fine_data = cell(1, length(an_data));
for i = 1:length(an_data)
    % extract envelope
    el_sig = abs(an_data{i}{1});
    er_sig = abs(an_data{i}{2});

    if size(el_sig,2) ~= size(er_sig,2)
        error("Length of left and right analyzed signal is not same!")
    end

    weights = repmat(par.weights,1,size(el_sig,2));
    el_sig = sqrt(weights.*(real(an_data{i}{1}).^2+imag(an_data{i}{1}).^2));
    er_sig = sqrt(weights.*(real(an_data{i}{2}).^2+imag(an_data{i}{2}).^2));
    
    env_data{i} = {el_sig, er_sig}; 
    
    % extract fine structure
    fl_sig = real(an_data{i}{1});
    fr_sig = real(an_data{i}{2});

    fine_data{i} = {fl_sig, fr_sig};

end

env_lp_data = cell(1, length(env_data));
for i = 1:length(env_data)
    l_sub_buf = nan(size(env_data{i}{1}));
    r_sub_buf = nan(size(env_data{i}{2}));

    for j = 1:size(env_data{1}{1},1)
        [b, a] = butter(1, (200./(fs/2)));
        l_sig = filter(b, a, env_data{i}{1}(j,:));
        r_sig = filter(b, a, env_data{i}{2}(j,:));

        l_sub_buf(j,:) = l_sig;
        r_sub_buf(j,:) = r_sig;
    end
    env_lp_data{i} = {l_sub_buf, r_sub_buf};

end

env_lp_kurt_data = cell(1, length(env_lp_data));
for i = 1:length(env_lp_data)
    l_sub_buf = nan(size(env_lp_data{i}{1}));
    r_sub_buf = nan(size(env_lp_data{i}{2}));

    for j = 1:size(env_data{1}{1},1)
        lseg = 400;
        theta = 1e-6;
        l_sig = kurt(env_data{i}{1}(j,:), lseg, theta);
        r_sig = kurt(env_data{i}{2}(j,:), lseg, theta);

        l_sub_buf(j,:) = l_sig;
        r_sub_buf(j,:) = r_sig;
    end
    env_lp_kurt_data{i} = {l_sub_buf, r_sub_buf};

end

