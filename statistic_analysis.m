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
%% vocoder preparation
% parameter
par.voc_sampling_frequency_hz = 48e3;
par.gamma_order_stimulation = 3;
par.gamma_order_auralisation = 3;
par.center_frequencies_hz_stimulation = [120 235 384 579 836 1175 1624 2222 3019 4084 5507 7410];
par.center_frequencies_hz_auralisation = [357 548 689 968 1483 2228 3319 4670 6630 9758 12530 15374];
par.bandwidth_factor = [1 1 1 1 1 1 1 1 1 1 1 1].*3;
par.weights = [0.98 0.98 0.98 0.68 0.68 0.45 0.45 0.2 0.2 0.15 0.15 0.15]';

%% pre-emphasis filter
pe_data = cell(1, length(raw_data));
for i = 1:length(raw_data)
    w     = 2*1200/fs;
    [b,a] = butter(1,w,'high');
    l_sig = filter(b, a, raw_data{i}(:,1));
    r_sig = filter(b, a, raw_data{i}(:,2));
    pe_data{i} = [l_sig, r_sig]; 
end

%% resampling 
if par.voc_sampling_frequency_hz ~= fs
    for i = 1:length(pe_data)
        pe_data{i} = resample(pe_data{i}, par.voc_sampling_frequency_hz, fs);
    end
    fs = par.voc_sampling_frequency_hz;
    max_len_data = max([length(pe_data{1}), length(pe_data{2}), length(pe_data{3})]);
    t = 0:(1/fs):(max_len_data-1)/fs;
end

%% apply full band LP analysis
% split audio into several frames
for i = 1:length(pe_data)
    lenFrames = 0.032*fs;
    lenOverlap = 0.5*lenFrames;
    nFrames = floor((length(pe_data{i}(:,1)) - lenOverlap)/lenOverlap) + 1;

    estl = zeros(1, length(pe_data{i}(:,1)));
    estr = zeros(1, length(pe_data{i}(:,1)));
    yresl = zeros(1, length(pe_data{i}(:,1)));
    yresr = zeros(1, length(pe_data{i}(:,1)));

    l = pe_data{i}(:,1);
    r = pe_data{i}(:,2);

    for n = 1:nFrames
        % define start and end index
        idxStart = 1 + (n-1) * (lenFrames - lenOverlap);
        idxEnd = idxStart + lenFrames -1;

        if idxEnd > length(pe_data{i}(:,1))
            idxEnd = length(pe_data{i}(:,1));
        end

        l_frame = l(idxStart:idxEnd);
        r_frame = r(idxStart:idxEnd);

        % apply LPC on left channel
        la = lpc(l_frame, 12);
        estl_frame = filter([0 -la(2:end)], 1, l_frame);
        resl_frame = l_frame - estl_frame;

        estl(idxStart:idxEnd) = estl_frame;
        yresl(idxStart:idxEnd) = resl_frame;

        % apply LPC on right channel
        ra = lpc(r_frame, 12);
        estr_frame = filter([0 -ra(2:end)], 1, r_frame);
        resr_frame = r_frame - estr_frame;

        estr(idxStart:idxEnd) = estr_frame;
        yresr(idxStart:idxEnd) = resr_frame;

    end

    pe_data{i} = [yresl', yresr']; 
end

%% decompose into subband
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
rms_env_data = cell(1, length(an_data));
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
    
    rmsl_sig = rms2(el_sig,2);
    rmsr_sig = rms2(er_sig,2);

    env_data{i} = {el_sig, er_sig}; 
    rms_env_data{i} = {rmsl_sig, rmsr_sig};
    
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

%% fine structure subband analysis
fine_cm_data = cell(1, length(fine_data));
for i = 1:length(fine_data)
    % split audio into several frames
    lenFrames = 0.032*fs;
    lenOverlap = 0.5*lenFrames;
    nFrames = floor((size(fine_data{i}{1},2) - lenOverlap)/lenOverlap) + 1;
    cmOrd = 1:6;

    l_sub_buf = cell(size(fine_data{i}{1}, 1), 1);
    r_sub_buf = cell(size(fine_data{i}{2}, 1), 1);
    
    for j = 1:size(fine_data{1}{1},1)
        l = fine_data{i}{1}(j,:);
        r = fine_data{i}{2}(j,:);

        cml = zeros(length(cmOrd), nFrames);
        cmr = zeros(length(cmOrd), nFrames);

        for n = 1:nFrames
            % define start and end index
            idxStart = 1 + (n-1) * (lenFrames - lenOverlap);
            idxEnd = idxStart + lenFrames -1;

            if idxEnd > length(l)
                idxEnd = length(l);
            end

            % conduct central moment analysis
            l_frame = l(idxStart:idxEnd);
            r_frame = r(idxStart:idxEnd);

            if or(length(l_frame) < lenFrames, length(r_frame) < lenFrames)
                l_frame = [l_frame zeros(1, lenFrames - length(l_frame))];
                r_frame = [r_frame zeros(1, lenFrames - length(r_frame))];
            end

            for k = 1:length(cmOrd)

                lm = moment(l_frame, cmOrd(k))/std(l_frame)^cmOrd(k);
                rm = moment(r_frame, cmOrd(k))/std(r_frame)^cmOrd(k);

                cml(k,n) = lm;
                cmr(k,n) = rm;

            end

        end

        % check whether NaN data appear during processing
        if or(anynan(cml) == 1, anynan(cmr) == 1)
            warning(strcat("NaN value appear on the channel of ", string(j)))
        end

        l_sub_buf{j} = cml;
        r_sub_buf{j} = cmr;

    end

    fine_cm_data{i} = {l_sub_buf, r_sub_buf};

end


%% envelope structure subband analysis
env_lp_cm_data = cell(1, length(env_lp_data));
for i = 1:length(env_lp_data)
    % split audio into several frames
    lenFrames = 0.032*fs;
    lenOverlap = 0.5*lenFrames;
    nFrames = floor((size(env_lp_data{i}{1},2) - lenOverlap)/lenOverlap) + 1;
    cmOrd = 1:6;

    l_sub_buf = cell(size(env_lp_data{i}{1}, 1), 1);
    r_sub_buf = cell(size(env_lp_data{i}{2}, 1), 1);
    
    for j = 1:size(env_lp_data{1}{1},1)
        l = env_lp_data{i}{1}(j,:);
        r = env_lp_data{i}{2}(j,:);

        cml = zeros(length(cmOrd), nFrames);
        cmr = zeros(length(cmOrd), nFrames);

        for n = 1:nFrames
            % define start and end index
            idxStart = 1 + (n-1) * (lenFrames - lenOverlap);
            idxEnd = idxStart + lenFrames -1;

            if idxEnd > length(l)
                idxEnd = length(l);
            end

            % conduct central moment analysis
            l_frame = l(idxStart:idxEnd);
            r_frame = r(idxStart:idxEnd);

            if or(length(l_frame) < lenFrames, length(r_frame) < lenFrames)
                l_frame = [l_frame zeros(1, lenFrames - length(l_frame))];
                r_frame = [r_frame zeros(1, lenFrames - length(r_frame))];
            end

            for k = 1:length(cmOrd)

                lm = moment(l_frame, cmOrd(k))/std(l_frame)^cmOrd(k);
                rm = moment(r_frame, cmOrd(k))/std(r_frame)^cmOrd(k);

                cml(k,n) = lm;
                cmr(k,n) = rm;

            end

        end

        % check whether NaN data appear during processing
        if or(anynan(cml) == 1, anynan(cmr) == 1)
            warning(strcat("NaN value appear on the channel of ", string(j)))
        end

        l_sub_buf{j} = cml;
        r_sub_buf{j} = cmr;

    end

    env_lp_cm_data{i} = {l_sub_buf, r_sub_buf};

end

%% apply LP-analysis on fine structure subband 
fine_lp_cm_data = cell(1, length(fine_data));
fine_lp_data = cell(1, length(fine_data));

for i = 1:length(fine_data)
    % split audio into several frames
    lenFrames = 0.032*fs;
    lenOverlap = 0.5*lenFrames;
    nFrames = floor((size(fine_data{i}{1},2) - lenOverlap)/lenOverlap) + 1;
    cmOrd = 1:6;

    l_sub_buf = cell(size(fine_data{i}{1}, 1), 1);
    r_sub_buf = cell(size(fine_data{i}{2}, 1), 1);

    l_res_buf = nan(size(fine_data{i}{1}));
    r_res_buf = nan(size(fine_data{i}{2}));
    
    for j = 1:size(fine_data{1}{1},1)
        l = fine_data{i}{1}(j,:);
        r = fine_data{i}{2}(j,:);

        cml = zeros(length(cmOrd), nFrames);
        cmr = zeros(length(cmOrd), nFrames);

        estl = zeros(size(l));
        estr = zeros(size(r));

        for n = 1:nFrames
            % define start and end index
            idxStart = 1 + (n-1) * (lenFrames - lenOverlap);
            idxEnd = idxStart + lenFrames -1;

            if idxEnd > length(l)
                idxEnd = length(l);
            end
            
            l_frame = l(idxStart:idxEnd);
            r_frame = r(idxStart:idxEnd);

            % apply LPC
            la = lpc(l_frame, 12);
            estl_frame = filter([0 -la(2:end)], 1, l_frame);% .* hann(length(idxStart:idxEnd))');
            l_frame = l_frame - estl_frame;

            ra = lpc(r_frame, 12);
            estr_frame = filter([0 -ra(2:end)], 1, r_frame);% .* hann(length(idxStart:idxEnd))');
            r_frame = r_frame - estr_frame;

            % conduct central moment analysis
            for k = 1:length(cmOrd)

                lm = moment(l_frame, cmOrd(k)) / (std(l_frame)^cmOrd(k));
                rm = moment(r_frame, cmOrd(k)) / (std(r_frame)^cmOrd(k));

                cml(k,n) = lm;
                cmr(k,n) = rm;

            end
            
            estl(idxStart:idxEnd) = l_frame(1:length(idxStart:idxEnd));
            estr(idxStart:idxEnd) = r_frame(1:length(idxStart:idxEnd));

        end

        % check whether NaN data appear during processing
        if or(anynan(cml) == 1, anynan(cmr) == 1)
            warning(strcat("NaN value appear on the channel of ", string(j)))
        end

        l_sub_buf{j} = cml;
        r_sub_buf{j} = cmr;

        l_res_buf(j,:) = estl;
        r_res_buf(j,:) = estr;

    end

    fine_lp_cm_data{i} = {l_sub_buf, r_sub_buf};
    fine_lp_data{i} = {l_res_buf, r_res_buf};

end

%% apply LP-analysis on envelope structure subband 
env_lp_lp_cm_data = cell(1, length(env_lp_data));
env_lp_lp_data = cell(1, length(env_lp_data));

for i = 1:length(env_lp_data)
    % split audio into several frames
    lenFrames = 0.032*fs;
    lenOverlap = 0.5*lenFrames;
    nFrames = floor((size(env_lp_data{i}{1},2) - lenOverlap)/lenOverlap) + 1;
    cmOrd = 1:6;

    l_sub_buf = cell(size(env_lp_data{i}{1}, 1), 1);
    r_sub_buf = cell(size(env_lp_data{i}{2}, 1), 1);

    l_res_buf = nan(size(env_lp_data{i}{1}));
    r_res_buf = nan(size(env_lp_data{i}{2}));
    
    for j = 1:size(env_lp_data{1}{1},1)
        l = env_lp_data{i}{1}(j,:);
        r = env_lp_data{i}{2}(j,:);

        cml = zeros(length(cmOrd), nFrames);
        cmr = zeros(length(cmOrd), nFrames);

        estl = zeros(size(l));
        estr = zeros(size(r));

        for n = 1:nFrames
            % define start and end index
            idxStart = 1 + (n-1) * (lenFrames - lenOverlap);
            idxEnd = idxStart + lenFrames -1;

            if idxEnd > length(l)
                idxEnd = length(l);
            end

            l_frame = l(idxStart:idxEnd);
            r_frame = r(idxStart:idxEnd);

            % apply LPC
            la = lpc(l_frame, 3);
            l_frame = filter(la, 1, l_frame);

            ra = lpc(r_frame, 3);
            r_frame = filter(ra, 1, r_frame);

            % conduct central moment analysis
            for k = 1:length(cmOrd)

                lm = moment(l_frame, cmOrd(k))/std(l_frame)^cmOrd(k);
                rm = moment(r_frame, cmOrd(k))/std(r_frame)^cmOrd(k);

                cml(k,n) = lm;
                cmr(k,n) = rm;

            end
            
            estl(idxStart:idxEnd) = l_frame(1:length(idxStart:idxEnd));
            estr(idxStart:idxEnd) = r_frame(1:length(idxStart:idxEnd));

        end

        % check whether NaN data appear during processing
        if or(anynan(cml) == 1, anynan(cmr) == 1)
            warning(strcat("NaN value appear on the channel of ", string(j)))
        end

        l_sub_buf{j} = cml;
        r_sub_buf{j} = cmr;

        l_res_buf(j,:) = estl;
        r_res_buf(j,:) = estr;

    end

    env_lp_lp_cm_data{i} = {l_sub_buf, r_sub_buf};
    env_lp_lp_data{i} = {l_res_buf, r_res_buf};

end

%% plot the fine data
data_label = ["clean speech", "low reverberant speech"];
%data_label = ["clean speech", "high reverberant speech"];
figure
for i = 1:2%length(env_lp_data)
    for j = 1:6%size(env_data{1}{1},1)
        subplot(6,1,j)
        if i == 1
            plot(t(1:length(fine_data{i}{1}(j,:))), fine_data{i}{1}(j,:), "Color", [0 0 0 1])
        else
            plot(t(1:length(fine_data{i}{1}(j,:))), fine_data{i}{1}(j,:), "Color", [0.5 0.5 0.5 0.8])
        end
        legend(data_label)
        hold on
        ylabel("Central Moment Value")
        if j == 6
            xlabel("Time (s)")
        end
    end
end

figure
for i = 1:2%length(env_lp_data)
    for j = 1:6%size(env_data{1}{1},1)
        subplot(6,1,j)
        if i == 1
            plot(t(1:length(fine_data{i}{1}(j+6,:))), fine_data{i}{1}(j+6,:), "Color", [0 0 0 1])
        else
            plot(t(1:length(fine_data{i}{1}(j+6,:))), fine_data{i}{1}(j+6,:), "Color", [0.5 0.5 0.5 0.8])
        end
        legend(data_label)
        hold on
        ylabel("Central Moment Value")
        if j == 6
            xlabel("Time (s)")
        end
    end
end

%% plot the kurtosis of fine data
data_label = ["clean speech", "low reverberant speech"];
%data_label = ["clean speech", "high reverberant speech"];
figure
chan = 1;
cm = 4;
t_frames = linspace(0, length(fine_data{1}{1})/fs, nFrames);
for i = 1:2%length(env_lp_data)
    for j = 1:6%size(env_data{1}{1},1)
        subplot(6,1,j)
        if i == 1
            plot(t(1:length(fine_cm_data{i}{chan}{j}(cm,:))), fine_cm_data{i}{chan}{j}(cm,:), "Color", [0 0 0 1])
        else
            plot(t(1:length(fine_cm_data{i}{chan}{j}(cm,:))), fine_cm_data{i}{chan}{j}(cm,:), "Color", [0.5 0.5 0.5 0.8])
        end
        legend(data_label)
        hold on
        ylabel("Central Moment Value")
        if j == 6
            xlabel("Time (s)")
        end
    end
end

figure
for i = 1:2%length(env_lp_data)
    for j = 1:6%size(env_data{1}{1},1)
        subplot(6,1,j)
        if i == 1
            plot(t(1:length(fine_cm_data{i}{chan}{j+6}(cm,:))), fine_cm_data{i}{chan}{j+6}(cm,:), "Color", [0 0 0 1])
        else
            plot(t(1:length(fine_cm_data{i}{chan}{j+6}(cm,:))), fine_cm_data{i}{chan}{j+6}(cm,:), "Color", [0.5 0.5 0.5 0.8])
        end
        legend(data_label)
        hold on
        ylabel("Central Moment Value")
        if j == 6
            xlabel("Time (s)")
        end
    end
end

%% plot the LP of fine data
data_label = ["clean speech", "low reverberant speech"];
%data_label = ["clean speech", "high reverberant speech"];
figure
for i = 1:2%length(env_lp_data)
    for j = 1:6%size(env_data{1}{1},1)
        subplot(6,1,j)
        if i == 1
            plot(t(1:length(fine_lp_data{i}{1}(j,:))), fine_lp_data{i}{1}(j,:), "Color", [0 0 0 1])
        else
            plot(t(1:length(fine_lp_data{i}{1}(j,:))), fine_lp_data{i}{1}(j,:), "Color", [0.5 0.5 0.5 0.8])
        end
        legend(data_label)
        hold on
        ylabel("Central Moment Value")
        if j == 6
            xlabel("Time (s)")
        end
    end
end

figure
for i = 1:2%length(env_lp_data)
    for j = 1:6%size(env_data{1}{1},1)
        subplot(6,1,j)
        if i == 1
            plot(t(1:length(fine_lp_data{i}{1}(j+6,:))), fine_lp_data{i}{1}(j+6,:), "Color", [0 0 0 1])
        else
            plot(t(1:length(fine_lp_data{i}{1}(j+6,:))), fine_lp_data{i}{1}(j+6,:), "Color", [0.5 0.5 0.5 0.8])
        end
        legend(data_label)
        hold on
        ylabel("Central Moment Value")
        if j == 6
            xlabel("Time (s)")
        end
    end
end

%% plot the kurtosis of LP of fine data
data_label = ["clean speech", "low reverberant speech"];
%data_label = ["clean speech", "high reverberant speech"];
chan = 1;
cm = 4;
t_frames = linspace(0, length(fine_lp_data{1}{1})/fs, nFrames);
figure
for i = 1:2%length(env_lp_data)
    for j = 1:6%size(env_data{1}{1},1)
        subplot(6,1,j)
        if i == 1
            plot(t_frames(1:length(fine_lp_cm_data{i}{chan}{j}(cm,:))), fine_lp_cm_data{i}{chan}{j}(cm,:), "Color", [0 0 0 1])
        else
            plot(t_frames(1:length(fine_lp_cm_data{i}{chan}{j}(cm,:))), fine_lp_cm_data{i}{chan}{j}(cm,:), "Color", [0.5 0.5 0.5 0.8])
        end
        legend(data_label)
        hold on
        ylabel("Central Moment Value")
        if j == 6
            xlabel("Time (s)")
        end
    end
end

figure
for i = 1:2%length(env_lp_data)
    for j = 1:6%size(env_data{1}{1},1)
        subplot(6,1,j)
        if i == 1
            plot(t_frames(1:length(fine_lp_cm_data{i}{chan}{j}(cm,:))), fine_lp_cm_data{i}{chan}{j+6}(cm,:), "Color", [0 0 0 1])
        else
            plot(t_frames(1:length(fine_lp_cm_data{i}{chan}{j}(cm,:))), fine_lp_cm_data{i}{chan}{j+6}(cm,:), "Color", [0.5 0.5 0.5 0.8])
        end
        legend(data_label)
        hold on
        ylabel("Central Moment Value")
        if j == 6
            xlabel("Time (s)")
        end
    end
end

%% plot the env data
data_label = ["clean speech", "low reverberant speech"];
%data_label = ["clean speech", "high reverberant speech"];
figure
for i = 1:2%length(env_lp_data)
    for j = 1:6%size(env_data{1}{1},1)
        subplot(6,1,j)
        if i == 1
            plot(t(1:length(env_lp_data{i}{1}(j,:))), env_lp_data{i}{1}(j,:), "Color", [0 0 0 1])
        else
            plot(t(1:length(env_lp_data{i}{1}(j,:))), env_lp_data{i}{1}(j,:), "Color", [0.5 0.5 0.5 0.8])
        end
        legend(data_label)
        hold on
        ylabel("Central Moment Value")
        if j == 6
            xlabel("Time (s)")
        end
    end
end

figure
for i = 1:2%length(env_lp_data)
    for j = 1:6%size(env_data{1}{1},1)
        subplot(6,1,j)
        if i == 1
            plot(t(1:length(env_lp_data{i}{1}(j+6,:))), env_lp_data{i}{1}(j+6,:), "Color", [0 0 0 1])
        else
            plot(t(1:length(env_lp_data{i}{1}(j+6,:))), env_lp_data{i}{1}(j+6,:), "Color", [0.5 0.5 0.5 0.8])
        end
        legend(data_label)
        hold on
        ylabel("Central Moment Value")
        if j == 6
            xlabel("Time (s)")
        end
    end
end

%% plot the kurtosis of env data
data_label = ["clean speech", "low reverberant speech"];
%data_label = ["clean speech", "high reverberant speech"];
figure
for i = 1:2%length(env_lp_data)
    for j = 1:6%size(env_data{1}{1},1)
        subplot(6,1,j)
        if i == 1
            plot(t(1:length(env_lp_cm_data{i}{chan}{j}(cm,:))), env_lp_cm_data{i}{chan}{j}(cm,:), "Color", [0 0 0 1])
        else
            plot(t(1:length(env_lp_cm_data{i}{chan}{j}(cm,:))), env_lp_cm_data{i}{chan}{j}(cm,:), "Color", [0.5 0.5 0.5 0.8])
        end
        legend(data_label)
        hold on
        ylabel("Central Moment Value")
        if j == 6
            xlabel("Time (s)")
        end
    end
end

figure
for i = 1:2%length(env_lp_data)
    for j = 1:6%size(env_data{1}{1},1)
        subplot(6,1,j)
        if i == 1
            plot(t(1:length(env_lp_cm_data{i}{chan}{j+6}(cm,:))), env_lp_cm_data{i}{chan}{j+6}(cm,:), "Color", [0 0 0 1])
        else
            plot(t(1:length(env_lp_cm_data{i}{chan}{j+6}(cm,:))), env_lp_cm_data{i}{chan}{j+6}(cm,:), "Color", [0.5 0.5 0.5 0.8])
        end
        legend(data_label)
        hold on
        ylabel("Central Moment Value")
        if j == 6
            xlabel("Time (s)")
        end
    end
end

%% plot the LP of env data
data_label = ["clean speech", "low reverberant speech"];
%data_label = ["clean speech", "high reverberant speech"];
figure
for i = 1:2%length(env_lp_data)
    for j = 1:6%size(env_data{1}{1},1)
        subplot(6,1,j)
        if i == 1
            plot(t(1:length(env_lp_lp_data{i}{1}(j,:))), env_lp_lp_data{i}{1}(j,:), "Color", [0 0 0 1])
        else
            plot(t(1:length(env_lp_lp_data{i}{1}(j,:))), env_lp_lp_data{i}{1}(j,:), "Color", [0.5 0.5 0.5 0.8])
        end
        legend(data_label)
        hold on
        ylabel("Central Moment Value")
        if j == 6
            xlabel("Time (s)")
        end
    end
end

figure
for i = 1:2%length(env_lp_data)
    for j = 1:6%size(env_data{1}{1},1)
        subplot(6,1,j)
        if i == 1
            plot(t(1:length(env_lp_lp_data{i}{1}(j+6,:))), env_lp_lp_data{i}{1}(j+6,:), "Color", [0 0 0 1])
        else
            plot(t(1:length(env_lp_lp_data{i}{1}(j+6,:))), env_lp_lp_data{i}{1}(j+6,:), "Color", [0.5 0.5 0.5 0.8])
        end
        legend(data_label)
        hold on
        ylabel("Central Moment Value")
        if j == 6
            xlabel("Time (s)")
        end
    end
end


%% plot the kurtosis of LP of fine data
data_label = ["clean speech", "low reverberant speech"];
%data_label = ["clean speech", "high reverberant speech"];
chan = 1;
cm = 5;
t_frames = linspace(0, length(env_lp_data{1}{1})/fs, nFrames);
figure
for i = 1:2:3%length(env_lp_data)
    for j = 1:6%size(env_data{1}{1},1)
        subplot(6,1,j)
        if i == 1
            plot(t_frames(1:length(env_lp_lp_cm_data{i}{chan}{j}(cm,:))), env_lp_lp_cm_data{i}{chan}{j}(cm,:), "Color", [0 0 0 1])
        else
            plot(t_frames(1:length(env_lp_lp_cm_data{i}{chan}{j}(cm,:))), env_lp_lp_cm_data{i}{chan}{j}(cm,:), "Color", [0.5 0.5 0.5 0.8])
        end
        legend(data_label)
        hold on
        ylabel("Central Moment Value")
        if j == 6
            xlabel("Time (s)")
        end
    end
end

figure
for i = 1:2:3%length(env_lp_data)
    for j = 1:6%size(env_data{1}{1},1)
        subplot(6,1,j)
        if i == 1
            plot(t_frames(1:length(env_lp_lp_cm_data{i}{chan}{j}(cm,:))), env_lp_lp_cm_data{i}{chan}{j+6}(cm,:), "Color", [0 0 0 1])
        else
            plot(t_frames(1:length(env_lp_lp_cm_data{i}{chan}{j}(cm,:))), env_lp_lp_cm_data{i}{chan}{j+6}(cm,:), "Color", [0.5 0.5 0.5 0.8])
        end
        legend(data_label)
        hold on
        ylabel("Central Moment Value")
        if j == 6
            xlabel("Time (s)")
        end
    end
end

%% customized plot
% fine low-rev
figure
subplot(3,2,1)
plot(t(1:length(fine_data{1}{1}(1,:))), fine_data{1}{1}(1,:), "Color", [0 0 0 1])
hold on 
plot(t(1:length(fine_data{2}{1}(1,:))), fine_data{3}{1}(1,:), "Color", [0.5 0.5 0.5 0.4])
ylabel("Amplitude")
title("Fine Structure of Subband (Left Channel)")

subplot(3,2,2)
plot(t(1:length(fine_data{1}{2}(1,:))), fine_data{1}{2}(1,:), "Color", [0 0 0 1])
hold on
plot(t(1:length(fine_data{2}{2}(1,:))), fine_data{3}{2}(1,:), "Color", [0.5 0.5 0.5 0.4])
ylabel("Amplitude")
title("Fine Structure of Subband (Right Channel)")

subplot(3,2,3)
plot(t(1:length(fine_data{1}{1}(2,:))), fine_data{1}{1}(2,:), "Color", [0 0 0 1])
hold on
plot(t(1:length(fine_data{2}{1}(2,:))), fine_data{3}{1}(2,:), "Color", [0.5 0.5 0.5 0.4])
ylabel("Amplitude")

subplot(3,2,4)
plot(t(1:length(fine_data{1}{2}(2,:))), fine_data{1}{2}(2,:), "Color", [0 0 0 1])
hold on
plot(t(1:length(fine_data{2}{2}(2,:))), fine_data{3}{2}(2,:), "Color", [0.5 0.5 0.5 0.4])
ylabel("Amplitude")

subplot(3,2,5)
plot(t(1:length(fine_data{1}{1}(3,:))), fine_data{1}{1}(3,:), "Color", [0 0 0 1])
hold on
plot(t(1:length(fine_data{2}{1}(3,:))), fine_data{3}{1}(3,:), "Color", [0.5 0.5 0.5 0.4])
xlabel("Time (s)")
ylabel("Amplitude")

subplot(3,2,6)
plot(t(1:length(fine_data{1}{2}(3,:))), fine_data{i}{2}(3,:), "Color", [0 0 0 1])
hold on
plot(t(1:length(fine_data{2}{2}(3,:))), fine_data{3}{2}(3,:), "Color", [0.5 0.5 0.5 0.4])
ylabel("Amplitude")
xlabel("Time (s)")

%% customized plot
% env_lp high-rev
figure
subplot(3,2,1)
plot(t(1:length(env_lp_data{1}{1}(1,:))), env_lp_data{1}{1}(1,:), "Color", [0 0 0 1])
hold on 
plot(t(1:length(env_lp_data{2}{1}(1,:))), env_lp_data{3}{1}(1,:), "Color", [0.5 0.5 0.5 0.4])
ylabel("Amplitude")
%ylim([0 0.01])
title("Envelope Structure of Subband (Left Channel)")

subplot(3,2,2)
plot(t(1:length(env_lp_data{1}{2}(1,:))), env_lp_data{1}{2}(1,:), "Color", [0 0 0 1])
hold on
plot(t(1:length(env_lp_data{2}{2}(1,:))), env_lp_data{3}{2}(1,:), "Color", [0.5 0.5 0.5 0.4])
ylabel("Amplitude")
%ylim([0 0.01])
title("Envelope Structure of Subband (Right Channel)")

subplot(3,2,3)
plot(t(1:length(env_lp_data{1}{1}(2,:))), env_lp_data{1}{1}(2,:), "Color", [0 0 0 1])
hold on
plot(t(1:length(env_lp_data{2}{1}(2,:))), env_lp_data{3}{1}(2,:), "Color", [0.5 0.5 0.5 0.4])
ylabel("Amplitude")
%ylim([0 0.02])

subplot(3,2,4)
plot(t(1:length(env_lp_data{1}{2}(2,:))), env_lp_data{1}{2}(2,:), "Color", [0 0 0 1])
hold on
plot(t(1:length(env_lp_data{2}{2}(2,:))), env_lp_data{3}{2}(2,:), "Color", [0.5 0.5 0.5 0.4])
ylabel("Amplitude")
%ylim([0 0.02])

subplot(3,2,5)
plot(t(1:length(env_lp_data{1}{1}(3,:))), env_lp_data{1}{1}(3,:), "Color", [0 0 0 1])
hold on
plot(t(1:length(env_lp_data{2}{1}(3,:))), env_lp_data{3}{1}(3,:), "Color", [0.5 0.5 0.5 0.4])
xlabel("Time (s)")
ylabel("Amplitude")
%ylim([0 0.02])

subplot(3,2,6)
plot(t(1:length(env_lp_data{1}{2}(3,:))), env_lp_data{i}{2}(3,:), "Color", [0 0 0 1])
hold on
plot(t(1:length(env_lp_data{2}{2}(3,:))), env_lp_data{3}{2}(3,:), "Color", [0.5 0.5 0.5 0.4])
ylabel("Amplitude")
xlabel("Time (s)")
%ylim([0 0.02])

%% customized plot
% env_lp high-rev
cm = 4;
t_frames = linspace(0, length(env_lp_data{1}{1})/fs, nFrames);
figure
subplot(3,2,1)
plot(t_frames(1:length(env_lp_cm_data{1}{1}{1}(cm,:))), env_lp_cm_data{1}{1}{1}(cm,:), "Color", 'black')
hold on 
plot(t_frames(1:length(env_lp_cm_data{2}{1}{1}(cm,:))), env_lp_cm_data{2}{1}{1}(cm,:), "Color", 'red')
hold on
plot(t_frames(1:length(env_lp_cm_data{3}{1}{1}(cm,:))), env_lp_cm_data{3}{1}{1}(cm,:), "Color", 'blue')
ylabel("Value")
%ylim([0 0.01])
title("High Order Statistic (HOS) of Subband (Left Channel)")

subplot(3,2,2)
plot(t_frames(1:length(env_lp_cm_data{1}{2}{1}(cm,:))), env_lp_cm_data{1}{2}{1}(cm,:), "Color", 'black')
hold on 
plot(t_frames(1:length(env_lp_cm_data{2}{2}{1}(cm,:))), env_lp_cm_data{2}{2}{1}(cm,:), "Color", 'red')
hold on
plot(t_frames(1:length(env_lp_cm_data{3}{2}{1}(cm,:))), env_lp_cm_data{3}{2}{1}(cm,:), "Color", 'blue')
ylabel("Value")
%ylim([0 0.01])
title("High Order Statistic (HOS) Structure of Subband (Right Channel)")

subplot(3,2,3)
plot(t_frames(1:length(env_lp_cm_data{1}{1}{2}(cm,:))), env_lp_cm_data{1}{1}{2}(cm,:), "Color", 'black')
hold on 
plot(t_frames(1:length(env_lp_cm_data{2}{1}{2}(cm,:))), env_lp_cm_data{2}{1}{2}(cm,:), "Color", 'red')
hold on
plot(t_frames(1:length(env_lp_cm_data{3}{1}{2}(cm,:))), env_lp_cm_data{3}{1}{2}(cm,:), "Color", 'blue')
ylabel("Value")
%ylim([0 0.02])

subplot(3,2,4)
plot(t_frames(1:length(env_lp_cm_data{1}{2}{2}(cm,:))), env_lp_cm_data{1}{2}{2}(cm,:), "Color", 'black')
hold on 
plot(t_frames(1:length(env_lp_cm_data{2}{2}{2}(cm,:))), env_lp_cm_data{2}{2}{2}(cm,:), "Color", 'red')
hold on
plot(t_frames(1:length(env_lp_cm_data{3}{2}{2}(cm,:))), env_lp_cm_data{3}{2}{2}(cm,:), "Color", 'blue')
ylabel("Value")
%ylim([0 0.02])

subplot(3,2,5)
plot(t_frames(1:length(env_lp_cm_data{1}{1}{3}(cm,:))), env_lp_cm_data{1}{1}{3}(cm,:), "Color", 'black')
hold on 
plot(t_frames(1:length(env_lp_cm_data{2}{1}{3}(cm,:))), env_lp_cm_data{2}{1}{3}(cm,:), "Color", 'red')
hold on
plot(t_frames(1:length(env_lp_cm_data{3}{1}{3}(cm,:))), env_lp_cm_data{3}{1}{3}(cm,:), "Color", 'blue')
xlabel("Time (s)")
ylabel("Value")
%ylim([0 0.02])

subplot(3,2,6)
plot(t_frames(1:length(env_lp_cm_data{1}{2}{3}(cm,:))), env_lp_cm_data{1}{2}{3}(cm,:), "Color", 'black')
hold on 
plot(t_frames(1:length(env_lp_cm_data{2}{2}{3}(cm,:))), env_lp_cm_data{2}{2}{3}(cm,:), "Color", 'red')
hold on
plot(t_frames(1:length(env_lp_cm_data{3}{2}{3}(cm,:))), env_lp_cm_data{3}{2}{3}(cm,:), "Color", 'blue')
ylabel("Value")
xlabel("Time (s)")
%ylim([0 0.02])



%% customized plot
% fine high-rev
cm = 6;
t_frames = linspace(0, length(fine_data{1}{1})/fs, nFrames);
figure
subplot(3,2,1)
plot(t_frames(1:length(fine_cm_data{1}{1}{1}(cm,:))), fine_cm_data{1}{1}{1}(cm,:), "Color", 'black')
hold on 
plot(t_frames(1:length(fine_cm_data{2}{1}{1}(cm,:))), fine_cm_data{2}{1}{1}(cm,:), "Color", 'red')
hold on
plot(t_frames(1:length(fine_cm_data{3}{1}{1}(cm,:))), fine_cm_data{3}{1}{1}(cm,:), "Color", 'blue')
ylabel("Value")
%ylim([0 0.01])
title("High Order Statistic (HOS) of Subband (Left Channel)")

subplot(3,2,2)
plot(t_frames(1:length(fine_cm_data{1}{2}{1}(cm,:))), fine_cm_data{1}{2}{1}(cm,:), "Color", 'black')
hold on 
plot(t_frames(1:length(fine_cm_data{2}{2}{1}(cm,:))), fine_cm_data{2}{2}{1}(cm,:), "Color", 'red')
hold on
plot(t_frames(1:length(fine_cm_data{3}{2}{1}(cm,:))), fine_cm_data{3}{2}{1}(cm,:), "Color", 'blue')
ylabel("Value")
%ylim([0 0.01])
title("High Order Statistic (HOS) Structure of Subband (Right Channel)")

subplot(3,2,3)
plot(t_frames(1:length(fine_cm_data{1}{1}{2}(cm,:))), fine_cm_data{1}{1}{2}(cm,:), "Color", 'black')
hold on 
plot(t_frames(1:length(fine_cm_data{2}{1}{2}(cm,:))), fine_cm_data{2}{1}{2}(cm,:), "Color", 'red')
hold on
plot(t_frames(1:length(fine_cm_data{3}{1}{2}(cm,:))), fine_cm_data{3}{1}{2}(cm,:), "Color", 'blue')
ylabel("Value")
%ylim([0 0.02])

subplot(3,2,4)
plot(t_frames(1:length(fine_cm_data{1}{2}{2}(cm,:))), fine_cm_data{1}{2}{2}(cm,:), "Color", 'black')
hold on 
plot(t_frames(1:length(fine_cm_data{2}{2}{2}(cm,:))), fine_cm_data{2}{2}{2}(cm,:), "Color", 'red')
hold on
plot(t_frames(1:length(fine_cm_data{3}{2}{2}(cm,:))), fine_cm_data{3}{2}{2}(cm,:), "Color", 'blue')
ylabel("Value")
%ylim([0 0.02])

subplot(3,2,5)
plot(t_frames(1:length(fine_cm_data{1}{1}{3}(cm,:))), fine_cm_data{1}{1}{3}(cm,:), "Color", 'black')
hold on 
plot(t_frames(1:length(fine_cm_data{2}{1}{3}(cm,:))), fine_cm_data{2}{1}{3}(cm,:), "Color", 'red')
hold on
plot(t_frames(1:length(fine_cm_data{3}{1}{3}(cm,:))), fine_cm_data{3}{1}{3}(cm,:), "Color", 'blue')
xlabel("Time (s)")
ylabel("Value")
%ylim([0 0.02])

subplot(3,2,6)
plot(t_frames(1:length(fine_cm_data{1}{2}{3}(cm,:))), fine_cm_data{1}{2}{3}(cm,:), "Color", 'black')
hold on 
plot(t_frames(1:length(fine_cm_data{2}{2}{3}(cm,:))), fine_cm_data{2}{2}{3}(cm,:), "Color", 'red')
hold on
plot(t_frames(1:length(fine_cm_data{3}{2}{3}(cm,:))), fine_cm_data{3}{2}{3}(cm,:), "Color", 'blue')
ylabel("Value")
xlabel("Time (s)")
%ylim([0 0.02])





