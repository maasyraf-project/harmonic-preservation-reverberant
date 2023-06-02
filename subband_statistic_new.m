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
    l_sub_buf = nan(size(fine_data{i}{1}));
    r_sub_buf = nan(size(fine_data{i}{2}));

    % split audio into several frames
    lenFrames = 0.03*fs;
    lenOverlap = 0.5*lenFrames;
    nFrames = floor(size(fine_data{i}{1},2)/(lenFrames-lenOverlap)); % the denumerator is false

    l_vad_buf = nan(size(fine_data{i}{1},1),nFrames);
    r_vad_buf = nan(size(fine_data{i}{2},1),nFrames);

    for j = 1:size(fine_data{1}{1},1)
        l = fine_data{i}{1}(j,:);
        r = fine_data{i}{2}(j,:);

        idxStart = 1;
        idxEnd = lenFrames;

        estl = zeros(size(l));
        estr = zeros(size(r));

        for n = 1:nFrames
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

            % calculate central moment analysis
            %l_frame = l_frame * hamming(lenFrames);
            %r_frame = r_frame * hamming(lenFrames);
            
            lm = moment(l_frame, 4)/std(l_frame)^4;
            rm = moment(r_frame, 4)/std(r_frame)^4;
            
            l_frame(idxStart:idxEnd) = lm;
            r_frame(idxStart:idxEnd) = rm;

            estl(idxStart:idxEnd) = estl(idxStart:idxEnd) + l_frame(idxStart:idxEnd);
            estr(idxStart:idxEnd) = estr(idxStart:idxEnd) + r_frame(idxStart:idxEnd);
            
            % update index for selecting next frame
            idxStart = idxStart + lenOverlap;
            idxEnd = idxEnd + lenOverlap;
            
        end
        
        % recap all central moment value
        l_sig = estl;
        r_sig = estr;

        % check whether NaN data appear during processing
        if anynan(l_sig) == 1
            warning(strcat("NaN value appear ", string(j)))
        end

        l_sub_buf(j,:) = l_sig;
        r_sub_buf(j,:) = r_sig;

    end

    fine_cm_data{i} = {l_sub_buf, r_sub_buf};

end


%% envelope subband analysis
env_lp_cm_data = cell(1, length(env_lp_data));
for i = 1:length(env_lp_data)
    l_sub_buf = nan(size(env_lp_data{i}{1}));
    r_sub_buf = nan(size(env_lp_data{i}{2}));

    % split audio into several frames
    lenFrames = 0.03*fs;
    lenOverlap = 0.5*lenFrames;
    nFrames = floor(size(fine_data{i}{1},2)/(lenFrames-lenOverlap));

    for j = 1:size(env_lp_data{1}{1},1)
        l = env_lp_data{i}{1}(j,:);
        r = env_lp_data{i}{2}(j,:);

        idxStart = 1;
        idxEnd = lenFrames;

%         estl = kurt(l, lenFrames, 1e-6);
%         estr = kurt(r, lenFrames, 1e-6);


        estl = zeros(size(l));
        estr = zeros(size(r));

        for n = 1:nFrames
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

            % calculate central moment analysis
            %l_frame = l_frame * hamming(lenFrames);
            %r_frame = r_frame * hamming(lenFrames);
            
            lm = moment(l_frame, 4)/std(l_frame)^4;
            rm = moment(r_frame, 4)/std(r_frame)^4;
            
            l_frame(idxStart:idxEnd) = lm;
            r_frame(idxStart:idxEnd) = rm;

            estl(idxStart:idxEnd) = estl(idxStart:idxEnd) + l_frame(idxStart:idxEnd);
            estr(idxStart:idxEnd) = estr(idxStart:idxEnd) + r_frame(idxStart:idxEnd);
            
            % update index for selecting next frame
            idxStart = idxStart + lenOverlap;
            idxEnd = idxEnd + lenOverlap;
            
        end
        
        % recap all central moment value
        l_sig = estl;
        r_sig = estr;

        % check whether NaN data appear during processing
        if anynan(l_sig) == 1
            warning(strcat("NaN value appear ", string(j)))
        end

        l_sub_buf(j,:) = l_sig;
        r_sub_buf(j,:) = r_sig;

    end

    env_lp_cm_data{i} = {l_sub_buf, r_sub_buf};

end

%% conduct LP-analysis on fine data
fine_lp_data = cell(1, length(fine_data));
fine_lp_cm_data = cell(1, length(fine_lp_data));
for i = 1:length(fine_data)
    l_sub_buf = nan(size(fine_data{i}{1}));
    r_sub_buf = nan(size(fine_data{i}{2}));

    l_cm_buf = nan(size(fine_data{i}{1}));
    r_cm_buf = nan(size(fine_data{i}{2}));

    % split audio into several frames
    lenFrames = 0.03*fs;
    lenOverlap = 0.5*lenFrames;
    nFrames = floor(size(fine_data{i}{1},2)/(lenFrames-lenOverlap)); 

    for j = 1:size(fine_data{1}{1},1)
        l = fine_data{i}{1}(j,:);
        r = fine_data{i}{2}(j,:);

        idxStart = 1;
        idxEnd = lenFrames;

        estl = zeros(size(l));
        estr = zeros(size(r));

        cml = zeros(size(l));
        cmr = zeros(size(r));

        for n = 1:nFrames
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


            % calculate central moment analysis
            l_frame = hann(lenFrames)' .* l_frame ;
            r_frame = hann(lenFrames)' .* r_frame ;
            
            la = lpc(l_frame, 10);
            estl_frame = filter([0 -la(2:end)], 1, l_frame);

            ra = lpc(r_frame, 10);
            estr_frame = filter([0 -ra(2:end)], 1, r_frame);
            
            % calculate residue
            if idxEnd == length(l)
                estl(idxStart:idxEnd) = estl(idxStart:idxEnd) + estl_frame(1:numel(idxStart:idxEnd)) - l_frame(1:numel(idxStart:idxEnd));
                estr(idxStart:idxEnd) = estr(idxStart:idxEnd) + estr_frame(1:numel(idxStart:idxEnd)) - r_frame(1:numel(idxStart:idxEnd));
            else
                estl(idxStart:idxEnd) = estl(idxStart:idxEnd) + estl_frame - l_frame;
                estr(idxStart:idxEnd) = estr(idxStart:idxEnd) + estr_frame - r_frame;
            end

            % calculate center moment of frame
            lm = moment(estl_frame, 4)/std(estl_frame)^4;
            rm = moment(estr_frame, 4)/std(estr_frame)^4;

            cml_frame = lm;
            cmr_frame = rm;

            cml(idxStart:idxEnd) = cml(idxStart:idxEnd) + cml_frame;
            cmr(idxStart:idxEnd) = cmr(idxStart:idxEnd) + cmr_frame;
            
            
            % update index for selecting next frame
            idxStart = idxStart + lenOverlap;
            idxEnd = idxEnd + lenOverlap;
            
        end
        
        % obtain residue of signal
        l_sig = estl;
        r_sig = estr;

        cml_sig = cml;
        cmr_sig = cmr;

        % check whether NaN data appear during processing
        if anynan(l_sig) == 1
            warning(strcat("NaN value appear ", string(j)))
        end

        l_sub_buf(j,:) = l_sig;
        r_sub_buf(j,:) = r_sig;

        l_cm_buf(j,:) = cml_sig;
        r_cm_buf(j,:) = cmr_sig;

    end

    fine_lp_data{i} = {l_sub_buf, r_sub_buf};
    fine_lp_cm_data{i} = {l_cm_buf, r_cm_buf};

end

%% conduct LP-analysis on env data
env_lp_lp_data = cell(1, length(env_lp_data));
env_lp_lp_cm_data = cell(1, length(env_lp_lp_data));
for i = 1:length(env_lp_data)
    l_sub_buf = nan(size(env_lp_data{i}{1}));
    r_sub_buf = nan(size(env_lp_data{i}{2}));
    
    l_cm_buf = nan(size(fine_data{i}{1}));
    r_cm_buf = nan(size(fine_data{i}{2}));

    % split audio into several frames
    lenFrames = 0.03*fs;
    lenOverlap = 0.5*lenFrames;
    nFrames = floor(size(fine_data{i}{1},2)/(lenFrames-lenOverlap)); 

    for j = 1:size(fine_data{1}{1},1)
        l = fine_data{i}{1}(j,:);
        r = fine_data{i}{2}(j,:);

        idxStart = 1;
        idxEnd = lenFrames;

        estl = zeros(size(l));
        estr = zeros(size(r));

        cml = zeros(size(l));
        cmr = zeros(size(r));

        for n = 1:nFrames
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


            % calculate central moment analysis
            l_frame = hann(lenFrames)' .* l_frame ;
            r_frame = hann(lenFrames)' .* r_frame ;
            
            la = lpc(l_frame, 10);
            estl_frame = filter([0 -la(2:end)], 1, l_frame);

            ra = lpc(r_frame, 10);
            estr_frame = filter([0 -ra(2:end)], 1, r_frame);
            
            % calculate residue
            if idxEnd == length(l)
                estl(idxStart:idxEnd) = estl(idxStart:idxEnd) + estl_frame(1:numel(idxStart:idxEnd)) - l_frame(1:numel(idxStart:idxEnd));
                estr(idxStart:idxEnd) = estr(idxStart:idxEnd) + estr_frame(1:numel(idxStart:idxEnd)) - r_frame(1:numel(idxStart:idxEnd));
            else
                estl(idxStart:idxEnd) = estl(idxStart:idxEnd) + estl_frame - l_frame;
                estr(idxStart:idxEnd) = estr(idxStart:idxEnd) + estr_frame - r_frame;
            end

            % calculate center moment of frame
            lm = moment(estl_frame, 4)/std(estl_frame)^4;
            rm = moment(estr_frame, 4)/std(estr_frame)^4;

            cml_frame = lm;
            cmr_frame = rm;

            cml(idxStart:idxEnd) = cml(idxStart:idxEnd) + cml_frame;
            cmr(idxStart:idxEnd) = cmr(idxStart:idxEnd) + cmr_frame;
            
            
            % update index for selecting next frame
            idxStart = idxStart + lenOverlap;
            idxEnd = idxEnd + lenOverlap;
            
        end
        
        % obtain residue of signal
        l_sig = estl;
        r_sig = estr;

        cml_sig = cml;
        cmr_sig = cmr;

        % check whether NaN data appear during processing
        if anynan(l_sig) == 1
            warning(strcat("NaN value appear ", string(j)))
        end

        l_sub_buf(j,:) = l_sig;
        r_sub_buf(j,:) = r_sig;

        l_cm_buf(j,:) = cml_sig;
        r_cm_buf(j,:) = cmr_sig;

    end

    env_lp_lp_data{i} = {l_sub_buf, r_sub_buf};
    env_lp_lp_cm_data{i} = {l_cm_buf, r_cm_buf};

end

%% implement inverse filtering to fine data
fine_inv_data = cell(1, length(fine_data));
for i = 1:length(fine_data)
    l_sub_buf = nan(size(fine_data{i}{1}));
    r_sub_buf = nan(size(fine_data{i}{2}));

    % split audio into several frames
    lenFrames = 0.03*fs;
    lenOverlap = 0.5*lenFrames;
    nFrames = floor(size(fine_data{i}{1},2)/(lenFrames-lenOverlap)); % the denumerator is false

    % initiate filter parameter
    nIter = 500;
    mu = 3e-9;
    p = lenFrames;              %filter order
    ss = lenOverlap;            %stepsize
    bs = ss+p-1;                  %blocksize
    
    gh = (1:p).';
    gh = gh./sqrt(sum(abs(gh).^2));
    Gh = fft([gh; zeros(ss-1,1)]);
    
    zkurt = zeros(niter,1);

    zr2 = zeros(bs,1);
    zr3 = zeros(bs,1);
    zr4 = zeros(bs,1);

    for j = 1:size(fine_data{1}{1},1)
        l = fine_data{i}{1}(j,:);
        r = fine_data{i}{2}(j,:);

        idxStart = 1;
        idxEnd = lenFrames;

        estl = zeros(size(l));
        estr = zeros(size(r));

        for m = 1:nIter
            yrn=zeros(bs,1);
            for n = 1:nFrames
                if idxEnd > length(l)
                    idxEnd = length(l);
                end

                yrn(1:p-1) = yrn(end-p+2:end);
                yrn(p:end) = l(k:k+ss-1);

                Yrn=fft(yrn);
                cYrn=conj(Yrn);
                zrn=ifft(Gh.*Yrn);

                zrn(1:p-1)=0;
                zr2(p:end)=zrn(p:end).^2;
                zr3(p:end)=zrn(p:end).^3;
                zr4(p:end)=zrn(p:end).^4;

                Z2=sum(zr2(p:end));
                Z4=sum(zr4(p:end));

                zkurt(m)=max(zkurt(m),Z4/(Z2^2+1e-15)*ss);

                z3y=fft(zr3).*cYrn;
                zy=fft(zrn).*cYrn;

                gJ=4*(Z2*z3y-Z4*zy)/(Z2^3+1e-20)*ss;
                Gh=Gh+mu*gJ;
                
                Gh=Gh./sqrt(sum(abs(Gh).^2)/bs);

                % update index for selecting next frame
                idxStart = idxStart + lenOverlap;
                idxEnd = idxEnd + lenOverlap;

            end
        
        end

        % recap all central moment value
        l_sig = estl;
        r_sig = estr;

        % check whether NaN data appear during processing
        if anynan(l_sig) == 1
            warning(strcat("NaN value appear ", string(j)))
        end

        l_sub_buf(j,:) = l_sig;
        r_sub_buf(j,:) = r_sig;

    end

    fine_inv_data{i} = {l_sub_buf, r_sub_buf};

end

%% implement inverse filtering to env data

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
for i = 1:2%length(env_lp_data)
    for j = 1:6%size(env_data{1}{1},1)
        subplot(6,1,j)
        if i == 1
            plot(t(1:length(fine_cm_data{i}{1}(j,:))), fine_cm_data{i}{1}(j,:), "Color", [0 0 0 1])
        else
            plot(t(1:length(fine_cm_data{i}{1}(j,:))), fine_cm_data{i}{1}(j,:), "Color", [0.5 0.5 0.5 0.8])
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
            plot(t(1:length(fine_cm_data{i}{1}(j+6,:))), fine_cm_data{i}{1}(j+6,:), "Color", [0 0 0 1])
        else
            plot(t(1:length(fine_cm_data{i}{1}(j+6,:))), fine_cm_data{i}{1}(j+6,:), "Color", [0.5 0.5 0.5 0.8])
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
% data_label = ["clean speech", "low reverberant speech"];
data_label = ["clean speech", "high reverberant speech"];
figure
for i = 1:2:3%length(env_lp_data)
    for j = 1:6%size(env_data{1}{1},1)
        subplot(6,1,j)
        if i == 1
            plot(t(1:length(fine_lp_cm_data{i}{1}(j,:))), fine_lp_cm_data{i}{1}(j,:), "Color", [0 0 0 1])
        else
            plot(t(1:length(fine_lp_cm_data{i}{1}(j,:))), fine_lp_cm_data{i}{1}(j,:), "Color", [0.5 0.5 0.5 0.8])
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
            plot(t(1:length(fine_lp_cm_data{i}{1}(j+6,:))), fine_lp_cm_data{i}{1}(j+6,:), "Color", [0 0 0 1])
        else
            plot(t(1:length(fine_lp_cm_data{i}{1}(j+6,:))), fine_lp_cm_data{i}{1}(j+6,:), "Color", [0.5 0.5 0.5 0.8])
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
            plot(t(1:length(env_lp_cm_data{i}{1}(j,:))), env_lp_cm_data{i}{1}(j,:), "Color", [0 0 0 1])
        else
            plot(t(1:length(env_lp_cm_data{i}{1}(j,:))), env_lp_cm_data{i}{1}(j,:), "Color", [0.5 0.5 0.5 0.8])
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
            plot(t(1:length(env_lp_cm_data{i}{1}(j+6,:))), env_lp_cm_data{i}{1}(j+6,:), "Color", [0 0 0 1])
        else
            plot(t(1:length(env_lp_cm_data{i}{1}(j+6,:))), env_lp_cm_data{i}{1}(j+6,:), "Color", [0.5 0.5 0.5 0.8])
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

%% plot the kurtosis of LP of env data
data_label = ["clean speech", "low reverberant speech"];
%data_label = ["clean speech", "high reverberant speech"];
figure
for i = 1:2%length(env_lp_data)
    for j = 1:6%size(env_data{1}{1},1)
        subplot(6,1,j)
        if i == 1
            plot(t(1:length(env_lp_lp_cm_data{i}{1}(j,:))), env_lp_lp_cm_data{i}{1}(j,:), "Color", [0 0 0 1])
        else
            plot(t(1:length(env_lp_lp_cm_data{i}{1}(j,:))), env_lp_lp_cm_data{i}{1}(j,:), "Color", [0.5 0.5 0.5 0.8])
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
            plot(t(1:length(env_lp_lp_cm_data{i}{1}(j+6,:))), env_lp_lp_cm_data{i}{1}(j+6,:), "Color", [0 0 0 1])
        else
            plot(t(1:length(env_lp_lp_cm_data{i}{1}(j+6,:))), env_lp_lp_cm_data{i}{1}(j+6,:), "Color", [0.5 0.5 0.5 0.8])
        end
        legend(data_label)
        hold on
        ylabel("Central Moment Value")
        if j == 6
            xlabel("Time (s)")
        end
    end
end