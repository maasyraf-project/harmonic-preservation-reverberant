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

pe_ref_data = pe_data; % for reference information

%% decompose into subband for obtain reverberant signal's information
% Create analyse -Filterbank
analyzer_stim = Gfb_Analyzer_new(par.voc_sampling_frequency_hz,...
    par.gamma_order_stimulation, ...
    par.center_frequencies_hz_stimulation,...
    par.bandwidth_factor);

an_ref_data = cell(1, length(pe_ref_data));
for i = 1:length(pe_ref_data)
    [l_sig, ~] = Gfb_Analyzer_process(analyzer_stim, pe_ref_data{i}(:,1)');
    [r_sig, ~] = Gfb_Analyzer_process(analyzer_stim, pe_ref_data{i}(:,2)');
    an_ref_data{i} = {l_sig, r_sig}; 
end

env_ref_data = cell(1, length(an_ref_data));
rms_env_ref_data = cell(1, length(an_ref_data));
fine_ref_data = cell(1, length(an_ref_data));
for i = 1:length(an_ref_data)
    % extract envelope
    el_sig = abs(an_ref_data{i}{1});
    er_sig = abs(an_ref_data{i}{2});

    if size(el_sig,2) ~= size(er_sig,2)
        error("Length of left and right analyzed signal is not same!")
    end

    weights = repmat(par.weights,1,size(el_sig,2));
    el_sig = sqrt(weights.*(real(an_ref_data{i}{1}).^2+imag(an_ref_data{i}{1}).^2));
    er_sig = sqrt(weights.*(real(an_ref_data{i}{2}).^2+imag(an_ref_data{i}{2}).^2));
    
    rmsl_sig = rms2(el_sig,2);
    rmsr_sig = rms2(er_sig,2);

    env_ref_data{i} = {el_sig, er_sig}; 
    rms_env_ref_data{i} = {rmsl_sig, rmsr_sig};
    
    % extract fine structure
    fl_sig = real(an_ref_data{i}{1});
    fr_sig = real(an_ref_data{i}{2});

    fine_ref_data{i} = {fl_sig, fr_sig};

end

env_lp_ref_data = cell(1, length(env_ref_data));
for i = 1:length(env_ref_data)
    l_sub_buf = nan(size(env_ref_data{i}{1}));
    r_sub_buf = nan(size(env_ref_data{i}{2}));

    for j = 1:size(env_ref_data{1}{1},1)
        [b, a] = butter(1, (200./(fs/2)));
        l_sig = filter(b, a, env_ref_data{i}{1}(j,:));
        r_sig = filter(b, a, env_ref_data{i}{2}(j,:));

        l_sub_buf(j,:) = l_sig;
        r_sub_buf(j,:) = r_sig;
    end
    env_lp_ref_data{i} = {l_sub_buf, r_sub_buf};

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

%% decompose into subband for obtain filter coefficient
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

%% conduct inverse filtering on fine data (LMS)
fine_inv_data = cell(1, length(an_data));
fine_rmse_data = cell(1, length(an_data));
for i = 3%:length(fine_lp_data)
    % split audio into several frames
    lenFrames = 250;%0.032*fs;
    lenOverlap = 0.2*lenFrames;
    nFrames = floor((size(fine_data{i}{1},2) - lenOverlap)/lenOverlap) + 1;
    
    l_sub_buf = nan(size(fine_data{i}{1}));
    r_sub_buf = nan(size(fine_data{i}{2}));

    l_rmse_buf = nan(size(fine_data{i}{1}, 1), 1);
    r_rmse_buf = nan(size(fine_data{i}{1}, 1), 1);

    hos = 3;

    for j = 1%:size(fine_lp_data{1}{1},1)
        lrev = fine_ref_data{i}{1}(j,:);
        rrev = fine_ref_data{i}{2}(j,:);

        lrms = rms(lrev);
        rrms = rms(lrev);

        l = fine_data{i}{1}(j,:);% ./ max(abs(fine_lp_data{i}{1}(j,:)));
        r = fine_data{i}{2}(j,:);% ./ max(abs(fine_lp_data{i}{1}(j,:)));

        % inverse filter application lies here
        lenFilter = lenFrames;
        nIter = 200;
        mu = 1e-6;
        
        filterOrder = lenFrames;
        stepSize = lenOverlap;
        blockSize = stepSize + filterOrder - 1;

        lgh = [1 zeros(1, length((1:filterOrder))-1)]';
        LGh = fft([lgh; zeros(stepSize-1, 1)]);
        
        lzKurt = zeros(nIter, 1);
        lzr2 = zeros(blockSize, 1);
        lzr3 = zeros(blockSize, 1);
        lzr4 = zeros(blockSize, 1);
        lzr5 = zeros(blockSize, 1);
        lzr6 = zeros(blockSize, 1);

        rgh = [1 zeros(1, length((1:filterOrder))-1)]';
        RGh = fft([rgh; zeros(stepSize-1, 1)]);
        
        rzKurt = zeros(nIter, 1);
        rzr2 = zeros(blockSize, 1);
        rzr3 = zeros(blockSize, 1);
        rzr4 = zeros(blockSize, 1);
        rzr5 = zeros(blockSize, 1);
        rzr6 = zeros(blockSize, 1);

        tic
        for m = 1:nIter
            lenh_fine = zeros(blockSize, 1);
            renh_fine = zeros(blockSize, 1);
            for n = 1:stepSize:length(l)
                %idxStart = 1 + (n-1) * (lenFrames - lenOverlap);
                %idxEnd = idxStart + lenOverlap - 1;

                if n+stepSize-1 > length(l)
                    li = [l(n:length(l)) zeros(1, 1)];
                    ri = [r(n:length(l)) zeros(1, 1)];
                else
                    li = l(n:n+stepSize-1);
                    ri = r(n:n+stepSize-1);
                end

                % processing left signal
                lenh_fine(1:filterOrder-1) = lenh_fine(end-filterOrder+2:end);
                lenh_fine(filterOrder:end) = li;

                LENH_FINE = fft(lenh_fine);
                cLENH_FINE = conj(LENH_FINE);
                lzrn = ifft(LGh.*LENH_FINE);

                lzrn(1:filterOrder-1) = 0;
                lzr2(filterOrder:end) = lzrn(filterOrder:end).^2;
                lzr3(filterOrder:end) = lzrn(filterOrder:end).^3;
                lzr4(filterOrder:end) = lzrn(filterOrder:end).^4;
                lzr5(filterOrder:end) = lzrn(filterOrder:end).^5;
                lzr6(filterOrder:end) = lzrn(filterOrder:end).^6;

                if hos == 3
                    LZ2 = mean(lzr2(filterOrder:end));
                    LZ3 = mean(lzr3(filterOrder:end));

                    LZ2t = fft(lzr2).*cLENH_FINE;
                    LZRN = fft(lzrn).*cLENH_FINE;

                    lzKurt(m) = max(lzKurt(m), LZ3/(LZ2^(3/2) + eps));

                    LgJ = 3 * (LZ2 * LZ2t - LZ3*LZRN) / (LZ2^3 + 1e-20);
                end

                if hos == 4
                    LZ2 = mean(lzr2(filterOrder:end));
                    LZ4 = mean(lzr4(filterOrder:end));

                    LZ3 = fft(lzr3).*cLENH_FINE;
                    LZRN = fft(lzrn).*cLENH_FINE;

                    lzKurt(m) = max(lzKurt(m), LZ4/(LZ2^2 + eps));

                    LgJ = 4 * (LZ2 * LZ3 - LZ4*LZRN) / (LZ2^3 + 1e-20);
                end

                if hos == 5
                    LZ2 = mean(lzr2(filterOrder:end));
                    LZ5 = mean(lzr5(filterOrder:end));

                    LZ4 = fft(lzr4).*cLENH_FINE;
                    LZRN = fft(lzrn).*cLENH_FINE;

                    lzKurt(m) = max(lzKurt(m), LZ5/(LZ2^(2.5) + eps));

                    LgJ = 5 * (LZ2 * LZ4 - LZ5*LZRN) / (LZ2^(3.5) + 1e-20);
                end

                if hos == 6
                    LZ2 = mean(lzr2(filterOrder:end));
                    LZ6 = mean(lzr6(filterOrder:end));

                    LZ5 = fft(lzr5).*cLENH_FINE;
                    LZRN = fft(lzrn).*cLENH_FINE;

                    lzKurt(m) = max(lzKurt(m), LZ6/(LZ2^3 + eps));

                    LgJ = 6 * (LZ2 * LZ5 - LZ6*LZRN) / (LZ2^4 + 1e-20);
                end

                % LGh = LGh + mu*LgJ / stepSize;
                LGh = LGh + (1/(eps+norm(LgJ)^2))*LgJ / stepSize;

                LGh=LGh./sqrt(sum(abs(LGh).^2)/blockSize);

                % processing right signal
                renh_fine(1:filterOrder-1) = renh_fine(end-filterOrder+2:end);
                renh_fine(filterOrder:end) = ri;

                RENH_FINE = fft(renh_fine);
                cRENH_FINE = conj(RENH_FINE);
                rzrn = ifft(RGh.*RENH_FINE);

                rzrn(1:filterOrder-1) = 0;
                rzr2(filterOrder:end) = rzrn(filterOrder:end).^2;
                rzr3(filterOrder:end) = rzrn(filterOrder:end).^3;
                rzr4(filterOrder:end) = rzrn(filterOrder:end).^4;
                rzr5(filterOrder:end) = rzrn(filterOrder:end).^5;
                rzr6(filterOrder:end) = rzrn(filterOrder:end).^6;

                if hos == 3
                    RZ2 = mean(rzr2(filterOrder:end));
                    RZ3 = mean(rzr3(filterOrder:end));

                    RZ2t = fft(rzr2).*cRENH_FINE;
                    RZRN = fft(rzrn).*cRENH_FINE;

                    rzKurt(m) = max(rzKurt(m), RZ3/(RZ2^(3/2) + eps));

                    RgJ = 3 * (RZ2 * RZ2t - RZ3*RZRN) / (RZ2^3 + 1e-20);
                end

                if hos == 4
                    RZ2 = mean(rzr2(filterOrder:end));
                    RZ4 = mean(rzr4(filterOrder:end));

                    RZ3 = fft(rzr3).*cRENH_FINE;
                    RZRN = fft(rzrn).*cRENH_FINE;

                    rzKurt(m) = max(rzKurt(m), RZ4/(RZ2^2 + eps));

                    RgJ = 4 * (RZ2 * RZ3 - RZ4*RZRN) / (RZ2^3 + 1e-20);
                end

                if hos == 5
                    RZ2 = mean(rzr2(filterOrder:end));
                    RZ5 = mean(rzr5(filterOrder:end));

                    RZ4 = fft(rzr4).*cRENH_FINE;
                    RZRN = fft(rzrn).*cRENH_FINE;

                    rzKurt(m) = max(rzKurt(m), RZ5/(RZ2^(2.5) + eps));

                    RgJ = 5 * (RZ2 * RZ4 - RZ5*RZRN) / (RZ2^(3.5) + 1e-20);
                end

                if hos == 6
                    RZ2 = mean(rzr2(filterOrder:end));
                    RZ6 = mean(rzr6(filterOrder:end));

                    RZ5 = fft(rzr5).*cLENH_FINE;
                    RZRN = fft(rzrn).*cLENH_FINE;

                    rzKurt(m) = max(rzKurt(m), RZ6/(RZ2^3 + eps));

                    RgJ = 6 * (RZ2 * RZ5 - RZ6*RZRN) / (RZ2^4 + 1e-20);
                end

                % RGh = RGh + mu*RgJ / stepSize;
                RGh = RGh + (1/(eps+norm(RgJ)^2))*RgJ / stepSize;

                RGh=RGh./sqrt(sum(abs(RGh).^2)/blockSize);

            end

        end
        toc

        lgh = ifft(LGh);
        ln = filter(lgh, 1, lrev);

        ln = ln * lrms/rms(ln);

        rgh = ifft(RGh);
        rn = filter(rgh, 1, rrev);

        rn = rn * lrms/rms(rn);

        l_sub_buf(j,:) = ln;
        r_sub_buf(j,:) = rn;

        % calculate error
        l_rmse_buf(j,:) = sqrt(mean((ln-fine_ref_data{1}{1}(j,:)).^2));
        r_rmse_buf(j,:) = sqrt(mean((rn-fine_ref_data{1}{2}(j,:)).^2));

    end
    fine_inv_data{i} = {l_sub_buf, r_sub_buf};
    fine_rmse_data{i} = {l_rmse_buf, r_rmse_buf};

end

%% conduct inverse filtering on env data (LMS)
env_inv_data = cell(1, length(an_data));
env_rmse_data = cell(1, length(an_data));
for i = 3%:length(fine_lp_data)
    % split audio into several frames
    lenFrames = 250;%0.032*fs;
    lenOverlap = 0.2*lenFrames;
    nFrames = floor((size(env_data{i}{1},2) - lenOverlap)/lenOverlap) + 1;
    
    l_sub_buf = nan(size(env_data{i}{1}));
    r_sub_buf = nan(size(env_data{i}{2}));

    l_rmse_buf = nan(size(env_data{i}{1}, 1), 1);
    r_rmse_buf = nan(size(env_data{i}{1}, 1), 1);

    hos = 6;

    for j = 1%:size(fine_lp_data{1}{1},1)
        lrev = env_ref_data{i}{1}(j,:);
        rrev = env_ref_data{i}{2}(j,:);

        lrms = rms(lrev);
        rrms = rms(lrev);

        l = env_data{i}{1}(j,:);% ./ max(abs(fine_lp_data{i}{1}(j,:)));
        r = env_data{i}{2}(j,:);% ./ max(abs(fine_lp_data{i}{1}(j,:)));

        % inverse filter application lies here
        lenFilter = lenFrames;
        nIter = 250;
        mu = 1e-6;
        
        filterOrder = lenFrames;
        stepSize = lenOverlap;
        blockSize = stepSize + filterOrder - 1;

        lgh = [1 zeros(1, length((1:filterOrder))-1)]';
        LGh = fft([lgh; zeros(stepSize-1, 1)]);
        
        lzKurt = zeros(nIter, 1);
        lzr2 = zeros(blockSize, 1);
        lzr3 = zeros(blockSize, 1);
        lzr4 = zeros(blockSize, 1);
        lzr5 = zeros(blockSize, 1);
        lzr6 = zeros(blockSize, 1);

        rgh = [1 zeros(1, length((1:filterOrder))-1)]';
        RGh = fft([rgh; zeros(stepSize-1, 1)]);
        
        rzKurt = zeros(nIter, 1);
        rzr2 = zeros(blockSize, 1);
        rzr3 = zeros(blockSize, 1);
        rzr4 = zeros(blockSize, 1);
        rzr5 = zeros(blockSize, 1);
        rzr6 = zeros(blockSize, 1);

        tic
        for m = 1:nIter
            lenh_fine = zeros(blockSize, 1);
            renh_fine = zeros(blockSize, 1);
            for n = 1:stepSize:length(l)
                %idxStart = 1 + (n-1) * (lenFrames - lenOverlap);
                %idxEnd = idxStart + lenOverlap - 1;

                if n+stepSize-1 > length(l)
                    li = [l(n:length(l)) zeros(1, 1)];
                    ri = [r(n:length(l)) zeros(1, 1)];
                else
                    li = l(n:n+stepSize-1);
                    ri = r(n:n+stepSize-1);
                end

                % processing left signal
                lenh_fine(1:filterOrder-1) = lenh_fine(end-filterOrder+2:end);
                lenh_fine(filterOrder:end) = li;

                LENH_FINE = fft(lenh_fine);
                cLENH_FINE = conj(LENH_FINE);
                lzrn = ifft(LGh.*LENH_FINE);

                lzrn(1:filterOrder-1) = 0;
                lzr2(filterOrder:end) = lzrn(filterOrder:end).^2;
                lzr3(filterOrder:end) = lzrn(filterOrder:end).^3;
                lzr4(filterOrder:end) = lzrn(filterOrder:end).^4;
                lzr5(filterOrder:end) = lzrn(filterOrder:end).^5;
                lzr6(filterOrder:end) = lzrn(filterOrder:end).^6;

                if hos == 3
                    LZ2 = mean(lzr2(filterOrder:end));
                    LZ3 = mean(lzr3(filterOrder:end));

                    LZ2t = fft(lzr2).*cLENH_FINE;
                    LZRN = fft(lzrn).*cLENH_FINE;

                    lzKurt(m) = max(lzKurt(m), LZ3/(LZ2^(3/2) + eps));

                    LgJ = 3 * (LZ2 * LZ2t - LZ3*LZRN) / (LZ2^3 + 1e-20);
                end

                if hos == 4
                    LZ2 = mean(lzr2(filterOrder:end));
                    LZ4 = mean(lzr4(filterOrder:end));

                    LZ3 = fft(lzr3).*cLENH_FINE;
                    LZRN = fft(lzrn).*cLENH_FINE;

                    lzKurt(m) = max(lzKurt(m), LZ4/(LZ2^2 + eps));

                    LgJ = 4 * (LZ2 * LZ3 - LZ4*LZRN) / (LZ2^3 + 1e-20);
                end

                if hos == 5
                    LZ2 = mean(lzr2(filterOrder:end));
                    LZ5 = mean(lzr5(filterOrder:end));

                    LZ4 = fft(lzr4).*cLENH_FINE;
                    LZRN = fft(lzrn).*cLENH_FINE;

                    lzKurt(m) = max(lzKurt(m), LZ5/(LZ2^(2.5) + eps));

                    LgJ = 5 * (LZ2 * LZ4 - LZ5*LZRN) / (LZ2^(3.5) + 1e-20);
                end

                if hos == 6
                    LZ2 = mean(lzr2(filterOrder:end));
                    LZ6 = mean(lzr6(filterOrder:end));

                    LZ5 = fft(lzr5).*cLENH_FINE;
                    LZRN = fft(lzrn).*cLENH_FINE;

                    lzKurt(m) = max(lzKurt(m), LZ6/(LZ2^3 + eps));

                    LgJ = 6 * (LZ2 * LZ5 - LZ6*LZRN) / (LZ2^4 + 1e-20);
                end

                % LGh = LGh + mu*LgJ / stepSize;
                LGh = LGh + (1/(eps+norm(LgJ)^2))*LgJ / stepSize;

                LGh=LGh./sqrt(sum(abs(LGh).^2)/blockSize);

                % processing right signal
                renh_fine(1:filterOrder-1) = renh_fine(end-filterOrder+2:end);
                renh_fine(filterOrder:end) = ri;

                RENH_FINE = fft(renh_fine);
                cRENH_FINE = conj(RENH_FINE);
                rzrn = ifft(RGh.*RENH_FINE);

                rzrn(1:filterOrder-1) = 0;
                rzr2(filterOrder:end) = rzrn(filterOrder:end).^2;
                rzr3(filterOrder:end) = rzrn(filterOrder:end).^3;
                rzr4(filterOrder:end) = rzrn(filterOrder:end).^4;
                rzr5(filterOrder:end) = rzrn(filterOrder:end).^5;
                rzr6(filterOrder:end) = rzrn(filterOrder:end).^6;

                if hos == 3
                    RZ2 = mean(rzr2(filterOrder:end));
                    RZ3 = mean(rzr3(filterOrder:end));

                    RZ2t = fft(rzr2).*cRENH_FINE;
                    RZRN = fft(rzrn).*cRENH_FINE;

                    rzKurt(m) = max(rzKurt(m), RZ3/(RZ2^(3/2) + eps));

                    RgJ = 3 * (RZ2 * RZ2t - RZ3*RZRN) / (RZ2^3 + 1e-20);
                end

                if hos == 4
                    RZ2 = mean(rzr2(filterOrder:end));
                    RZ4 = mean(rzr4(filterOrder:end));

                    RZ3 = fft(rzr3).*cRENH_FINE;
                    RZRN = fft(rzrn).*cRENH_FINE;

                    rzKurt(m) = max(rzKurt(m), RZ4/(RZ2^2 + eps));

                    RgJ = 4 * (RZ2 * RZ3 - RZ4*RZRN) / (RZ2^3 + 1e-20);
                end

                if hos == 5
                    RZ2 = mean(rzr2(filterOrder:end));
                    RZ5 = mean(rzr5(filterOrder:end));

                    RZ4 = fft(rzr4).*cRENH_FINE;
                    RZRN = fft(rzrn).*cRENH_FINE;

                    rzKurt(m) = max(rzKurt(m), RZ5/(RZ2^(2.5) + eps));

                    RgJ = 5 * (RZ2 * RZ4 - RZ5*RZRN) / (RZ2^(3.5) + 1e-20);
                end

                if hos == 6
                    RZ2 = mean(rzr2(filterOrder:end));
                    RZ6 = mean(rzr6(filterOrder:end));

                    RZ5 = fft(rzr5).*cLENH_FINE;
                    RZRN = fft(rzrn).*cLENH_FINE;

                    rzKurt(m) = max(rzKurt(m), RZ6/(RZ2^3 + eps));

                    RgJ = 6 * (RZ2 * RZ5 - RZ6*RZRN) / (RZ2^4 + 1e-20);
                end

                % RGh = RGh + mu*RgJ / stepSize;
                RGh = RGh + (1/(eps+norm(RgJ)^2))*RgJ / stepSize;

                RGh=RGh./sqrt(sum(abs(RGh).^2)/blockSize);

            end

        end
        toc

        lgh = ifft(LGh);
        ln = filter(lgh, 1, lrev);

        ln = ln * lrms/rms(ln);

        rgh = ifft(RGh);
        rn = filter(rgh, 1, rrev);

        rn = rn * lrms/rms(rn);

        l_sub_buf(j,:) = ln;
        r_sub_buf(j,:) = rn;

        % calculate error
        l_rmse_buf(j,:) = sqrt(mean((ln-env_ref_data{1}{1}(j,:)).^2));
        r_rmse_buf(j,:) = sqrt(mean((rn-env_ref_data{1}{2}(j,:)).^2));

    end
    env_inv_data{i} = {l_sub_buf, r_sub_buf};
    env_rmse_data{i} = {l_rmse_buf, r_rmse_buf};

end

%% conduct fine spectral substraction
fine_ss_data = cell(1, length(an_data));
for i = 3%:length(fine_lp_data)
    % split audio into several frames
    wlen = 2048;
    hop = 1024;
    nfft = 2048;

    l_sub_buf = zeros(size(fine_data{i}{1}));
    r_sub_buf = zeros(size(fine_data{i}{2}));

    for j = 1%:size(fine_lp_data{1}{1},1)
        l = fine_inv_data{i}{1}(j,:);% ./ max(abs(fine_lp_data{i}{1}(j,:)));
        r = fine_inv_data{i}{2}(j,:);% ./ max(abs(fine_lp_data{i}{1}(j,:)));

        ch = [l; r];
        ss = [];
        
        for c = 1:size(ch,1)
        % apply spectral substraction
        [Sz, f_stft, t_stft] = stft(ch(c,:), fs, Window=hamming(wlen), OverlapLength=hop, FFTLength=nfft);
        Pz=abs(Sz).^2;
        phz=angle(Sz);

        ro_w=7;
        a_w=5;
        i_w=-ro_w:15;
        wS=(i_w+a_w)/(a_w^2).*exp(-.5*(i_w/a_w+1).^2);
        wS(i_w<-a_w)=0;
        % figure;plot(i_w,wS)

        gS=.32;
        epS=1e-3;
        Pl=gS*filter(wS,1,Pz.').';

        P_att=(Pz-Pl)./Pz;
        P_att(P_att<epS)=epS;

        Pxh=Pz.*P_att;
        Sxh=sqrt(Pxh).*exp(1i*phz);

        nu1 = .05;
        nu2 = 4;
        Ez=sum(abs(Sz).^2, 1)/size(Sz,1);
        Exh=sum(abs(Sxh).^2, 1)/size(Sxh,1);

        P_att=ones(size(Sz,2),1);
        P_att(Ez<nu1 & Ez./Exh>nu2)=1e-3;
        Sxh = Sxh*spdiags(P_att,0,length(P_att),length(P_att));

        [xh, t_istft] = istft(Sxh, fs, Window=hamming(wlen), OverlapLength=hop, FFTLength=nfft);
            
        ss = [ss; xh'];
        
        end

        l_sub_buf(j,1:length(ss)) = real(ss(1,:));
        r_sub_buf(j,1:length(ss)) = real(ss(2,:));

    end
    fine_ss_data{i} = {l_sub_buf, r_sub_buf};

end

%% conduct env spectral substraction
env_ss_data = cell(1, length(an_data));
for i = 3%:length(fine_lp_data)
    % split audio into several frames
    wlen = 2048;
    hop = 1024;
    nfft = 2048;

    l_sub_buf = zeros(size(fine_data{i}{1}));
    r_sub_buf = zeros(size(fine_data{i}{2}));

    for j = 1%:size(fine_lp_data{1}{1},1)
        l = env_inv_data{i}{1}(j,:);% ./ max(abs(fine_lp_data{i}{1}(j,:)));
        r = env_inv_data{i}{2}(j,:);% ./ max(abs(fine_lp_data{i}{1}(j,:)));

        ch = [l; r];
        ss = [];
        
        for c = 1:size(ch,1)
        % apply spectral substraction
        [Sz, f_stft, t_stft] = stft(ch(c,:), fs, Window=hamming(wlen), OverlapLength=hop, FFTLength=nfft);
        Pz=abs(Sz).^2;
        phz=angle(Sz);

        ro_w=7;
        a_w=5;
        i_w=-ro_w:15;
        wS=(i_w+a_w)/(a_w^2).*exp(-.5*(i_w/a_w+1).^2);
        wS(i_w<-a_w)=0;
        % figure;plot(i_w,wS)

        gS=.32;
        epS=1e-3;
        Pl=gS*filter(wS,1,Pz.').';

        P_att=(Pz-Pl)./Pz;
        P_att(P_att<epS)=epS;

        Pxh=Pz.*P_att;
        Sxh=sqrt(Pxh).*exp(1i*phz);

        nu1 = .05;
        nu2 = 4;
        Ez=sum(abs(Sz).^2, 1)/size(Sz,1);
        Exh=sum(abs(Sxh).^2, 1)/size(Sxh,1);

        P_att=ones(size(Sz,2),1);
        P_att(Ez<nu1 & Ez./Exh>nu2)=1e-3;
        Sxh = Sxh*spdiags(P_att,0,length(P_att),length(P_att));

        [xh, t_istft] = istft(Sxh, fs, Window=hamming(wlen), OverlapLength=hop, FFTLength=nfft);
            
        ss = [ss; xh'];
        
        end

        l_sub_buf(j,1:length(ss)) = real(ss(1,:));
        r_sub_buf(j,1:length(ss)) = real(ss(2,:));

    end
    fine_ss_data{i} = {l_sub_buf, r_sub_buf};

end

%% customized plot
% fine low-rev
figure
subplot(3,2,1)
plot(t(1:length(fine_data{1}{1}(1,:))), fine_ref_data{1}{1}(1,:), "Color", [0 0 0 1])
hold on 
plot(t(1:length(fine_data{2}{1}(1,:))), fine_ss_data{3}{1}(1,:), "Color", [0.5 0.5 0.5 0.4])
ylabel("Amplitude")
title("Fine Structure of Subband (Left Channel)")

subplot(3,2,2)
plot(t(1:length(fine_data{1}{2}(1,:))), fine_ref_data{1}{2}(1,:), "Color", [0 0 0 1])
hold on
plot(t(1:length(fine_data{2}{2}(1,:))), fine_ss_data{3}{2}(1,:), "Color", [0.5 0.5 0.5 0.4])
ylabel("Amplitude")
title("Fine Structure of Subband (Right Channel)")

subplot(3,2,3)
plot(t(1:length(fine_data{1}{1}(2,:))), fine_ref_data{1}{1}(2,:), "Color", [0 0 0 1])
hold on
plot(t(1:length(fine_data{2}{1}(2,:))), fine_inv_data{3}{1}(2,:), "Color", [0.5 0.5 0.5 0.4])
ylabel("Amplitude")

subplot(3,2,4)
plot(t(1:length(fine_data{1}{2}(2,:))), fine_ref_data{1}{2}(2,:), "Color", [0 0 0 1])
hold on
plot(t(1:length(fine_data{2}{2}(2,:))), fine_inv_data{3}{2}(2,:), "Color", [0.5 0.5 0.5 0.4])
ylabel("Amplitude")

subplot(3,2,5)
plot(t(1:length(fine_data{1}{1}(3,:))), fine_ref_data{1}{1}(3,:), "Color", [0 0 0 1])
hold on
plot(t(1:length(fine_data{2}{1}(3,:))), fine_inv_data{3}{1}(3,:), "Color", [0.5 0.5 0.5 0.4])
xlabel("Time (s)")
ylabel("Amplitude")

subplot(3,2,6)
plot(t(1:length(fine_data{1}{2}(3,:))), fine_ref_data{i}{2}(3,:), "Color", [0 0 0 1])
hold on
plot(t(1:length(fine_data{2}{2}(3,:))), fine_inv_data{3}{2}(3,:), "Color", [0.5 0.5 0.5 0.4])
ylabel("Amplitude")
xlabel("Time (s)")
