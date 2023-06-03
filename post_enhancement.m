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
            la = lpc(l_frame, 10);
            estl_frame = filter([0 -la(2:end)], 1, l_frame .* hann(length(idxStart:idxEnd))');
            l_frame = l_frame - estl_frame;

            ra = lpc(r_frame, 10);
            estr_frame = filter([0 -ra(2:end)], 1, r_frame .* hann(length(idxStart:idxEnd))');
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

%% conduct spectral substraction on fine data
fine_inv_data = cell(1, length(an_data));
for i = 1:length(fine_lp_data)
    % split audio into several frames
    lenFrames = 0.032*fs;
    lenOverlap = 0.5*lenFrames;
    nFrames = floor((size(fine_lp_data{i}{1},2) - lenOverlap)/lenOverlap) + 1;
    
    l_sub_buf = nan(size(fine_lp_data{i}{1}));
    r_sub_buf = nan(size(fine_lp_data{i}{2}));

    for j = 1:size(fine_lp_data{1}{1},1)
        lrev = fine_data{i}{1}(j,:);
        rrev = fine_data{i}{1}(j,:);

        l = fine_lp_data{i}{1}(j,:);
        r = fine_lp_data{i}{2}(j,:);

        % inverse filter application lies here
        lenFilter = lenFrames;
        nIter = 500;
        mu = 3e-9;
        
        filterOrder = lenFrames;
        stepSize = lenOverlap;
        blockSize = stepSize + filterOrder - 1;

        lgh = ones(1,length((1:filterOrder)))';
        LGh = fft([lgh; zeros(stepSize-1, 1)]);
        
        lzKurt = zeros(nIter, 1);
        lzr2 = zeros(blockSize, 1);
        lzr3 = zeros(blockSize, 1);
        lzr4 = zeros(blockSize, 1);

        rgh = (1:filterOrder)';
        RGh = fft([rgh; zeros(stepSize-1, 1)]);
        
        rzKurt = zeros(nIter, 1);
        rzr2 = zeros(blockSize, 1);
        rzr3 = zeros(blockSize, 1);
        rzr4 = zeros(blockSize, 1);

        tic
        for m = 1:nIter
            lenh_fine = zeros(blockSize, 1);
            renh_fine = zeros(blockSize, 1);
            for n = 1:nFrames
                idxStart = 1 + (n-1) * (lenFrames - lenOverlap);
                idxEnd = idxStart + lenOverlap - 1;
                
                % processing left signal
                lenh_fine(1:filterOrder-1) = lenh_fine(end-filterOrder+2:end);
                lenh_fine(filterOrder:end) = l(idxStart:idxEnd);

                LENH_FINE = fft(lenh_fine);
                cLENH_FINE = conj(LENH_FINE);
                lzrn = ifft(LGh.*LENH_FINE);

                lzrn(1:filterOrder-1) = 0;
                lzr2(filterOrder:end) = lzrn(filterOrder:end).^2;
                lzr3(filterOrder:end) = lzrn(filterOrder:end).^3;
                lzr4(filterOrder:end) = lzrn(filterOrder:end).^4;

                LZ2 = sum(lzr2(filterOrder:end));
                LZ4 = sum(lzr4(filterOrder:end));

                LZ3 = fft(lzr3).*cLENH_FINE;
                LZRN = fft(lzrn).*cLENH_FINE;

                lzKurt(m) = max(lzKurt(m), LZ4/(LZ2^2 + 1e-15)*stepSize);
                
                LgJ = 4 * (LZ2 * LZ3 - LZ4*LZRN) / ((LZ2^3 + 1e-20)*stepSize);
                LGh = LGh + mu*LgJ;

                LGh=LGh./sqrt(sum(abs(LGh).^2)/blockSize);

                % processing right signal
                renh_fine(1:filterOrder-1) = renh_fine(end-filterOrder+2:end);
                renh_fine(filterOrder:end) = r(idxStart:idxEnd);

                RENH_FINE = fft(renh_fine);
                cRENH_FINE = conj(RENH_FINE);
                rzrn = ifft(RGh.*RENH_FINE);

                rzrn(1:filterOrder-1) = 0;
                rzr2(filterOrder:end) = rzrn(filterOrder:end).^2;
                rzr3(filterOrder:end) = rzrn(filterOrder:end).^3;
                rzr4(filterOrder:end) = rzrn(filterOrder:end).^4;

                RZ2 = sum(rzr2(filterOrder:end));
                RZ4 = sum(rzr4(filterOrder:end));

                RZ3 = fft(rzr3).*cRENH_FINE;
                RZRN = fft(rzrn).*cRENH_FINE;

                rzKurt(m) = max(rzKurt(m), RZ4/(RZ2^2 + 1e-15)*stepSize);
                
                RgJ = 4 * (RZ2 * RZ3 - RZ4*RZRN) / ((RZ2^3 + 1e-20)*stepSize);
                RGh = RGh + mu*RgJ;

                RGh = RGh./sqrt(sum(abs(RGh).^2)/blockSize);

            end

        end
        toc

        lgh = ifft(LGh);
        ln = filter(lgh, 1, lrev);

        rgh = ifft(RGh);
        rn = filter(rgh, 1, rrev);

        l_sub_buf(j,:) = ln;
        r_sub_buf(j,:) = rn;
    end
    fine_inv_data{i} = {l_sub_buf, r_sub_buf};

end

%% conduct spectral substractin on fine data
fine_inv_ss_data = cell(1, length(env_data));
for i = 1:length(fine_inv_data)
    l_sub_buf = nan(size(fine_inv_data{i}{1}));
    r_sub_buf = nan(size(fine_inv_data{i}{2}));
    
    % divide signal into several frames
    lenFrames = floor(0.032*fs);
    lenOverlap = floor(0.5*lenFrames);
    nFFT = lenFrames;
    nFrames = floor((size(fine_inv_data{i}{1},2) - lenOverlap)/lenOverlap) + 1;

    % smoothing function
    ro_w=7;
    a_w=5;
    i_w=-ro_w:15;
    wS=(i_w+a_w)/(a_w^2).*exp(-.5*(i_w/a_w+1).^2);
    wS(i_w<-a_w)=0;

    gS=.32;
    epS=1e-3;

    for j = 1:size(fine_inv_data{1}{1},1)
        % spectral substraction implementation lies below
        l = fine_inv_data{i}{1}(j,:);
        r = fine_inv_data{i}{2}(j,:);

        l = l/max(abs(l));
        r = r/max(abs(r));

        [SL, ~, ~] = stft(l, fs, "Window", hamming(lenFrames), "OverlapLength", lenOverlap, "FFTLength", nFFT);
        PSL = abs(SL).^2;
        phSL = angle(SL);
        
        [SR, ~, ~] = stft(r, fs, "Window", hamming(lenFrames), "OverlapLength", lenOverlap, "FFTLength", nFFT);
        PSR = abs(SR).^2;
        phSR = angle(SR);
        
        PLL = gS * filter(wS, 1, PSL.').';
        PLR = gS * filter(wS, 1, PSR.').';

        PL_att = (PSL - PLL) ./ PSL;
        PR_att = (PSR - PLR) ./ PSR;

        PxhL = PSL .* PL_att;
        SxhL = sqrt(PxhL).*exp(1i*phSL);
        PxhR = PSR .* PR_att;
        SxhR = sqrt(PxhR).*exp(1i*phSR);

        nu1 = .05;
        nu2 = 4;

        El = sum(abs(SL).^2, 1)/size(SL,1);
        Exhl = sum(abs(SxhL).^2, 1)/size(SxhL,1);
        PL_att = ones(size(SL,2),1);
        PL_att(El<nu1 & El./Exhl>nu2)=1e-3;
        SxhL = SxhL*spdiags(PL_att,0,length(PL_att),length(PL_att));

        Er = sum(abs(SR).^2, 1)/size(SR,1);
        Exhr = sum(abs(SxhR).^2, 1)/size(SxhR,1);
        PR_att = ones(size(SR,2),1);
        PR_att(Er<nu1 & El./Exhr>nu2)=1e-3;
        SxhR = SxhR*spdiags(PR_att,0,length(PR_att),length(PR_att));
        
        [l_sig, ~] = istft(SxhL, fs, 'Window', hann(lenFrames), 'OverlapLength', lenOverlap, 'FFTLength', nFFT);
        [r_sig, ~] = istft(SxhR, fs, 'Window', hann(lenFrames), 'OverlapLength', lenOverlap, 'FFTLength', nFFT);

        l_sub_buf(j,:) = [real(l_sig)' zeros(1, length(l) - length(l_sig))];
        r_sub_buf(j,:) = [real(r_sig)' zeros(1, length(l) - length(l_sig))];

    end
    fine_inv_ss_data{i} = {l_sub_buf, r_sub_buf};

end