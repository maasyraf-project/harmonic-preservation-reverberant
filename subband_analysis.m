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


%% extract linear prediction of envelope
env_lp_lpc_data = cell(1, length(env_lp_data));
sig_vad = cell(1, length(raw_data));
for i = 1:length(env_data)
    l_sub_buf = nan(size(env_lp_data{i}{1}));
    r_sub_buf = nan(size(env_lp_data{i}{2}));

    % split audio into several frames
    lenFrames = 0.03*fs;
    lenOverlap = 0.5*lenFrames;
    nFrames = ceil(size(env_lp_data{i}{1},2)/lenFrames);

    l_vad_buf = nan(size(env_lp_data{i}{1},1),nFrames);
    r_vad_buf = nan(size(env_lp_data{i}{2},1),nFrames);

    for j = 1:size(env_data{1}{1},1)
        lv = fine_data{i}{1}(j,:);
        rv = fine_data{i}{2}(j,:);

        l = env_lp_data{i}{1}(j,:);
        r = env_lp_data{i}{2}(j,:);

        % conduct VAD
        VADPar = InitVADPar(fs,lenFrames,lenOverlap);
        idxStart = 1;
        idxEnd = lenOverlap;

        l_vad = [];
        r_vad = [];

        estl = zeros(size(l));
        estr = zeros(size(r));

        % y_vad = zeros(size(l));

        for n = 1:nFrames
            % calculate VAD
            lv_frame = lv(idxStart:idxEnd);
            rv_frame = rv(idxStart:idxEnd);
            
            if or(length(lv_frame) < lenOverlap, length(rv_frame) < lenOverlap)
                lv_frame = [lv_frame zeros(1, lenOverlap - length(lv_frame))];
                rv_frame = [rv_frame zeros(1, lenOverlap - length(rv_frame))];
            end
           
            [l_voiced,VADPar] = VAD(lv_frame', VADPar);
            [r_voiced,VADPar] = VAD(rv_frame', VADPar);

            l_vad = [l_vad, l_voiced];
            r_vad = [r_vad, r_voiced];
            
            % conduct LP-analysis
            l_frame = l(idxStart:idxEnd);
            r_frame = r(idxStart:idxEnd);

            if or(length(l_frame) < lenOverlap, length(r_frame) < lenOverlap)
                l_frame = [l_frame zeros(1, lenOverlap - length(l_frame))];
                r_frame = [r_frame zeros(1, lenOverlap - length(r_frame))];
            end

            if l_voiced == 1 
                % apply LPC on voiced segments
                la = lpc(l_frame, 12);
                estl_frame = filter([0 -la(2:end)], 1, l_frame);
                estl(idxStart:idxEnd) = estl(idxStart:idxEnd) + estl_frame;
            else 
                estl(idxStart:idxEnd) = estl(idxStart:idxEnd) + l(idxStart:idxEnd);
            end

            if r_voiced == 1 
                % apply LPC on voiced segments
                ra = lpc(r_frame, 12);
                estr_frame = filter([0 -ra(2:end)], 1, r_frame);
                estr(idxStart:idxEnd) = estr(idxStart:idxEnd) + estr_frame;

            else
                estr(idxStart:idxEnd) = estr(idxStart:idxEnd) + r(idxStart:idxEnd);
            end
            
            % update index for selecting next frame
            idxStart = idxStart + lenOverlap;
            idxEnd = idxEnd + lenOverlap;
        
        end

        l_vad_buf(j,:) = l_vad;
        r_vad_buf(j,:) = r_vad;
        
        % calculating LP residual signal
        l_sig = l - estl;
        r_sig = r - estr;

        % check whether NaN data appear during processing
        if anynan(l_sig) == 1
            warning(strcat("NaN value appear ", string(j)))
        end

        l_sub_buf(j,:) = l_sig;
        r_sub_buf(j,:) = r_sig;

    end
    
    sig_vad{i} = {l_vad_buf, r_vad_buf};
    env_lp_lpc_data{i} = {l_sub_buf, r_sub_buf};

end

%% extract kurtosis data of LP-residual of envelope
env_lp_kurt_data = cell(1, length(env_lp_lpc_data));
for i = 1:length(env_lp_lpc_data)
    l_sub_buf = nan(size(env_lp_data{i}{1}));
    r_sub_buf = nan(size(env_lp_data{i}{2}));

    for j = 1:size(env_data{1}{1},1)
        lseg = 0.03*fs;
        theta = 1e-6;
        l_sig = kurt(env_lp_lpc_data{i}{1}(j,:), lseg, theta);
        r_sig = kurt(env_lp_lpc_data{i}{2}(j,:), lseg, theta);

        l_sub_buf(j,:) = l_sig;
        r_sub_buf(j,:) = r_sig;
    end
    env_lp_kurt_data{i} = {l_sub_buf, r_sub_buf};

end

%%
env_lp_lpc_kurt_data = cell(1, length(env_data));
sig_vad = cell(1, length(raw_data));
for i = 1:length(env_data)
    l_sub_buf = nan(size(env_lp_lpc_data{i}{1}));
    r_sub_buf = nan(size(env_lp_lpc_data{i}{2}));

    % split audio into several frames
    lenFrames = 0.03*fs;
    lenOverlap = 0.5*lenFrames;
    nFrames = ceil(size(env_lp_lpc_data{i}{1},2)/lenFrames);

    l_vad_buf = nan(size(env_lp_lpc_data{i}{1},1),nFrames);
    r_vad_buf = nan(size(env_lp_lpc_data{i}{2},1),nFrames);

    for j = 1:size(env_lp_lpc_data{1}{1},1)
        lv = env_lp_lpc_data{i}{1}(j,:);
        rv = env_lp_lpc_data{i}{2}(j,:);

        l = env_lp_lpc_data{i}{1}(j,:);
        r = env_lp_lpc_data{i}{2}(j,:);

        % conduct VAD
        VADPar = InitVADPar(fs,lenFrames,lenOverlap);
        idxStart = 1;
        idxEnd = lenOverlap;

        l_vad = [];
        r_vad = [];

        estl = zeros(size(l));
        estr = zeros(size(r));

        % y_vad = zeros(size(l));

        for n = 1:nFrames
            % calculate VAD
            lv_frame = lv(idxStart:idxEnd);
            rv_frame = rv(idxStart:idxEnd);
            
            if or(length(lv_frame) < lenOverlap, length(rv_frame) < lenOverlap)
                lv_frame = [lv_frame zeros(1, lenOverlap - length(lv_frame))];
                rv_frame = [rv_frame zeros(1, lenOverlap - length(rv_frame))];
            end
           
            [l_voiced,VADPar] = VAD(lv_frame', VADPar);
            [r_voiced,VADPar] = VAD(rv_frame', VADPar);

            l_vad = [l_vad, l_voiced];
            r_vad = [r_vad, r_voiced];
            
            % extract kurtosis
            l_frame = l(idxStart:idxEnd);
            r_frame = r(idxStart:idxEnd);

            if or(length(l_frame) < lenOverlap, length(r_frame) < lenOverlap)
                l_frame = [l_frame zeros(1, lenOverlap - length(l_frame))];
                r_frame = [r_frame zeros(1, lenOverlap - length(r_frame))];
            end

            if l_voiced == 1 
                % extract kurtosis on voiced segments
                estl_frame = kurtosis(l_frame);
                estl_frame(isnan(estl_frame)) = 0;
                estl(idxStart:idxEnd) = estl_frame;
            else 
                estl(idxStart:idxEnd) = estl(idxStart:idxEnd) + l(idxStart:idxEnd);
            end

            if r_voiced == 1 
                % extract kurtosis on voiced segments
                estr_frame = kurtosis(r_frame);
                estr_frame(isnan(estr_frame)) = 0;
                estr(idxStart:idxEnd) = estr_frame;
            else
                estr(idxStart:idxEnd) = estr(idxStart:idxEnd) + r(idxStart:idxEnd);
            end
            
            % update index for selecting next frame
            idxStart = idxStart + lenOverlap;
            idxEnd = idxEnd + lenOverlap;
        
        end

        l_vad_buf(j,:) = l_vad;
        r_vad_buf(j,:) = r_vad;
        
        % calculating LP residual signal
        l_sig = l - estl;
        r_sig = r - estr;

        % check whether NaN data appear during processing
        if anynan(l_sig) == 1
            warning(strcat("NaN value appear ", string(j)))
        end

        l_sub_buf(j,:) = l_sig;
        r_sub_buf(j,:) = r_sig;

    end
    
    sig_vad = {l_vad_buf, r_vad_buf};
    env_lp_lpc_kurt_data{i} = {l_sub_buf, r_sub_buf};

end

%% conduct adaptive filtering on envelope
env_filtered_data = cell(1, length(env_lp_data));
rms_env_filtered_data = cell(1, length(env_lp_data));
for i = 2:length(env_lp_kurt_data)
    l_sub_buf = nan(size(env_lp_data{i}{1}));
    r_sub_buf = nan(size(env_lp_data{i}{2}));

    for j = 1:size(env_data{1}{1},1)
        % adjust parameter
        Lf=250;
        niter=150;
        mu=3e-6;
        p=Lf;                       %filter order
        ss=Lf;                      %stepsize
        bs=ss+p-1;                  %blocksize
        % gh=[1;zeros(p-1,1)];
        gh=(1:p).';
        gh=gh./sqrt(sum(abs(gh).^2));
        Gh=fft([gh; zeros(ss-1,1)]);
        zkurt=zeros(niter,1);

        zr2=zeros(bs,1);
        zr3=zeros(bs,1);
        zr4=zeros(bs,1);

        disp("processing Left Channel")
        tic
        for m=1:niter
            yrn=zeros(bs,1);
            for k=1:ss:size(env_lp_kurt_data{i}{1},2)
                yrn(1:p-1)=yrn(end-p+2:end);

                if k+ss-1 < size(env_lp_kurt_data{i}{1},2)
                    yrn(p:end) = env_lp_kurt_data{i}{1}(j,k:k+ss-1);
                else
                    env_padded = [env_lp_kurt_data{i}{1}(j,k:end) zeros(1, size(yrn(p:end),1)-size(env_lp_kurt_data{i}{1}(j,k:end),2))];
                    yrn(p:end) = env_padded';
                end

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

                z3y=ifft(fft(zr3).*cYrn);
                z3y(p+1:end)=0;
                zy=ifft(fft(zrn).*cYrn);
                zy(p+1:end)=0;

                gJ=4*(Z2*z3y-Z4*zy)/(Z2^3+1e-20)*ss;
                Gh=Gh+mu*fft(gJ);


                Gh=Gh./sqrt(sum(abs(Gh).^2)/bs);
            end
        end
        toc

        gh = ifft(Gh);
        gh = gh(1:p);

        lz = filter(gh, 1, env_lp_data{i}{1}(j,:));
        l_sub_buf(j,:) = lz;

        % processing right channel
        clear gh Gh
        gh=(1:p).';
        gh=gh./sqrt(sum(abs(gh).^2));
        Gh=fft([gh; zeros(ss-1,1)]);
        zkurt=zeros(niter,1);

        zr2=zeros(bs,1);
        zr3=zeros(bs,1);
        zr4=zeros(bs,1);

        disp("processing Right Channel")
        tic
        for m=1:niter
            yrn=zeros(bs,1);
            for k=1:ss:size(env_lp_kurt_data{i}{2},2)
                yrn(1:p-1)=yrn(end-p+2:end);

                if k+ss-1 < size(env_lp_kurt_data{i}{2},2)
                    yrn(p:end) = env_lp_kurt_data{i}{2}(j,k:k+ss-1);
                else
                    env_padded = [env_lp_kurt_data{i}{2}(j,k:end) zeros(1, size(yrn(p:end),1)-size(env_lp_kurt_data{i}{2}(j,k:end),2))];
                    yrn(p:end) = env_padded';
                end

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

                z3y=ifft(fft(zr3).*cYrn);
                z3y(p+1:end)=0;
                zy=ifft(fft(zrn).*cYrn);
                zy(p+1:end)=0;

                gJ=4*(Z2*z3y-Z4*zy)/(Z2^3+1e-20)*ss;
                Gh=Gh+mu*fft(gJ);


                Gh=Gh./sqrt(sum(abs(Gh).^2)/bs);
            end
        end
        toc

        gh = ifft(Gh);
        gh = gh(1:p);

        rz = filter(gh, 1, env_lp_data{i}{2}(j,:));
        r_sub_buf(j,:) = rz;

    end
    env_filtered_data{i} = {l_sub_buf, r_sub_buf};

    rmsl = rms2(l_sub_buf, 2);
    rmsr = rms2(r_sub_buf, 2);
    rms_env_filtered_data{i} = {rmsl, rmsr};

end

% adjust RMS of signal
for i = 2:length(env_filtered_data)
    l_sub_buf = nan(size(env_filtered_data{i}{1}));
    r_sub_buf = nan(size(env_filtered_data{i}{2}));

    for j = 1:size(env_data{1}{1},1)
        l_sub_buf(j,:) = (rms_env_data{i}{1}(j) / rms_env_filtered_data{i}{1}(j)) * env_filtered_data{i}{1}(j,:);
        r_sub_buf(j,:) = (rms_env_data{i}{2}(j) / rms_env_filtered_data{i}{2}(j)) * env_filtered_data{i}{2}(j,:);
    end
    env_filtered_data{i} = {l_sub_buf, r_sub_buf};

end

% conduct spectral substraction
env_ss_data = cell(1, length(env_data));
for i = 2:length(env_filtered_data)
    l_sub_buf = nan(size(env_lp_data{i}{1}));
    r_sub_buf = nan(size(env_lp_data{i}{2}));

    for j = 1:size(env_filtered_data{2}{1},1)
        % define parameter
        [L, f, t] = stft(env_filtered_data{i}{1}(j,:),fs, Window=hanning(400), OverlapLength=200, FFTLength=400);
        [R, f, t] = stft(env_filtered_data{i}{2}(j,:),fs, Window=hanning(400), OverlapLength=200, FFTLength=400);
        
        L_pow = abs(L).^2;
        R_pow = abs(R).^2;

        L_phase = angle(L);
        R_phase = angle(R);
        
        ro_w=7;
        a_w=5;
        i_w=-ro_w:15;
        wS=(i_w+a_w)/(a_w^2).*exp(-.5*(i_w/a_w+1).^2);
        wS(i_w<-a_w)=0;

        gS=.32;
        epS=1e-3;
        L_sub = gS*filter(wS,1,L_pow.').';
        R_sub = gS*filter(wS,1,R_pow.').';

        LP_att=(L_pow-L_sub)./(L_pow+eps);
        LP_att(LP_att<epS)=epS;

        RP_att=(R_pow-R_sub)./(R_pow+eps);
        RP_att(RP_att<epS)=epS;
        
        % predict enhanced specstrum
        SL_pow = L_pow .* LP_att;
        SR_pow = R_pow .* RP_att;

        SL = sqrt(SL_pow).*exp(1i*L_phase);
        SR = sqrt(SR_pow).*exp(1i*R_phase);
        
        % inverse transfrom to time domain
        nu1 = .05;
        nu2 = 4;

        EL = sum(abs(L).^2, 1)/size(L,1);
        ESL = sum(abs(SL).^2, 1)/size(SL,1);

        ER = sum(abs(R).^2, 1)/size(R,1);
        ESR = sum(abs(SR).^2, 1)/size(SR,1);

        LP_att=ones(size(L,2),1);
        LP_att(EL<nu1 & EL./ESL>nu2) = 1e-3;

        RP_att=ones(size(R,2),1);
        RP_att(EL<nu1 & EL./ESL>nu2) = 1e-3;

        SL = SL*spdiags(LP_att,0,length(LP_att),length(LP_att));
        SR = SR*spdiags(RP_att,0,length(RP_att),length(RP_att));

        [sl, t_istft] = istft(SL, fs, Window=hanning(400), OverlapLength=200, FFTLength=400);
        [sr, t_istft] = istft(SR, fs, Window=hanning(400), OverlapLength=200, FFTLength=400);

        l_sig = sl/max(abs(sl));
        r_sig = sr/max(abs(sr));

        % spectrum analysis
        l_sub_buf(j,:) = [l_sig' zeros(1, length(env_filtered_data{i}{1}(j,:))-length(l_sig))];
        r_sub_buf(j,:) = [r_sig' zeros(1, length(env_filtered_data{i}{2}(j,:))-length(r_sig))];
    end
    env_ss_data{i} = {l_sub_buf, r_sub_buf};

end

lags_env_data = cell(1, length(env_lp_data));
for i = 2:length(env_lp_kurt_data)
    l_sub_buf = nan(size(env_lp_data{i}{1},1),1);
    r_sub_buf = nan(size(env_lp_data{i}{2},1),1);
    for k = 1:size(env_data{1}{1},1)
        [coef, lags] = xcorr(env_filtered_data{i}{1}(k,:), env_lp_data{i}{1}(k,:));
        [~, idx] = max((coef));
        l = lags(idx)/fs;

        [coef, lags] = xcorr(env_filtered_data{i}{2}(k,:), env_lp_data{i}{2}(k,:));
        [~, idx] = max((coef));
        disp(idx)
        r = lags(idx)/fs;
        
        l_sub_buf(k) = l;
        r_sub_buf(k) = r;

    end
    lags_env_data{i} = {l_sub_buf, r_sub_buf};

end

%% calculate similarity between each envelope


%% reconstruction / synthesis part
% Create analyse -Filterbank
analyzer_aural = Gfb_Analyzer_new(par.voc_sampling_frequency_hz,...
    par.gamma_order_auralisation,...
    par.center_frequencies_hz_auralisation,...
    par.bandwidth_factor);

enh_data = cell(1, length(env_filtered_data));
for i = 2:length(pe_data)
    l_sig = Gfb_Analyzer_process(analyzer_aural, env_filtered_data{i}{1});
    r_sig = Gfb_Analyzer_process(analyzer_aural, env_filtered_data{i}{2});
    enh_data{i} = {l_sig, r_sig}; 
end

synthesizer_aural = Gfb_Synthesizer_new (analyzer_aural, 1/100);

for i = 2:length(pe_data)
    [l_sig, ~] = Gfb_Synthesizer_process(synthesizer_aural, env_filtered_data{i}{1});
    [r_sig, ~] = Gfb_Synthesizer_process(synthesizer_aural, env_filtered_data{i}{2});
    enh_data{i} = {l_sig, r_sig}; 
end


%% plot the figure
data_label = ["clean speech", "low reverberant speech"];
% data_label = ["clean speech", "low reverberant speech", "high reverberant speech"];

% figure
% for i = 1:length(env_lp_kurt_data)
%     for j = 1:size(env_data{1}{1},1)
%         subplot(12,3,i+(j-1)*3)
%         plot(env_lp_kurt_data{i}{1}(j,:))
%     end
% end

%% plot fine structure of data
data_label = ["clean speech", "low reverberant speech"];
%data_label = ["clean speech", "high reverberant speech"];
figure
for i = 1:2%linspace(2,1,2)
    for j = 1:6%size(env_data{1}{1},1)
        subplot(6,1,j)
        if i == 1
            plot(t(1:length(fine_data{i}{1}(j,:))), fine_data{i}{1}(j,:), "Color", [0 0 0 1])
        else
            plot(t(1:length(fine_data{i}{1}(j,:))), fine_data{i}{1}(j,:), "Color", [0.5 0.5 0.5 0.8])
        end
        legend(data_label)
        hold on
        ylabel("Subband Amplitude")
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
        ylabel("Subband Amplitude")
        if j == 6
            xlabel("Time (s)")
        end
    end
end

%% plot envelope of data
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
        ylabel("Envelope Amplitude")
        if j == 6
            xlabel("Sample")
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
        ylabel("Envelope Amplitude")
        if j == 6
            xlabel("Sample")
        end
    end
end

%% plot LP-analysis results
%data_label = ["clean speech", "low reverberant speech"];
data_label = ["clean speech", "high reverberant speech"];
figure
for i = 1:2:3%length(env_lp_data)
    for j = 1:6%size(env_data{1}{1},1)
        subplot(6,1,j)
        if i == 1
            plot(t(1:length(env_lp_lpc_data{i}{1}(j,:))), env_lp_lpc_data{i}{1}(j,:), "Color", [0 0 0 1])
        else
            plot(t(1:length(env_lp_lpc_data{i}{1}(j,:))), env_lp_lpc_data{i}{1}(j,:), "Color", [0.5 0.5 0.5 0.8])
        end
        legend(data_label)
        hold on
        ylabel("Envelope Amplitude")
        if j == 6
            xlabel("Sample")
        end
    end
end

figure
for i = 1:2:3%length(env_lp_data)
    for j = 1:6%size(env_data{1}{1},1)
        subplot(6,1,j)
        if i == 1
            plot(t(1:length(env_lp_lpc_data{i}{1}(j+6,:))), env_lp_lpc_data{i}{1}(j+6,:), "Color", [0 0 0 1])
        else
            plot(t(1:length(env_lp_lpc_data{i}{1}(j+6,:))), env_lp_lpc_data{i}{1}(j+6,:), "Color", [0.5 0.5 0.5 0.8])
        end
        legend(data_label)
        hold on
        ylabel("Envelope Amplitude")
        if j == 6
            xlabel("Sample")
        end
    end
end

%% plot kurtosis of LP-residual of envelope
%data_label = ["clean speech", "low reverberant speech"];
data_label = ["clean speech", "high reverberant speech"];
figure
for i = 1:2:3%length(env_lp_data)
    for j = 1:6%size(env_data{1}{1},1)
        subplot(6,1,j)
        if i == 1
            plot(t(1:length(env_lp_kurt_data{i}{1}(j,:))), env_lp_kurt_data{i}{1}(j,:), "Color", [0 0 0 1])
        else
            plot(t(1:length(env_lp_kurt_data{i}{1}(j,:))), env_lp_kurt_data{i}{1}(j,:), "Color", [0.5 0.5 0.5 0.8])
        end
        legend(data_label)
        hold on
        ylabel("Envelope Amplitude")
        if j == 6
            xlabel("Sample")
        end
    end
end

figure
for i = 1:2:3%length(env_lp_data)
    for j = 1:6%size(env_data{1}{1},1)
        subplot(6,1,j)
        if i == 1
            plot(t(1:length(env_lp_kurt_data{i}{1}(j+6,:))), env_lp_kurt_data{i}{1}(j+6,:), "Color", [0 0 0 1])
        else
            plot(t(1:length(env_lp_kurt_data{i}{1}(j+6,:))), env_lp_kurt_data{i}{1}(j+6,:), "Color", [0.5 0.5 0.5 0.8])
        end
        legend(data_label)
        hold on
        ylabel("Envelope Amplitude")
        if j == 6
            xlabel("Sample")
        end
    end
end

%%
data_label = ["clean speech", "low reverberant speech"];
%data_label = ["clean speech", "high reverberant speech"];

chan = 1;

figure
for i = 1%length(env_lp_data)
    for j = 1:6%size(env_data{1}{1},1)
        subplot(6,1,j)
        if i == 1
            plot(t(1:length(env_lp_lpc_kurt_data{i}{chan}(j,:))), env_lp_lpc_kurt_data{i}{chan}(j,:), "Color", [0 0 0 1])
        else
            plot(t(1:length(env_lp_lpc_kurt_data{i}{chan}(j,:))), env_lp_lpc_kurt_data{i}{chan}(j,:), "Color", [0.5 0.5 0.5 0.8])
        end
        legend(data_label)
        hold on
        ylabel("Envelope Amplitude")
        if j == 6
            xlabel("Sample")
        end
    end
end

figure
for i = 1%length(env_lp_data)
    for j = 1:6%size(env_data{1}{1},1)
        subplot(6,1,j)
        if i == 1
            plot(t(1:length(env_lp_lpc_kurt_data{i}{chan}(j+6,:))), env_lp_lpc_kurt_data{i}{chan}(j+6,:), "Color", [0 0 0 1])
        else
            plot(t(1:length(env_lp_lpc_kurt_data{i}{chan}(j+6,:))), env_lp_lpc_kurt_data{i}{chan}(j+6,:), "Color", [0.5 0.5 0.5 0.8])
        end
        legend(data_label)
        hold on
        ylabel("Envelope Amplitude")
        if j == 6
            xlabel("Sample")
        end
    end
end

%%

figure
for i = 1:length(env_filtered_data)
    for j = 1:4%size(env_data{1}{1},1)
        subplot(4,1,j)
        if i == 1
            plot(env_lp_data{i}{1}(j,:))
        else
            plot(env_filtered_data{i}{1}(j,:))
        end
        legend(data_label)
        hold on
        ylabel("Enhanced Envelope Amplitude")
        if j == 4
            xlabel("Sample")
        end
    end
end

figure
for i = 1:length(env_ss_data)
    for j = 1:4%size(env_data{1}{1},1)
        subplot(4,1,j)
        if i == 1
            plot(env_lp_data{i}{1}(j,:))
        else
            plot(abs(env_ss_data{i}{1}(j,:)))
        end
        legend(data_label)
        hold on
        ylabel("Enhanced Envelope Amplitude")
        if j == 4
            xlabel("Sample")
        end
    end
end