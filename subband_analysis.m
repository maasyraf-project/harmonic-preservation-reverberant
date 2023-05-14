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
par.gamma_order_auralisation = 3;
par.center_frequencies_hz_stimulation = [120 235 384 579 836 1175 1624 2222 3019 4084 5507 7410];
par.center_frequencies_hz_auralisation = [357 548 689 968 1483 2228 3319 4670 6630 9758 12530 15374];
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

env_lp_lpc_data = cell(1, length(env_data));
sig_vad = cell(1, length(raw_data));
for i = 1:length(env_data)
    l_sub_buf = nan(size(env_data{i}{1}));
    r_sub_buf = nan(size(env_data{i}{2}));

    % split audio into several frames
    lenFrames = 400;
    lenOverlap = 200;
    nFrames = ceil(size(env_data{i}{1},2)/lenFrames);

    l_vad_buf = nan(size(env_data{i}{1},1),nFrames);
    r_vad_buf = nan(size(env_data{i}{2},1),nFrames);

    for j = 1:size(env_data{1}{1},1)
        l = env_data{i}{1}(j,:);
        r = env_data{i}{2}(j,:);

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
            l_frame = l(idxStart:idxEnd);
            r_frame = r(idxStart:idxEnd);

            if or(length(l_frame) < lenOverlap, length(r_frame) < lenOverlap)
                l_frame = [l_frame zeros(1, lenOverlap - length(l_frame))];
                r_frame = [r_frame zeros(1, lenOverlap - length(r_frame))];
            end
            
            % calculate VAD
            [lv,VADPar] = VAD(l_frame', VADPar);
            [rv,VADPar] = VAD(r_frame', VADPar);

            l_vad = [l_vad, lv];
            r_vad = [r_vad, rv];
            
            % y_vad(idxStart:idxEnd) = lv;

            if lv == 1 
                % apply LPC on voiced segments
                la = lpc(l_frame, 12);
                estl_frame = filter([0 -la(2:end)], 1, l_frame);
                estl(idxStart:idxEnd) = estl(idxStart:idxEnd) + estl_frame;
            else 
                estl(idxStart:idxEnd) = estl(idxStart:idxEnd) + l(idxStart:idxEnd);
            end

            if rv == 1 
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

        l_sub_buf(j,:) = l_sig;
        r_sub_buf(j,:) = r_sig;

    end
    
    sig_vad = {l_vad_buf, r_vad_buf};
    env_lp_lpc_data{i} = {l_sub_buf, r_sub_buf};

end

env_lp_kurt_data = cell(1, length(env_lp_data));
for i = 1:length(env_lp_data)
    l_sub_buf = nan(size(env_lp_data{i}{1}));
    r_sub_buf = nan(size(env_lp_data{i}{2}));

    for j = 1:size(env_data{1}{1},1)
        lseg = 400;
        theta = 1e-6;
        l_sig = kurt(env_lp_lpc_data{i}{1}(j,:), lseg, theta);
        r_sig = kurt(env_lp_lpc_data{i}{2}(j,:), lseg, theta);

        l_sub_buf(j,:) = l_sig;
        r_sub_buf(j,:) = r_sig;
    end
    env_lp_kurt_data{i} = {l_sub_buf, r_sub_buf};

end

% conduct adaptive filtering
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
            for k=1:ss:size(env_kurt_data{i}{1},2)
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

                %         z3y=fft(zr3).*cYrn;
                %         zy=fft(zrn).*cYrn;
                %
                %         gJ=4*(Z2*z3y-Z4*zy)/(Z2^3+1e-20)*ss;
                %         Gh=Gh+mu*gJ;


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

        disp("processing Right Channel")
        tic
        for m=1:niter
            yrn=zeros(bs,1);
            for k=1:ss:size(env_kurt_data{i}{2},2)
                yrn(1:p-1)=yrn(end-p+2:end);

                if k+ss-1 < size(env_lp_kurt_data{i}{2},2)
                    yrn(p:end) = env_lp_kurt_data{i}{2}(j,k:k+ss-1);
                else
                    env_padded = [env_lp_kurt_data{i}{2}(j,k:end) zeros(1, size(yrn(p:end),1)-size(env_lp_kurt_data{i}{1}(j,k:end),2))];
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

                %         z3y=fft(zr3).*cYrn;
                %         zy=fft(zrn).*cYrn;
                %
                %         gJ=4*(Z2*z3y-Z4*zy)/(Z2^3+1e-20)*ss;
                %         Gh=Gh+mu*gJ;


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

% conduct spectral substraction
env_ss_data = cell(1, length(env_data));
for i = 1:length(env_filtered_data)
    l_sub_buf = nan(size(env_lp_data{i}{1}));
    r_sub_buf = nan(size(env_lp_data{i}{2}));

    for j = 1:size(env_filtered_data{1}{1},1)
        % define parameter
        [L, f, t] = stft(env_filtered_data{i}{1}(j,:), 512, 256, 512, fs);
        [R, f, t] = stft(env_filtered_data{i}{2}(j,:), 512, 256, 512, fs);
        
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

        LP_att=(L_pow-L_sub)./L_pow;
        LP_att(LP_att<epS)=epS;

        RP_att=(R_pow-R_sub)./R_pow;
        RP_att(RP_att<epS)=epS;
        
        % predict enhanced specstrum
        SL_pow = L_pow .* LP_att;
        SR_pow = R_pow .* RP_att;

        SL_phase = sqrt(SL_pow).*exp(1i*L_phase);
        SR_phase = sqrt(SR_pow).*exp(1i*R_phase);

        % spectrum analysis
        l_sub_buf(j,:) = l_sig;
        r_sub_buf(j,:) = r_sig;
    end
    env_lp_kurt_data{i} = {l_sub_buf, r_sub_buf};

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
% analyzer_aural = Gfb_Analyzer_new(par.voc_sampling_frequency_hz,...
%     par.gamma_order_auralisation,...
%     par.center_frequencies_hz_auralisation,...
%     par.bandwidth_factor);
% 
% enh_data = cell(1, length(env_filtered_data));
% for i = 2:length(pe_data)
%     [l_sig, ~] = Gfb_Analyzer_process(analyzer_aural, env_filtered_data{i}{1});
%     [r_sig, ~] = Gfb_Analyzer_process(analyzer_aural, env_filtered_data{i}{2});
%     enh_data{i} = {l_sig, r_sig}; 
% end
% 
% synthesizer_aural = Gfb_Synthesizer_new (analyzer_aural, 1/100);
% 
% for i = 2:length(pe_data)
%     [l_sig, ~] = Gfb_Synthesizer_process(synthesizer_aural, env_filtered_data{i}{1});
%     [r_sig, ~] = Gfb_Synthesizer_process(synthesizer_aural, env_filtered_data{i}{2});
%     enh_data{i} = {l_sig, r_sig}; 
% end


%% plot the figure
data_label = ["clean speech", "low reverberant speech", "high reverberant speech"];

% figure
% for i = 1:length(env_lp_kurt_data)
%     for j = 1:size(env_data{1}{1},1)
%         subplot(12,3,i+(j-1)*3)
%         plot(env_lp_kurt_data{i}{1}(j,:))
%     end
% end

figure
for i = 1:length(env_lp_data)
    for j = 1:4%size(env_data{1}{1},1)
        subplot(4,1,j)
        plot(env_lp_data{i}{1}(j,:))
        legend(data_label)
        hold on
        ylabel("Envelope Amplitude")
        if j == 4
            xlabel("Sample")
        end
    end
end

figure
for i = 1:length(env_lp_lpc_data)
    for j = 1:4%size(env_data{1}{1},1)
        subplot(4,1,j)
        plot(env_lp_lpc_data{i}{1}(j,:))
        legend(data_label)
        hold on
        ylabel("LPC Amplitude")
        if j == 4
            xlabel("Sample")
        end
    end
end

figure
for i = 1:length(env_lp_kurt_data)
    for j = 1:4%size(env_data{1}{1},1)
        subplot(4,1,j)
        plot(env_lp_kurt_data{i}{1}(j,:))
        legend(data_label)
        hold on
        ylabel("Kurtosis Amplitude")
        if j == 4
            xlabel("Sample")
        end
    end
end

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