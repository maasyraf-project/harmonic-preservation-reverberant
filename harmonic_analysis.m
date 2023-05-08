clc
clear all
close all
addpath(genpath('additional-packages'));

%% import and auralize full audio
[h, fs] = audioread("sample_rev_high_min90.wav");
hl = h(:,1);
hr = h(:,2);
t = (0:size(h,1)-1)/fs;

[l, fs] = audioread("sample_rev_low_min90.wav");
l = [l' zeros(2, length(h)-length(l))]';
ll = l(:,1);
lr = l(:,2);

[x, fs] = audioread("sample_anechoic_min90.wav");
x = [x' zeros(2, length(h)-length(x))]';
xl = x(:,1);
xr = x(:,2);

signal_data = {x'; l'; h'};
raw_data = [xl'; xr'; ll'; lr'; hl'; hr'];

figure(1)
for i = 1:size(raw_data,1)
    subplot(size(raw_data,1)/2,2,i)
    plot(t, raw_data(i,:))
    xlim([0 3.2])
    ylim([-0.15 0.15])
    hold on
end

figure(2)
for i = 1:size(raw_data,1)
    subplot(size(raw_data,1)/2,2,i)
    spectrogram(raw_data(i,:), 512, 256, 512, fs, 'yaxis')
    xlim([0 3.2])
    set(colorbar, 'visible', 'on')
end

%% speech analysis
% conduct VAD on clean speech
Nw = 512;
Nsh = 256;
yl = 0.05*G729(xl, fs, Nw, Nsh);
yr = 0.05*G729(xr, fs, Nw, Nsh);

figure(1)
for i = 1:size(raw_data,1)
    subplot(size(raw_data,1)/2,2,i)
    if mod(i,2) == 0
        plot(t, yr)
    else
        plot(t, yl)
    end
    xlim([0 3.2])
    ylim([-0.15 0.15])
    hold on
end

% select a voiced segment from clean and reverberant speech based on VAD
% results on clean speech
xl_voiced = xl(yl~=0);
xr_voiced = xr(yr~=0);
ll_voiced = ll(yl~=0);
lr_voiced = lr(yr~=0);
hl_voiced = hl(yl~=0);
hr_voiced = hr(yr~=0);

voiced_left_data  = [xl_voiced'; ll_voiced'; hl_voiced'];
voiced_right_data  = [xr_voiced'; lr_voiced'; hr_voiced'];

shortFrame = 30e-3 * fs;
longFrame = 80e-3 * fs;

idx = 20;

idx_start_short = 1 + idx*shortFrame;
idx_end_short = idx_start_short + shortFrame - 1;

idx_start_long = 1 + idx*longFrame;
idx_end_long = idx_start_long + longFrame - 1;

figure(3)
label = ["anechoic", "low reverberant", "high reverberant"];
for i = 1:size(voiced_left_data,1)
    l_voiced = voiced_left_data(i,:);
    r_voiced = voiced_right_data(i,:);

    l_voiced_shortFrame = l_voiced(idx_start_short:idx_end_short);
    r_voiced_shortFrame = r_voiced(idx_start_short:idx_end_short);

    l_voiced_longFrame = l_voiced(idx_start_long:idx_end_long);
    r_voiced_longFrame = r_voiced(idx_start_long:idx_end_long);

    % make fft observation on short frame
    nfft_short = 2^nextpow2(shortFrame);
    F_short = fs/2 * linspace(0,1,shortFrame/2+1);
    win_short = hamming(shortFrame);

    L_short = (1/shortFrame) * fft([zeros(1,(nfft_short-shortFrame)/2) l_voiced_shortFrame zeros(1,(nfft_short-shortFrame)/2)], nfft_short);
    R_short = (1/shortFrame) * fft([zeros(1,(nfft_short-shortFrame)/2) r_voiced_shortFrame zeros(1,(nfft_short-shortFrame)/2)], nfft_short);
    
    L_short = abs(L_short(1:shortFrame/2+1));
    R_short = abs(R_short(1:shortFrame/2+1));

    % make fft observation on long frame
    nfft_long = 2^nextpow2(longFrame);
    F_long = fs/2 * linspace(0,1,longFrame/2+1);
    win_long = hamming(longFrame);

    L_long = (1/longFrame) * fft([l_voiced_longFrame zeros(1,nfft_short-shortFrame)], nfft_long);
    R_long = (1/longFrame) * fft([r_voiced_longFrame zeros(1,nfft_short-shortFrame)], nfft_long);
    
    L_long = abs(L_long(1:longFrame/2+1));
    R_long = abs(R_long(1:longFrame/2+1));
    
    % plot the result
    subplot(2,2,1)
    plot(F_short, L_short)
    xlim([0 3000])
    legend(label)
    hold on

    subplot(2,2,2)
    plot(F_short, R_short)
    xlim([0 3000])
    legend(label)
    hold on

    subplot(2,2,3)
    plot(F_long, L_long)
    xlim([0 3000])
    legend(label)
    hold on

    subplot(2,2,4)
    plot(F_long, R_long)
    xlim([0 3000])
    legend(label)
    hold on

end

% estimate F0 of voiced audio sample

%% Speech enhancement implementation
XL = stft(xl, fs, "Window",hanning(512),"OverlapLength",256, "FFTLength", 512, "FrequencyRange","onesided");

XR = stft(xr, fs, "Window",hanning(512),"OverlapLength",256, "FFTLength", 512, "FrequencyRange","onesided");

LL = stft(ll, fs, "Window",hanning(512),"OverlapLength",256, "FFTLength", 512, "FrequencyRange","onesided");

LR = stft(lr, fs, "Window",hanning(512),"OverlapLength",256, "FFTLength", 512, "FrequencyRange","onesided");

HL = stft(hl, fs, "Window",hanning(512),"OverlapLength",256, "FFTLength", 512, "FrequencyRange","onesided");

HR = stft(hr, fs, "Window",hanning(512),"OverlapLength",256, "FFTLength", 512, "FrequencyRange","onesided");

% low-reverberant enhancement
IRM_LL = XL.^2 ./ (LL.^2 + XL.^2 + eps);
IRM_LR = XR.^2 ./ (LR.^2 + XR.^2 + eps);

SLL = LL .* IRM_LL;
SLR = LR .* IRM_LR;

sll = istft(SLL, fs, "Window", hanning(512), "OverlapLength", 256, "FFTLength", 512, "FrequencyRange", "onesided");
slr = istft(SLR, fs, "Window", hanning(512), "OverlapLength", 256, "FFTLength", 512, "FrequencyRange", "onesided");

% high-reverberant enhancement
IRM_HL = XL.^2 ./ (HL.^2 + XL.^2 + eps);
IRM_HR = XR.^2 ./ (HR.^2 + XR.^2 + eps);

SHL = HL .* IRM_HL;
SHR = HR .* IRM_HR;

shl = istft(SHL, fs, "Window", hanning(512), "OverlapLength", 256, "FFTLength", 512, "FrequencyRange", "onesided");
shr = istft(SHR, fs, "Window", hanning(512), "OverlapLength", 256, "FFTLength", 512, "FrequencyRange", "onesided");

enhanced_signal_data = {[xl, xr]', [sll, slr]', [shl, shr]'};

enhanced_data = [xl(1:length(sll))';xr(1:length(sll))';sll';slr';shl';shr'];
te = (0:size(enhanced_data,2)-1)/fs;

% plot enhanced data
figure(1)
for i = 1:size(enhanced_data,1)
    subplot(size(enhanced_data,1)/2,2,i)
    plot(te, enhanced_data(i,:))
    xlim([0 3.2])
    ylim([-0.15 0.15])
    hold on
end

figure(4)
for i = 1:size(enhanced_data,1)
    subplot(size(enhanced_data,1)/2,2,i)
    spectrogram(enhanced_data(i,:), 512, 256, 512, fs, 'yaxis')
    xlim([0 3.2])
    set(colorbar, 'visible', 'on')
end


%% enhanced speech analysis
% select a voiced segment from clean and enhanced speech based on VAD
% results on clean speech
sll_voiced = sll(yl(1:length(sll))~=0);
slr_voiced = slr(yr(1:length(sll))~=0);
shl_voiced = shl(yl(1:length(sll))~=0);
shr_voiced = shr(yr(1:length(sll))~=0);

enhLow_voiced_left_data  = [xl_voiced'; ll_voiced'; sll_voiced'];
enhLow_voiced_right_data  = [xr_voiced'; lr_voiced'; slr_voiced'];

enhHigh_voiced_left_data  = [xl_voiced'; hl_voiced'; shl_voiced'];
enhHigh_voiced_right_data  = [xr_voiced'; hr_voiced'; shr_voiced'];

shortFrame = 30e-3 * fs;
longFrame = 80e-3 * fs;

idx = 20;

idx_start_short = 1 + idx*shortFrame;
idx_end_short = idx_start_short + shortFrame - 1;

idx_start_long = 1 + idx*longFrame;
idx_end_long = idx_start_long + longFrame - 1;

figure(5)
label = ["anechoic", "low reverberant", "enhanced low reverberant"];
for i = 1:size(enhLow_voiced_left_data,1)
    l_voiced = enhLow_voiced_left_data(i,:);
    r_voiced = enhLow_voiced_right_data(i,:);

    l_voiced_shortFrame = l_voiced(idx_start_short:idx_end_short);
    r_voiced_shortFrame = r_voiced(idx_start_short:idx_end_short);

    l_voiced_longFrame = l_voiced(idx_start_long:idx_end_long);
    r_voiced_longFrame = r_voiced(idx_start_long:idx_end_long);

    % make fft observation on short frame
    F_short = fs/2 * linspace(0,1,shortFrame/2+1);
    L_short = fft(l_voiced_shortFrame, shortFrame);
    R_short = fft(r_voiced_shortFrame, shortFrame);
    
    L_short = abs(L_short(1:shortFrame/2+1));
    R_short = abs(R_short(1:shortFrame/2+1));

    % make fft observation on long frame
    F_long = fs/2 * linspace(0,1,longFrame/2+1);
    L_long = fft(l_voiced_longFrame, longFrame);
    R_long = fft(r_voiced_longFrame, longFrame);
    
    L_long = abs(L_long(1:longFrame/2+1));
    R_long = abs(R_long(1:longFrame/2+1));
    
    % plot the result
    subplot(2,2,1)
    plot(F_short, L_short)
    xlim([0 3000])
    legend(label)
    hold on

    subplot(2,2,2)
    plot(F_short, R_short)
    xlim([0 3000])
    legend(label)
    hold on

    subplot(2,2,3)
    plot(F_long, L_long)
    xlim([0 3000])
    legend(label)
    hold on

    subplot(2,2,4)
    plot(F_long, R_long)
    xlim([0 3000])
    legend(label)
    hold on

end

figure(6)
label = ["anechoic", "high reverberant", "enhanced high reverberant"];
for i = 1:size(enhHigh_voiced_left_data,1)
    l_voiced = enhHigh_voiced_left_data(i,:);
    r_voiced = enhHigh_voiced_right_data(i,:);

    l_voiced_shortFrame = l_voiced(idx_start_short:idx_end_short);
    r_voiced_shortFrame = r_voiced(idx_start_short:idx_end_short);

    l_voiced_longFrame = l_voiced(idx_start_long:idx_end_long);
    r_voiced_longFrame = r_voiced(idx_start_long:idx_end_long);

    % make fft observation on short frame
    F_short = fs/2 * linspace(0,1,shortFrame/2+1);
    L_short = fft(l_voiced_shortFrame, shortFrame);
    R_short = fft(r_voiced_shortFrame, shortFrame);
    
    L_short = abs(L_short(1:shortFrame/2+1));
    R_short = abs(R_short(1:shortFrame/2+1));

    % make fft observation on long frame
    F_long = fs/2 * linspace(0,1,longFrame/2+1);
    L_long = fft(l_voiced_longFrame, longFrame);
    R_long = fft(r_voiced_longFrame, longFrame);
    
    L_long = abs(L_long(1:longFrame/2+1));
    R_long = abs(R_long(1:longFrame/2+1));
    
    % plot the result
    subplot(2,2,1)
    plot(F_short, L_short)
    xlim([0 3000])
    legend(label)
    hold on

    subplot(2,2,2)
    plot(F_short, R_short)
    xlim([0 3000])
    legend(label)
    hold on

    subplot(2,2,3)
    plot(F_long, L_long)
    xlim([0 3000])
    legend(label)
    hold on

    subplot(2,2,4)
    plot(F_long, R_long)
    xlim([0 3000])
    legend(label)
    hold on

end

%% spatial cues analysis
% Parameters of the auditory filterbank processor
fb_type       = 'gammatone';
fb_lowFreqHz  = 80;
fb_highFreqHz = 8000;
fb_nChannels  = 512;  

% Parameters of innerhaircell processor
ihc_method    = 'dau';

% Parameters of crosscorrelation processor (for obtain ITD)
cc_wSizeSec   = 0.02;
cc_hSizeSec   = 0.01;
cc_wname      = 'hann';

% Parameters of ILD processor
ild_wSizeSec  = 0.02;
ild_hSizeSec  = 0.01;
ild_wname     = 'hann';

% Summary of parameters 
itd_par = genParStruct('fb_type',fb_type,'fb_lowFreqHz',fb_lowFreqHz,...
                   'fb_highFreqHz',fb_highFreqHz,'fb_nChannels',fb_nChannels,...
                   'ihc_method',ihc_method,'cc_wSizeSec',cc_wSizeSec,...
                   'cc_hSizeSec',cc_hSizeSec,'cc_wname',cc_wname);

ild_par = genParStruct('fb_type',fb_type,'fb_lowFreqHz',fb_lowFreqHz,...
                   'fb_highFreqHz',fb_highFreqHz,'fb_nChannels',fb_nChannels,...
                   'ihc_method',ihc_method,'ild_wSizeSec',ild_wSizeSec,...
                   'cc_hSizeSec',ild_hSizeSec,'cc_wname',ild_wname);

% raw signal
for i = 1%1:size(signal_data,1)
    x = signal_data{i}';

    %% ITD Calculation
    % Create a data object based on the ear signals
    itd_dObj = dataObject(x(1:end,:),fs);

    % Request interaural time differences (ITDs)
    itd_requests = {'itd'};

    % Create a manager
    itd_mObj = manager(itd_dObj,itd_requests,itd_par);

    % Request processing
    itd_mObj.processSignal();

    % Plot results
    itd_dObj.plot([],[],'bGray',1,'decimateRatio',3,'bSignal',1);
    ylim([-1.25 1.25]);

    % Plot ITDs
    itd_dObj.itd{1}.plot;

    %% ILD Calculation
    % Create a data object based on the ear signals
    ild_dObj = dataObject(x(1:end,:),fs);

    % Request interaural level differences (ILDs)
    ild_requests = {'ild'};

    % Create a manager
    ild_mObj = manager(ild_dObj,ild_requests,ild_par);

    % Request processing
    ild_mObj.processSignal();

    % Plot ILDs
    ild_dObj.ild{1}.plot;

end

% enhanced signal
for i = 1%1:size(enhanced_signal_data,1)
    x = enhanced_signal_data{i}';

    %% ITD Calculation
    % Create a data object based on the ear signals
    itd_dObj = dataObject(x(1:end,:),fs);

    % Request interaural time differences (ITDs)
    itd_requests = {'itd'};

    % Create a manager
    itd_mObj = manager(itd_dObj,itd_requests,itd_par);

    % Request processing
    itd_mObj.processSignal();

    % Plot results
    itd_dObj.plot([],[],'bGray',1,'decimateRatio',3,'bSignal',1);
    ylim([-1.25 1.25]);

    % Plot ITDs
    itd_dObj.itd{1}.plot;

    %% ILD Calculation
    % Create a data object based on the ear signals
    ild_dObj = dataObject(x(1:end,:),fs);

    % Request interaural level differences (ILDs)
    ild_requests = {'ild'};

    % Create a manager
    ild_mObj = manager(ild_dObj,ild_requests,ild_par);

    % Request processing
    ild_mObj.processSignal();

    % Plot ILDs
    ild_dObj.ild{1}.plot;

end

% spatial cues calculation
clean_itd = estimate_ITD_Broadband(x, fs)*1000;
rev_low_itd = estimate_ITD_Broadband(l, fs)*1000;
rev_high_itd = estimate_ITD_Broadband(h, fs)*1000;
enh_low_itd = estimate_ITD_Broadband([sll, slr], fs)*1000;
enh_high_itd = estimate_ITD_Broadband([shl, shr], fs)*1000;

clean_ild = 20*log10(rms(xr)/rms(xl));
rev_low_ild = 20*log10(rms(lr)/rms(ll));
rev_high_ild = 20*log10(rms(hr)/rms(hl));
enh_low_ild = 20*log10(rms(slr)/rms(sll));
enh_high_ild = 20*log10(rms(shr)/rms(shl));

clean_sii = mbstoi(xl, xr, xl, xr, fs);
rev_low_sii = mbstoi(xl, xr, ll, lr, fs);
rev_high_sii = mbstoi(xl, xr, hl, hr, fs);
enh_low_sii = mbstoi(xl, xr, sll, slr, fs);
enh_high_sii = mbstoi(xl, xr, shl, shr, fs);