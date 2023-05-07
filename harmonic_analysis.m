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

idx = 10;

idx_start_short = 1 + idx*shortFrame;
idx_end_short = idx_start_short + shortFrame - 1;

idx_start_long = 1 + idx*longFrame;
idx_end_long = idx_start_long + longFrame - 1;

figure(4)
label = ["anechoic", "low reverberant", "high reverberant"];
for i = 1:size(voiced_left_data,1)
    l_voiced = voiced_left_data(i,:);
    r_voiced = voiced_right_data(i,:);

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

% estimate F0 of voiced audio sample

