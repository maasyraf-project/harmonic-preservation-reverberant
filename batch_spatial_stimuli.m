clc
clear all
close all

% load all packages
addpath(genpath('additional-packages'));
target_fs = 16e3;

startTwoEars  % optional

% load impulse response
irFilename = 'impulse_responses/surrey_cortex_rooms/UniS_Anechoic_BRIR_16k.sofa';
irName = split(irFilename,'/');

irFilenameDir = char(strcat('\', irName(1),'\', irName(2),'\', irName(3)));
irDir = char(strcat(pwd, '\additional-packages\TwoEars\BinauralSimulator\tmp',irFilenameDir));

irName = char(irName(end));
irName = irName(1:end-5);

if ~exist(irDir, "file")
    irFilename = db.downloadFile(irFilename);
    ir = SOFAload(irFilename);
else
    ir = SOFAload(strcat(irName,'.sofa'));
end

% link to clean stimuli directory
audioInputDir = strcat(pwd,'\stimuli\clean');
audioInputNames = dir(fullfile(audioInputDir, '*.wav'));

% create directory path for spatialized reverberant stimmuli
audioOutputDir = strcat(pwd,'\stimuli\unvocoded\',irName);
if ~exist(audioOutputDir, 'dir')
    mkdir(audioOutputDir);
end

% directory path for reference stimuli
refDir = strcat(pwd,'\stimuli\unvocoded\UniS_Anechoic_BRIR_16k');
if ~exist(refDir, 'dir')
    mkdir(refDir);
end

% create directory path for result
resultsDir = strcat(pwd,'\results\unvocoded\');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

% initiate metric for all stimuli
itd = [];
ild = [];
sii = [];

% spatialize the stimuli
for i = 1:length(audioInputNames)
    % import audio
    [audioInput, fs] = audioread(fullfile(audioInputDir, audioInputNames(i).name));
    
    % normalize the input audio
    audioInput = audioInput ./ max(abs(audioInput));
    
    % resample audio
    if fs == target_fs
        audioInput = audioInput;
    else
        [P, Q] = rat(target_fs/fs);
        audioInput = resample(audioInput, P, Q);
    end
    
    % initiate metric for a stimuli for all degree
    itdDegree = [];
    ildDegree = [];
    siiDegree = [];
    
    for j = 1:2:size(ir.Data.IR,1)
        % make output audio filenames
%         if j < 19
%             OutputFilenames = strcat(audioInputNames(i).name(1:end-4),'_', 'min', string(abs(-(ir.SourceView(j,1) - 180))),'.wav');
%         else
%             OutputFilenames = strcat(audioInputNames(i).name(1:end-4),'_', string(-(ir.SourceView(j,1) - 180)),'.wav');
%         end
        
        % make output audio filenames (spesific for anechoic surrey room
        % SOFA format)
        if j < 19
            OutputFilenames = strcat(audioInputNames(i).name(1:end-4),'_','min', string(abs((ir.SourcePosition(j,1) - 360))),'.wav');
        else
            OutputFilenames = strcat(audioInputNames(i).name(1:end-4),'_',string((ir.SourcePosition(j,1))),'.wav');
        end

        audioOutputNames = fullfile(audioOutputDir, OutputFilenames);
        
        if ir.Data.SamplingRate ~= target_fs
            ir_left = resample(squeeze(ir.Data.IR(j, 1, :)), target_fs, ir.Data.SamplingRate);
            ir_right = resample(squeeze(ir.Data.IR(j, 2, :)), target_fs, ir.Data.SamplingRate);
        else
            ir_left = squeeze(ir.Data.IR(j, 1, :));
            ir_right = squeeze(ir.Data.IR(j, 2, :));
        end

        % make binaural stimuli
        audioOutput = [conv(squeeze(ir_left), audioInput) ...
               conv(squeeze(ir_right), audioInput)];
        
        % avoid clipping on audio output
        audioOutput(audioOutput > 1) = 1;
        audioOutput(audioOutput < -1) = -1;
        
        % write output audio file
        audiowrite(audioOutputNames, audioOutput, target_fs);
        
        % make left channel stimuli
        audioLeft = audioOutput(:,1);
        
        % make right channel stimuli
        audioRight = audioOutput(:,2);
        
        % input reference audio
        [ref, fs] = audioread(fullfile(refDir, OutputFilenames));
        refLeft = ref(:,1);
        refRight = ref(:,2);

        % calculate metric
        itdValue = estimate_ITD_Broadband(audioOutput, target_fs)*1000;             % in ms
        ildValue = 20*log10(rms(audioRight)/rms(audioLeft));                        % in dB
        siiValue = mbstoi(refLeft, refRight, audioLeft, audioRight, target_fs);     % in range of 0 to 1

        
        itdDegree = [itdDegree, itdValue];
        ildDegree = [ildDegree, ildValue];
        siiDegree = [siiDegree, siiValue];

        disp(strcat(OutputFilenames, ' is created, ILD: ', string(ildValue), ' - ITD: ', string(itdValue), ' - SII: ', string(siiValue)))
        
    end
    
    % recap all calculated metric for all stimuli
    itd = [itd; itdDegree];
    ild = [ild; ildDegree];
    sii = [sii; siiDegree];
    
end

% save the calculated output
outputMetricFilename = strcat(resultsDir, irName, '_metric.mat');
save(outputMetricFilename, 'itd', 'ild', 'sii');