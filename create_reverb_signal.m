clc
clear all 
close all
addpath(genpath('additional-packages'));

% load mono-clean audio signal
[x, fs] = audioread("full-speech-sample.wav");

% load impulse response
irFilename = 'impulse_responses/surrey_cortex_rooms/UniS_Room_D_BRIR_16k.sofa';
irName = split(irFilename,'/');

irFilenameDir = char(strcat('\', irName(1),'\', irName(2),'\', irName(3)));
irDir = char(strcat(pwd, '\additional-packages\two-ears\BinauralSimulator\tmp',irFilenameDir));

irName = char(irName(end));
irName = irName(1:end-5);

if ~exist(irDir, "file")
    irFilename = db.downloadFile(irFilename);
    ir = SOFAload(irFilename);
else
    ir = SOFAload(strcat(irName,'.sofa'));
end

ir_left = resample(squeeze(ir.Data.IR(1, 1, :)), fs, ir.Data.SamplingRate);
ir_right = resample(squeeze(ir.Data.IR(1, 2, :)), fs, ir.Data.SamplingRate);

y = [conv(squeeze(ir_left), x) ...
    conv(squeeze(ir_right), x)];

audioOutputNames = 'sample_rev_high_min90.wav';
audiowrite(audioOutputNames, y, fs);