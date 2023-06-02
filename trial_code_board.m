clc
clear all
close all

% Load the speech signal (replace 'speech.wav' with the path to your audio file)
[speech, fs] = audioread('sample_anechoic_min90.wav');
speech = speech(:,1);

start = find(speech, 1, "first");
last = find(speech, 1, "last");

speech = speech(start:last);

% Define frame duration and overlap
frameDuration = 0.032;  % 32ms
frameShift = 0.016;  % 16ms
frameLength = round(frameDuration * fs);
frameShiftLength = round(frameShift * fs);
overlap = frameLength - frameShiftLength;

% Calculate the number of frames
numFrames = floor((length(speech) - overlap) / frameShiftLength);

% Initialize the array for storing kurtosis values
kurtosisVals = zeros(numFrames, 1);

% Calculate kurtosis for each frame
for i = 1:numFrames
    startIndex = (i-1) * (frameLength - overlap) + 1;
    endIndex = startIndex + frameLength - 1;
    frame = speech(startIndex:endIndex);
    
    % Check for zero variance
    if var(frame) > 0
        kurtosisVals(i) = kurtosis(frame);
    else
        kurtosisVals(i) = 0;  % Set kurtosis to zero for frames with zero variance
    end
end

% Calculate the duration of the speech signal
signalDuration = length(speech) / fs;

% Create time vector for plotting
frameTimes = linspace(0, signalDuration, numFrames);

% Plot the kurtosis values
figure;
plot(frameTimes, kurtosisVals);
xlabel('Time (s)');
ylabel('Kurtosis');
title('Kurtosis of Speech Signal');

%% Load speech signal
% Load the speech signal (replace 'speech.wav' with the path to your audio file)
[speech, Fs] = audioread('sample_anechoic_min90.wav');
speech = speech(:,1);

start = find(speech, 1, "first");
last = find(speech, 1, "last");

speech = speech(start:last);

% Preprocessing
speech = speech/max(abs(speech)); % Normalize the speech signal

% LP analysis
order = 12; % LP order
windowLength = 20e-3; % Analysis window length in seconds
windowLengthSamples = round(windowLength*Fs); % Convert window length to samples
overlap = 0.5; % Overlap percentage between windows
overlapSamples = round(windowLengthSamples*overlap); % Calculate overlap samples

residual = []; % LP residual signal

% Perform LP analysis on overlapping windows
for i = 1:overlapSamples:(length(speech)-windowLengthSamples)
    window = speech(i:i+windowLengthSamples-1); % Extract current window
    [A, ~] = lpc(window, order); % Calculate LP coefficients using Levinson-Durbin algorithm
    currentResidual = filter(A, 1, window); % Obtain the LP residual signal for the current window
    residual = [residual; currentResidual]; % Append current residual to the overall residual signal
end

% Plot the LP residual
time = (0:length(residual)-1)/Fs; % Time axis in seconds
figure;
plot(time, residual);
xlabel('Time (s)');
ylabel('Amplitude');
title('LP Residual');

% Play the LP residual
soundsc(residual, Fs);