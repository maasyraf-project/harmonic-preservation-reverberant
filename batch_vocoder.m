clc
clear all
close all

% load all packages
addpath(genpath('additional-packages'));
addpath(genpath('stimuli'));
addpath(genpath('results'));

resultsDir = strcat(pwd,'\results\statistical_analysis\');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

% create directory path for new result
saveDir = strcat(pwd,'\results\room_statistical_analysis\');
if ~exist(saveDir , 'dir')
    mkdir(saveDir);
end

resultsFilenames = dir(resultsDir);
resultsFilenames = resultsFilenames(3:end);

numSentences = 5;
resultsFilenames = resultsFilenames(1:9*numSentences);

degree = ["0", "30", "45", "60", "90", "min30", "min45", "min60", "min90"];

fine_statistic = cell(1, length(degree)); % the structure is degree - sentence - room - channel - subband
env_statistic = cell(1, length(degree));

for d = 1:length(degree)

    FINE = cell(1, numSentences);
    %FINE_LP = cell(1, numSentences);
    ENV = cell(1, numSentences);
    %ENV_LP = cell(1, numSentences);

    idx = 1;
    for i = d:9:length(resultsFilenames) % indicator for sentences
        load(resultsFilenames(i).name)
        disp(strcat("Processing ", resultsFilenames(i).name));

        fine = cell(1, 5);
        %fine_lp = cell(1, 5);
        env = cell(1, 5);
        %env_lp = cell(1, 5);

        for j = 1:5 % indicator for numbers of room 

            

        end

        FINE{idx} = fine;
        %FINE_LP{idx} = fine_lp;
        ENV{idx} = env;
        %ENV_LP{idx} = env_lp;

        idx = idx + 1;

    end

    fine_statistic{d} = FINE;
    %fine_lp_statistic{d} = FINE_LP;
    env_statistic{d} = ENV;
    %env_lp_statistic{d} = ENV_LP;

    disp(strcat("Completing processing data for sound comes from ", degree(d)));

end