clc
clear all
close all

% load all packages
addpath(genpath('additional-packages'));
addpath(genpath('stimuli'));
addpath(genpath('results'));

resultsDir = strcat(pwd,'\results\room_statistical_analysis\');
roomDataName = strcat(pwd,'\results\room_statistical_analysis\room_data.mat');

load(roomDataName)

degree = ["0", "30", "45", "60", "90", "min30", "min45", "min60", "min90"];
room = ["Anechoic", "Room A", "Room B", "Room C", "Room D"];

%% obtain mean data
cm = 1;
mean_data = cell(1, length(degree));
for d = 1:length(degree)
    lroom_data = zeros(12, 5);
    rroom_data = zeros(12, 5);
    
    for ro = 1:5
        lsentence_data = zeros(12, 5);
        rsentence_data = zeros(12, 5);

        for s = 1:5 
            lsubband_data = cell(12, 1);
            rsubband_data = cell(12, 1);

            for j = 1:12
                l = fine_statistic{d}{s}{ro}{1}{j}(cm);
                r = fine_statistic{d}{s}{ro}{2}{j}(cm);

                lsubband_data{j} = l;
                rsubband_data{j} = r;

            end

            lsentence_data(:,s) = cell2mat(lsubband_data);
            rsentence_data(:,s) = cell2mat(rsubband_data);
        
        end
        % average values across all used sentences
        lsentence_data = mean(lsentence_data, 2);
        rsentence_data = mean(rsentence_data, 2);
        
        lroom_data(:,ro) = lsentence_data;
        rroom_data(:,ro) = rsentence_data;

    end
    mean_data{d} = {lroom_data, rroom_data};
    
end

%% obtain variance data
cm = 2;
variance_data = cell(1, length(degree));
for d = 1:length(degree)
    lroom_data = zeros(12, 5);
    rroom_data = zeros(12, 5);
    
    for ro = 1:5
        lsentence_data = zeros(12, 5);
        rsentence_data = zeros(12, 5);

        for s = 1:5 
            lsubband_data = cell(12, 1);
            rsubband_data = cell(12, 1);

            for j = 1:12
                l = fine_statistic{d}{s}{ro}{1}{j}(cm);
                r = fine_statistic{d}{s}{ro}{2}{j}(cm);

                lsubband_data{j} = l;
                rsubband_data{j} = r;

            end

            lsentence_data(:,s) = cell2mat(lsubband_data);
            rsentence_data(:,s) = cell2mat(rsubband_data);
        
        end
        % average values across all used sentences
        lsentence_data = mean(lsentence_data, 2);
        rsentence_data = mean(rsentence_data, 2);
        
        lroom_data(:,ro) = lsentence_data;
        rroom_data(:,ro) = rsentence_data;

    end
    variance_data{d} = {lroom_data, rroom_data};
    
end

%% obtain variance data
cm = 3;
skewness_data = cell(1, length(degree));
for d = 1:length(degree)
    lroom_data = zeros(12, 5);
    rroom_data = zeros(12, 5);
    
    for ro = 1:5
        lsentence_data = zeros(12, 5);
        rsentence_data = zeros(12, 5);

        for s = 1:5 
            lsubband_data = cell(12, 1);
            rsubband_data = cell(12, 1);

            for j = 1:12
                l = fine_statistic{d}{s}{ro}{1}{j}(cm);
                r = fine_statistic{d}{s}{ro}{2}{j}(cm);

                lsubband_data{j} = l;
                rsubband_data{j} = r;

            end

            lsentence_data(:,s) = cell2mat(lsubband_data);
            rsentence_data(:,s) = cell2mat(rsubband_data);
        
        end
        % average values across all used sentences
        lsentence_data = mean(lsentence_data, 2);
        rsentence_data = mean(rsentence_data, 2);
        
        lroom_data(:,ro) = lsentence_data;
        rroom_data(:,ro) = rsentence_data;

    end
    skewness_data{d} = {lroom_data, rroom_data};
    
end

%% obtain kurtosis data
cm = 4;
kurtosis_data = cell(1, length(degree));
for d = 1:length(degree)
    lroom_data = zeros(12, 5);
    rroom_data = zeros(12, 5);
    
    for ro = 1:5
        lsentence_data = zeros(12, 5);
        rsentence_data = zeros(12, 5);

        for s = 1:5 
            lsubband_data = cell(12, 1);
            rsubband_data = cell(12, 1);

            for j = 1:12
                l = fine_statistic{d}{s}{ro}{1}{j}(cm);
                r = fine_statistic{d}{s}{ro}{2}{j}(cm);

                lsubband_data{j} = l;
                rsubband_data{j} = r;

            end

            lsentence_data(:,s) = cell2mat(lsubband_data);
            rsentence_data(:,s) = cell2mat(rsubband_data);
        
        end
        % average values across all used sentences
        lsentence_data = mean(lsentence_data, 2);
        rsentence_data = mean(rsentence_data, 2);
        
        lroom_data(:,ro) = lsentence_data;
        rroom_data(:,ro) = rsentence_data;

    end
    kurtosis_data{d} = {lroom_data, rroom_data};
    
end

%% obtain kewkurtosis data
cm = 5;
kewkurtosis_data = cell(1, length(degree));
for d = 1:length(degree)
    lroom_data = zeros(12, 5);
    rroom_data = zeros(12, 5);
    
    for ro = 1:5
        lsentence_data = zeros(12, 5);
        rsentence_data = zeros(12, 5);

        for s = 1:5 
            lsubband_data = cell(12, 1);
            rsubband_data = cell(12, 1);

            for j = 1:12
                l = fine_statistic{d}{s}{ro}{1}{j}(cm);
                r = fine_statistic{d}{s}{ro}{2}{j}(cm);

                lsubband_data{j} = l;
                rsubband_data{j} = r;

            end

            lsentence_data(:,s) = cell2mat(lsubband_data);
            rsentence_data(:,s) = cell2mat(rsubband_data);
        
        end
        % average values across all used sentences
        lsentence_data = mean(lsentence_data, 2);
        rsentence_data = mean(rsentence_data, 2);
        
        lroom_data(:,ro) = lsentence_data;
        rroom_data(:,ro) = rsentence_data;

    end
    kewkurtosis_data{d} = {lroom_data, rroom_data};
    
end

%% obtain hexakikurtosis data
cm = 6;
hexakikurtosis_data = cell(1, length(degree));
for d = 1:length(degree)
    lroom_data = zeros(12, 5);
    rroom_data = zeros(12, 5);
    
    for ro = 1:5
        lsentence_data = zeros(12, 5);
        rsentence_data = zeros(12, 5);

        for s = 1:5 
            lsubband_data = cell(12, 1);
            rsubband_data = cell(12, 1);

            for j = 1:12
                l = fine_statistic{d}{s}{ro}{1}{j}(cm);
                r = fine_statistic{d}{s}{ro}{2}{j}(cm);

                lsubband_data{j} = l;
                rsubband_data{j} = r;

            end

            lsentence_data(:,s) = cell2mat(lsubband_data);
            rsentence_data(:,s) = cell2mat(rsubband_data);
        
        end
        % average values across all used sentences
        lsentence_data = mean(lsentence_data, 2);
        rsentence_data = mean(rsentence_data, 2);
        
        lroom_data(:,ro) = lsentence_data;
        rroom_data(:,ro) = rsentence_data;

    end
    hexakikurtosis_data {d} = {lroom_data, rroom_data};
    
end

data_name = "hos_fine.mat";
save(fullfile(resultsDir, data_name), "mean_data", "variance_data", "skewness_data", "kurtosis_data", "kewkurtosis_data", "hexakikurtosis_data")