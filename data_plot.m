clc
clear all
close all

% load all packages
addpath(genpath('additional-packages'));
addpath(genpath('stimuli'));
addpath(genpath('results'));

resultsDir = strcat(pwd,'\results\room_statistical_analysis\');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

resultsFilenames = dir(resultsDir);
resultsFilenames = resultsFilenames(3:end);

degree = ["0", "30", "45", "60", "90", "min30", "min45", "min60", "min90"];
room = ["Anechoic", "Room A", "Room B", "Room C", "Room D"];

for i = 1:length(resultsFilenames)
    load(resultsFilenames(1).name)
    

end
