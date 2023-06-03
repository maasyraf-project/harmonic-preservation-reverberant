clc
clear all
close all

% load all packages
addpath(genpath('additional-packages'));
addpath(genpath('stimuli'));
addpath(genpath('results'));

resultsDir = strcat(pwd,'\results\room_statistical_analysis\');
resultsFilenames = dir(resultsDir);

hosDataFilenames = {};
for i = 1:length(resultsFilenames)
    if contains(resultsFilenames(i).name, "hos_", 'IgnoreCase', true)
        hosDataFilenames{end+1} = resultsFilenames(i).name;
    end
end

hosDataLabel = [];
for i = 1:length(hosDataFilenames)
    label = erase(hosDataFilenames(i), "hos_");
    label = erase(label, ".mat");
    hosDataLabel = [hosDataLabel, string(label)]; 
end

degree = ["0", "30", "45", "60", "90", "min30", "min45", "min60", "min90"];
room = ["Anechoic", "Room A", "Room B", "Room C", "Room D"];

%% plot kurtosis data
for deg = 1%:length(degree)
    lowfig = figure;
    highfig = figure;
    for i = 1:length(hosDataFilenames)
        hosDataName = strcat(pwd,'\results\room_statistical_analysis\', hosDataFilenames{i});
        load(hosDataName)
        data = kurtosis_data;

        for j = 1:12
            l = data{deg}{1}(j,:);
            r = data{deg}{2}(j,:);
            
            if j <= 6
                set(0, "CurrentFigure", lowfig)
                subplot(6,2,j+(j-1))
                plot(l, "-x", "MarkerSize", 5)
                ylabel("Value")
                hold on
                xticklabels(room);
                grid on

                if j+(j-1) == 1
                    axP = get(gca,'Position');
                    legend(hosDataLabel, Orientation="horizontal", Location="northoutside")
                    set(gca, 'Position', axP)
                end

                if j+(j-1) == 11
                    xticklabels(room)
                else
                    xticklabels(["","","","",""])
                end
                

                subplot(6,2,2*j)
                plot(r, "-x", "MarkerSize", 5)
                ylabel("Value")
                hold on
                xticklabels(room)
                grid on

                if 2*j == 2
                    axP = get(gca,'Position');
                    legend(hosDataLabel, Orientation="horizontal", Location="northoutside")
                    set(gca, 'Position', axP)
                end

                if 2*j == 12
                    xticklabels(room)
                else
                    xticklabels(["","","","",""])
                end

            end
            

            if j > 6
                set(0, "CurrentFigure", highfig)
                subplot(6,2,2*j-13)
                plot(l, "-x", "MarkerSize", 5)
                ylabel("Value")
                hold on
                xticklabels(room)
                grid on

                if 2*j-13 == 1
                    axP = get(gca,'Position');
                    legend(hosDataLabel, Orientation="horizontal", Location="northoutside")
                    set(gca, 'Position', axP)
                end

                if 2*j-13 == 11
                    xticklabels(room)
                else
                    xticklabels(["","","","",""])
                end

                subplot(6,2,2*j-12)
                plot(r, "-x", "MarkerSize", 5)
                ylabel("Value")
                hold on
                xticklabels(room)
                grid on

                if 2*j-12 == 2
                    axP = get(gca,'Position');
                    legend(hosDataLabel, Orientation="horizontal", Location="northoutside")
                    set(gca, 'Position', axP)
                end

                if 2*j-12 == 12
                    xticklabels(room)
                else
                    xticklabels(["","","","",""])
                end

            end
            

        end
        clear kurtosis_data

    end
end