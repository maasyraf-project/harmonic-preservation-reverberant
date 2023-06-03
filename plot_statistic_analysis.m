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


for d = 1:length(degree)

    FINE = cell(1, numSentences);
    FINE_LP = cell(1, numSentences);
    ENV = cell(1, numSentences);
    ENV_LP = cell(1, numSentences);

    for i = d:9:length(resultsFilenames)
        load(resultsFilenames(i).name)
        disp(strcat("Processing ", resultsFilenames(i).name));

        fine = cell(1, 5);
        fine_lp = cell(1, 5);
        env = cell(1, 5);
        env_lp = cell(1, 5);

        for j = 1:5

            l_fine = cell(12, 1);
            r_fine = cell(12, 1);

            l_fine_lp = cell(12, 1);
            r_fine_lp = cell(12, 1);

            l_env = cell(12, 1);
            r_env = cell(12, 1);

            l_env_lp = cell(12, 1);
            r_env_lp = cell(12, 1);

            for k = 1:12
                % obtain mean across channel
                lcm_fine = mean(fine_cm_data{j}{1}{k}, 2);
                rcm_fine = mean(fine_cm_data{j}{2}{k}, 2);

                lcm_fine_lp = mean(fine_lp_cm_data{j}{1}{k}, 2);
                rcm_fine_lp = mean(fine_lp_cm_data{j}{2}{k}, 2);

                lcm_env = mean(env_lp_cm_data{j}{1}{k}, 2);
                rcm_env = mean(env_lp_cm_data{j}{2}{k}, 2);

                lcm_env_lp = mean(env_lp_lp_cm_data{j}{1}{k}, 2);
                rcm_env_lp = mean(env_lp_lp_cm_data{j}{2}{k}, 2);

                % add on the provided array
                l_fine{k} = lcm_fine;
                r_fine{k} = rcm_fine;

                l_fine_lp{k} = lcm_fine_lp;
                r_fine_lp{k} = rcm_fine_lp;

                l_env{k} = lcm_env;
                r_env{k} = rcm_env;

                l_env_lp{k} = lcm_env_lp;
                r_env_lp{k} = rcm_env_lp;

            end

            fine{j} = {l_fine, r_fine};
            fine_lp{j} = {l_fine_lp, r_fine_lp};
            env{j} = {l_env, r_env};
            env_lp{j} = {l_env_lp, r_env_lp};

        end

        FINE{i} = fine;
        FINE_LP{i} = fine_lp;
        ENV{i} = env;
        ENV_LP{i} = env_lp;

    end

    % save the variable
    dataName = strcat(string(numSentences), "_", degree(d), ".mat");
    save(fullfile(saveDir, dataName), "FINE", "FINE_LP", "ENV", "ENV_LP");
    disp(strcat("Saving ... ", dataName));
    disp(" ")

end