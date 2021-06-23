%% Analysis of the EEG signals - topoplots.

clear; clc; close all;
addpath('C:\Users\maria\OneDrive\Documentos\GitHub\Combined-EEG-fNIRS-system\Analysis');

laptop = 'laptopMariana';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);
eeglab;
ft_defaults;
results_path = 'C:\Users\maria\OneDrive\Ambiente de Trabalho\Automaticity Results\Topoplots';

subrec = ["28" "04"];

% Loop through every subject.
for subject = 1:size(subrec, 1)
    sub = subrec(subject, 1);
    rec = subrec(subject, 2);
  
    % Load the subject's EEG signals.
    load([mainpath_in, '\pre-processed\sub-', char(sub), '\eeg\sub-',...
        char(sub), '_rec-', char(rec), '_eeg_divided.mat']);
    
    % Separate into the four different tasks.
    EEG_AutoUncued = EEG_divided.EEG_AutoNoCue;
    EEG_NonAutoUncued = EEG_divided.EEG_NonAutoNoCue;
    EEG_AutoCued = EEG_divided.EEG_AutoCue;
    EEG_NonAutoCued = EEG_divided.EEG_NonAutoCue;
    
    %% Auto Uncued.

    event_samp  = [EEG_AutoUncued.event.latency];
    startBaseline = find(strcmp({EEG_AutoUncued.event.type}, 'boundary')==1);
    startBaseline(12) = [];
    startBaseline(1) = [];
    startTask = find(strcmp({EEG_AutoUncued.event.type}, 's1703')==1);
    endTask = find(strcmp({EEG_AutoUncued.event.type}, 's1711')==1);

    % Get the power spectrum density (PSD) averaged over all trials.
    % For the baseline.
    [power_base_theta, power_base_alpha, power_base_beta,...
        power_base_gamma, freq_base_theta, freq_base_alpha,...
        freq_base_beta, freq_base_gamma] =...
        calculateAveragePowerAllTrials(EEG_AutoUncued, event_samp,...
        startBaseline, startTask);
    % For the task.
    [power_theta, power_alpha, power_beta, power_gamma, freq_theta,...
        freq_alpha, freq_beta, freq_gamma] =...
        calculateAveragePowerAllTrials(EEG_AutoUncued, event_samp,...
        startTask, endTask);
    
    % Calculate the ERD/ERS for each of the frequency bands above.
    ERD_ERS_theta = (power_theta - power_base_theta)./power_base_theta; 
    ERD_ERS_alpha = (power_alpha - power_base_alpha)./power_base_alpha; 
    ERD_ERS_beta = (power_beta - power_base_beta)./power_base_beta; 
    ERD_ERS_gamma = (power_gamma - power_base_gamma)./power_base_gamma; 

    % Topographic distribution of the frequency bands over the head
    % (topoplot).
    figure;
    subplot(2, 2, 1);
    text(-0.13, 0.7, 'Theta', 'FontSize', 18);
    topoplot(ERD_ERS_theta, EEG_AutoUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 2);
    text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
    topoplot(ERD_ERS_alpha, EEG_AutoUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 3);
    text(-0.1, 0.7, 'Beta', 'FontSize', 18)
    topoplot(ERD_ERS_beta, EEG_AutoUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 4);
    text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
    topoplot(ERD_ERS_gamma, EEG_AutoUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    
    % Save the values onto a allSubjects variable.
    autouncued_ERD_ERS_theta_allSubjects(:, subject) = ERD_ERS_theta;
    autouncued_ERD_ERS_alpha_allSubjects(:, subject) = ERD_ERS_alpha;
    autouncued_ERD_ERS_beta_allSubjects(:, subject) = ERD_ERS_beta;
    autouncued_ERD_ERS_gamma_allSubjects(:, subject) = ERD_ERS_gamma;
    
    %% Non-Auto Uncued.

    event_samp  = [EEG_NonAutoUncued.event.latency];
    startBaseline = find(strcmp({EEG_NonAutoUncued.event.type}, 'boundary')==1);
    startBaseline(12) = [];
    startBaseline(1) = [];
    startTask = find(strcmp({EEG_NonAutoUncued.event.type}, 's1705')==1);
    endTask = find(strcmp({EEG_NonAutoUncued.event.type}, 's1713')==1);
    
    % Get the power spectrum density (PSD) averaged over all trials.
    % For the baseline.
    [power_base_theta, power_base_alpha, power_base_beta,...
        power_base_gamma, freq_base_theta, freq_base_alpha,...
        freq_base_beta, freq_base_gamma] =...
        calculateAveragePowerAllTrials(EEG_NonAutoUncued, event_samp,...
        startBaseline, startTask);
    % For the task.
    [power_theta, power_alpha, power_beta, power_gamma, freq_theta,...
        freq_alpha, freq_beta, freq_gamma] =...
        calculateAveragePowerAllTrials(EEG_NonAutoUncued, event_samp,...
        startTask, endTask);
    
    % Calculate the ERD/ERS for each of the frequency bands above.
    ERD_ERS_theta = (power_theta - power_base_theta)./power_base_theta; 
    ERD_ERS_alpha = (power_alpha - power_base_alpha)./power_base_alpha; 
    ERD_ERS_beta = (power_beta - power_base_beta)./power_base_beta; 
    ERD_ERS_gamma = (power_gamma - power_base_gamma)./power_base_gamma; 
    
    % Topographic distribution of the frequency bands over the head
    % (topoplot).
    figure;
    subplot(2, 2, 1);
    text(-0.13, 0.7, 'Theta', 'FontSize', 18);
    topoplot(ERD_ERS_theta, EEG_NonAutoUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 2);
    text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
    topoplot(ERD_ERS_alpha, EEG_NonAutoUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 3);
    text(-0.1, 0.7, 'Beta', 'FontSize', 18)
    topoplot(ERD_ERS_beta, EEG_NonAutoUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 4);
    text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
    topoplot(ERD_ERS_gamma, EEG_NonAutoUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    
    % Save the values onto a allSubjects variable.
    nonautouncued_ERD_ERS_theta_allSubjects(:, subject) = ERD_ERS_theta;
    nonautouncued_ERD_ERS_alpha_allSubjects(:, subject) = ERD_ERS_alpha;
    nonautouncued_ERD_ERS_beta_allSubjects(:, subject) = ERD_ERS_beta;
    nonautouncued_ERD_ERS_gamma_allSubjects(:, subject) = ERD_ERS_gamma;
    
    %% Auto Cued.

    event_samp  = [EEG_AutoCued.event.latency];
    startBaseline = find(strcmp({EEG_AutoCued.event.type}, 'boundary')==1);
    startBaseline(12) = [];
    startBaseline(1) = [];
    startTask = find(strcmp({EEG_AutoCued.event.type}, 's1702')==1);
    endTask = find(strcmp({EEG_AutoCued.event.type}, 's1710')==1);

    % Get the power spectrum density (PSD) averaged over all trials.
    % For the baseline.
    [power_base_theta, power_base_alpha, power_base_beta,...
        power_base_gamma, freq_base_theta, freq_base_alpha,...
        freq_base_beta, freq_base_gamma] =...
        calculateAveragePowerAllTrials(EEG_AutoCued, event_samp,...
        startBaseline, startTask);
    % For the task.
    [power_theta, power_alpha, power_beta, power_gamma, freq_theta,...
        freq_alpha, freq_beta, freq_gamma] =...
        calculateAveragePowerAllTrials(EEG_AutoCued, event_samp,...
        startTask, endTask);
    
    % Calculate the ERD/ERS for each of the frequency bands above.
    ERD_ERS_theta = (power_theta - power_base_theta)./power_base_theta; 
    ERD_ERS_alpha = (power_alpha - power_base_alpha)./power_base_alpha; 
    ERD_ERS_beta = (power_beta - power_base_beta)./power_base_beta; 
    ERD_ERS_gamma = (power_gamma - power_base_gamma)./power_base_gamma; 
    
    % Topographic distribution of the frequency bands over the head
    % (topoplot).
    figure;
    subplot(2, 2, 1);
    text(-0.13, 0.7, 'Theta', 'FontSize', 18);
    topoplot(ERD_ERS_theta, EEG_AutoCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 2);
    text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
    topoplot(ERD_ERS_alpha, EEG_AutoCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 3);
    text(-0.1, 0.7, 'Beta', 'FontSize', 18)
    topoplot(ERD_ERS_beta, EEG_AutoCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;   
    subplot(2, 2, 4);
    text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
    topoplot(ERD_ERS_gamma, EEG_AutoCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    
    % Save the values onto a allSubjects variable.
    autocued_ERD_ERS_theta_allSubjects(:, subject) = ERD_ERS_theta;
    autocued_ERD_ERS_alpha_allSubjects(:, subject) = ERD_ERS_alpha;
    autocued_ERD_ERS_beta_allSubjects(:, subject) = ERD_ERS_beta;
    autocued_ERD_ERS_gamma_allSubjects(:, subject) = ERD_ERS_gamma;
    
    %% Non-Auto Cued.

    event_samp  = [EEG_NonAutoCued.event.latency];
    startBaseline = find(strcmp({EEG_NonAutoCued.event.type}, 'boundary')==1);
    startBaseline(12) = [];
    startBaseline(1) = [];
    startTask = find(strcmp({EEG_NonAutoCued.event.type}, 's1704')==1);
    endTask = find(strcmp({EEG_NonAutoCued.event.type}, 's1712')==1);

    % Get the power spectrum density (PSD) averaged over all trials.  
    % For the baseline.
    [power_base_theta, power_base_alpha, power_base_beta,...
        power_base_gamma, freq_base_theta, freq_base_alpha,...
        freq_base_beta, freq_base_gamma] =...
        calculateAveragePowerAllTrials(EEG_NonAutoCued, event_samp,...
        startBaseline, startTask);
    % For the task.
    [power_theta, power_alpha, power_beta, power_gamma, freq_theta,...
        freq_alpha, freq_beta, freq_gamma] =...
        calculateAveragePowerAllTrials(EEG_NonAutoCued, event_samp,...
        startTask, endTask);
    
    % Calculate the ERD/ERS for each of the frequency bands above.
    ERD_ERS_theta = (power_theta - power_base_theta)./power_base_theta; 
    ERD_ERS_alpha = (power_alpha - power_base_alpha)./power_base_alpha; 
    ERD_ERS_beta = (power_beta - power_base_beta)./power_base_beta; 
    ERD_ERS_gamma = (power_gamma - power_base_gamma)./power_base_gamma; 

    % Topographic distribution of the frequency bands over the head
    % (topoplot).
    figure;
    subplot(2, 2, 1);
    text(-0.13, 0.7, 'Theta', 'FontSize', 18);
    topoplot(ERD_ERS_theta, EEG_NonAutoCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 2);
    text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
    topoplot(ERD_ERS_alpha, EEG_NonAutoCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 3);
    text(-0.1, 0.7, 'Beta', 'FontSize', 18)
    topoplot(ERD_ERS_beta, EEG_NonAutoCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar; 
    subplot(2, 2, 4);
    text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
    topoplot(ERD_ERS_gamma, EEG_NonAutoCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    
    % Save the values onto a allSubjects variable.
    nonautocued_ERD_ERS_theta_allSubjects(:, subject) = ERD_ERS_theta;
    nonautocued_ERD_ERS_alpha_allSubjects(:, subject) = ERD_ERS_alpha;
    nonautocued_ERD_ERS_beta_allSubjects(:, subject) = ERD_ERS_beta;
    nonautocued_ERD_ERS_gamma_allSubjects(:, subject) = ERD_ERS_beta;
    
    disp(['These are the topoplots for subject ', char(sub), '.']);
    disp('Press any key to move onto the next subject.');
    pause;
    close all;
    
end

% Get the power spectrum density (PSD) averaged over all subjects.
% Auto Uncued.
autouncued_ERD_ERS_theta = mean(autouncued_ERD_ERS_theta_allSubjects, 2);
autouncued_ERD_ERS_alpha = mean(autouncued_ERD_ERS_alpha_allSubjects, 2);
autouncued_ERD_ERS_beta = mean(autouncued_ERD_ERS_beta_allSubjects, 2);
autouncued_ERD_ERS_gamma = mean(autouncued_ERD_ERS_gamma_allSubjects, 2);
% Non-Auto Uncued.
nonautouncued_ERD_ERS_theta = mean(nonautouncued_ERD_ERS_theta_allSubjects, 2);
nonautouncued_ERD_ERS_alpha = mean(nonautouncued_ERD_ERS_alpha_allSubjects, 2);
nonautouncued_ERD_ERS_beta = mean(nonautouncued_ERD_ERS_beta_allSubjects, 2);
nonautouncued_ERD_ERS_gamma = mean(nonautouncued_ERD_ERS_gamma_allSubjects, 2);
% Auto Cued.
autocued_ERD_ERS_theta = mean(autocued_ERD_ERS_theta_allSubjects, 2);
autocued_ERD_ERS_alpha = mean(autocued_ERD_ERS_alpha_allSubjects, 2);
autocued_ERD_ERS_beta = mean(autocued_ERD_ERS_beta_allSubjects, 2);
autocued_ERD_ERS_gamma = mean(autocued_ERD_ERS_gamma_allSubjects, 2);
% Non-Auto Cued.
nonautocued_ERD_ERS_theta = mean(nonautocued_ERD_ERS_theta_allSubjects, 2);
nonautocued_ERD_ERS_alpha = mean(nonautocued_ERD_ERS_alpha_allSubjects, 2);
nonautocued_ERD_ERS_beta = mean(nonautocued_ERD_ERS_beta_allSubjects, 2);
nonautocued_ERD_ERS_gamma = mean(nonautocued_ERD_ERS_gamma_allSubjects, 2);

% Topographic distribution of the frequency bands over the head
% (topoplot).

% Auto Uncued.
figure;
subplot(2, 2, 1);
text(-0.13, 0.7, 'Theta', 'FontSize', 18);
topoplot(autouncued_ERD_ERS_theta, EEG_AutoUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 2);
text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
topoplot(autouncued_ERD_ERS_alpha, EEG_AutoUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 3);
text(-0.1, 0.7, 'Beta', 'FontSize', 18)
topoplot(autouncued_ERD_ERS_beta, EEG_AutoUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar; 
subplot(2, 2, 4);
text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
topoplot(autouncued_ERD_ERS_gamma, EEG_AutoUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar;

% Save figure.
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, fullfile(results_path, 'autouncued_erders'),'png');

% Non-Auto Uncued.
figure;
subplot(2, 2, 1);
text(-0.13, 0.7, 'Theta', 'FontSize', 18);
topoplot(nonautouncued_ERD_ERS_theta, EEG_NonAutoUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 2);
text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
topoplot(nonautouncued_ERD_ERS_alpha, EEG_NonAutoUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 3);
text(-0.1, 0.7, 'Beta', 'FontSize', 18)
topoplot(nonautouncued_ERD_ERS_beta, EEG_NonAutoUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar; 
subplot(2, 2, 4);
text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
topoplot(nonautouncued_ERD_ERS_gamma, EEG_NonAutoUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar;

% Save figure.
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, fullfile(results_path, 'nonautouncued_erders'),'png');

% Auto Cued.
figure;
subplot(2, 2, 1);
text(-0.13, 0.7, 'Theta', 'FontSize', 18);
topoplot(autocued_ERD_ERS_theta, EEG_AutoCued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 2);
text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
topoplot(autocued_ERD_ERS_alpha, EEG_AutoCued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 3);
text(-0.1, 0.7, 'Beta', 'FontSize', 18)
topoplot(autocued_ERD_ERS_beta, EEG_AutoCued.chanlocs, 'electrodes', 'ptslabels');
colorbar; 
subplot(2, 2, 4);
text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
topoplot(autocued_ERD_ERS_gamma, EEG_AutoCued.chanlocs, 'electrodes', 'ptslabels');
colorbar;

% Save figure.
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, fullfile(results_path, 'autocued_erders'),'png');

% Non-Auto Cued.
figure;
subplot(2, 2, 1);
text(-0.13, 0.7, 'Theta', 'FontSize', 18);
topoplot(nonautocued_ERD_ERS_theta, EEG_NonAutoCued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 2);
text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
topoplot(nonautocued_ERD_ERS_alpha, EEG_NonAutoCued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 3);
text(-0.1, 0.7, 'Beta', 'FontSize', 18)
topoplot(nonautocued_ERD_ERS_beta, EEG_NonAutoCued.chanlocs, 'electrodes', 'ptslabels');
colorbar; 
subplot(2, 2, 4);
text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
topoplot(nonautocued_ERD_ERS_gamma, EEG_NonAutoCued.chanlocs, 'electrodes', 'ptslabels');
colorbar;

% Save figure.
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, fullfile(results_path, 'nonautocued_erders'),'png');

disp('This was the end of individual subjects.');
disp('These are the topoplots for the average of all subjects.');

%% Functions

% Loop through the power from the individual trials and average them.
function [power_theta, power_alpha, power_beta, power_gamma, freq_theta,...
    freq_alpha, freq_beta, freq_gamma] =...
    calculateAveragePowerAllTrials(EEG, event_samp, startTask, endTask)

    for trial=1:length(startTask)
    
        title = char(strcat('Trial_', string(trial)));
        startTask_times = event_samp(startTask(trial));
        endTask_times = event_samp(endTask(trial));

        EEG_trial = pop_select(EEG, 'point', [startTask_times endTask_times]);
        trial_data = EEG_trial.data;

        [power_theta_oneTrial, power_alpha_oneTrial, power_beta_oneTrial,...
            power_gamma_oneTrial, freq_theta, freq_alpha, freq_beta,...
            freq_gamma] = calculatePowerPerTrial(EEG_trial, trial_data);
            
        power_theta_allTrials(:, trial) = power_theta_oneTrial;
        power_alpha_allTrials(:, trial) = power_alpha_oneTrial;
        power_beta_allTrials(:, trial) = power_beta_oneTrial;
        power_gamma_allTrials(:, trial) = power_gamma_oneTrial;
    
    end
    
    power_theta = mean(power_theta_allTrials, 2);
    power_alpha = mean(power_alpha_allTrials, 2);
    power_beta = mean(power_beta_allTrials, 2);
    power_gamma = mean(power_gamma_allTrials, 2);

end

% From the trial data, calculate the power over each frequency band for all
% electrodes.
% Theta - 4 to 8 Hz; Alpha - 8 to 13 Hz; Beta - 13 to 32 Hz; Gamma - 32 to
% 48 Hz.
function [power_theta, power_alpha, power_beta, power_gamma, freq_theta,...
    freq_alpha, freq_beta, freq_gamma] =...
    calculatePowerPerTrial(EEG_trial, trial_data)
    
    % Using a sliding Hann window.
    window_id = 1;
    window = 1:1*EEG_trial.srate;
    while window(end) <= size(trial_data, 2)
        % Select the data of this specific window [channel x time].
        data_window = trial_data(:, window);
        
        % Channel loop.
        for channel = 1:size(data_window, 1) 
            % If window is NOT removed because of badchannel (=NaN)
            if isempty(find(isnan(data_window(channel, :))))
                % Calculate PSD
                [P, f] = periodogram(data_window(channel, :),...
                    hann(size(data_window, 2)),...
                    2^(2 + nextpow2(size(data_window, 2))), EEG_trial.srate);
                % Save the power for the frequencies of interest in the 
                % different pow variables (all windows will be saved
                % here)
                pow_theta(:, channel, window_id) = P((f(:,1)>=4 & f(:,1)<=8),1);
                pow_alpha(:, channel, window_id) = P((f(:,1)>=8 & f(:,1)<=13),1);
                pow_beta(:, channel, window_id) = P((f(:,1)>=13 & f(:,1)<=32),1);
                pow_gamma(:, channel, window_id) = P((f(:,1)>=32 & f(:,1)<=48),1);
            else
                pow_theta(:, channel, window_id) = NaN;
                pow_alpha(:, channel, window_id) = NaN;
                pow_beta(:, channel, window_id) = NaN;
                pow_gamma(:, channel, window_id) = NaN;
            end
        end
        % Increase indices and window (increase sliding window with
        % 0.5*fs).
        window_id = window_id + 1;
        window = window+0.5*EEG_trial.srate;
    end
    
    % Change frequency variable for frequencies of interest.
    freq_theta = f(f(:,1)>=4 & f(:,1)<=8);
    freq_alpha = f(f(:,1)>=8 & f(:,1)<=13);
    freq_beta = f(f(:,1)>=13 & f(:,1)<=32);
    freq_gamma = f(f(:,1)>=32 & f(:,1)<=48);
    % Average power per channel over windows and then average over the
    % different channels.
    power_theta = mean(mean(pow_theta,3,'omitnan'));
    power_alpha = mean(mean(pow_alpha,3,'omitnan'));
    power_beta = mean(mean(pow_beta,3,'omitnan'));
    power_gamma = mean(mean(pow_gamma,3,'omitnan'));

end