%% Analysis of the EEG signals - topoplots.

clear; clc; close all;
addpath('C:\Users\maria\OneDrive\Documentos\GitHub\Combined-EEG-fNIRS-system\Analysis');

laptop = 'laptopMariana';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);
eeglab;
ft_defaults;

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
    startTask = find(strcmp({EEG_AutoUncued.event.type}, 's1703')==1);
    endTask = find(strcmp({EEG_AutoUncued.event.type}, 's1711')==1);

    % Get the power spectrum density (PSD) averaged over all trials.
    [power_theta, power_alpha, power_beta, power_gamma, freq_theta,...
        freq_alpha, freq_beta, freq_gamma] =...
        calculateAveragePowerAllTrials(EEG_AutoUncued, event_samp,...
        startTask, endTask);

    % Topographic distribution of the frequency bands over the head
    % (topoplot).
    figure;
    subplot(2, 2, 1);
    text(-0.13, 0.7, 'Theta', 'FontSize', 18);
    topoplot(power_theta, EEG_AutoUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 2);
    text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
    topoplot(power_alpha, EEG_AutoUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 3);
    text(-0.1, 0.7, 'Beta', 'FontSize', 18)
    topoplot(power_beta, EEG_AutoUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 4);
    text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
    topoplot(power_gamma, EEG_AutoUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    
    % Save the values onto a allSubjects variable.
    autouncued_power_theta_allSubjects(:, subject) = power_theta;
    autouncued_power_alpha_allSubjects(:, subject) = power_alpha;
    autouncued_power_beta_allSubjects(:, subject) = power_beta;
    autouncued_power_gamma_allSubjects(:, subject) = power_gamma;
    
    %% Non-Auto Uncued.

    event_samp  = [EEG_NonAutoUncued.event.latency];
    startTask = find(strcmp({EEG_NonAutoUncued.event.type}, 's1705')==1);
    endTask = find(strcmp({EEG_NonAutoUncued.event.type}, 's1713')==1);
    
    % Get the power spectrum density (PSD) averaged over all trials.
    [power_theta, power_alpha, power_beta, power_gamma, freq_theta,...
        freq_alpha, freq_beta, freq_gamma] =...
        calculateAveragePowerAllTrials(EEG_NonAutoUncued, event_samp,...
        startTask, endTask);
    
    % Topographic distribution of the frequency bands over the head
    % (topoplot).
    figure;
    subplot(2, 2, 1);
    text(-0.13, 0.7, 'Theta', 'FontSize', 18);
    topoplot(power_theta, EEG_NonAutoUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 2);
    text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
    topoplot(power_alpha, EEG_NonAutoUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 3);
    text(-0.1, 0.7, 'Beta', 'FontSize', 18)
    topoplot(power_beta, EEG_NonAutoUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 4);
    text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
    topoplot(power_gamma, EEG_NonAutoUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    
    % Save the values onto a allSubjects variable.
    nonautouncued_power_theta_allSubjects(:, subject) = power_theta;
    nonautouncued_power_alpha_allSubjects(:, subject) = power_alpha;
    nonautouncued_power_beta_allSubjects(:, subject) = power_beta;
    nonautouncued_power_gamma_allSubjects(:, subject) = power_gamma;
    
    %% Auto Cued.

    event_samp  = [EEG_AutoCued.event.latency];
    startTask = find(strcmp({EEG_AutoCued.event.type}, 's1702')==1);
    endTask = find(strcmp({EEG_AutoCued.event.type}, 's1710')==1);

    % Get the power spectrum density (PSD) averaged over all trials.
    [power_theta, power_alpha, power_beta, power_gamma, freq_theta,...
        freq_alpha, freq_beta, freq_gamma] =...
        calculateAveragePowerAllTrials(EEG_AutoCued, event_samp,...
        startTask, endTask);

    % Topographic distribution of the frequency bands over the head
    % (topoplot).
    figure;
    subplot(2, 2, 1);
    text(-0.13, 0.7, 'Theta', 'FontSize', 18);
    topoplot(power_theta, EEG_AutoCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 2);
    text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
    topoplot(power_alpha, EEG_AutoCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 3);
    text(-0.1, 0.7, 'Beta', 'FontSize', 18)
    topoplot(power_beta, EEG_AutoCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;   
    subplot(2, 2, 4);
    text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
    topoplot(power_gamma, EEG_AutoCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    
    % Save the values onto a allSubjects variable.
    autocued_power_theta_allSubjects(:, subject) = power_theta;
    autocued_power_alpha_allSubjects(:, subject) = power_alpha;
    autocued_power_beta_allSubjects(:, subject) = power_beta;
    autocued_power_gamma_allSubjects(:, subject) = power_gamma;
    
    %% Non-Auto Cued.

    event_samp  = [EEG_NonAutoCued.event.latency];
    startTask = find(strcmp({EEG_NonAutoCued.event.type}, 's1704')==1);
    endTask = find(strcmp({EEG_NonAutoCued.event.type}, 's1712')==1);

    % Get the power spectrum density (PSD) averaged over all trials.   
    [power_theta, power_alpha, power_beta, power_gamma, freq_theta,...
        freq_alpha, freq_beta, freq_gamma] =...
        calculateAveragePowerAllTrials(EEG_NonAutoCued, event_samp,...
        startTask, endTask);

    % Topographic distribution of the frequency bands over the head
    % (topoplot).
    figure;
    subplot(2, 2, 1);
    text(-0.13, 0.7, 'Theta', 'FontSize', 18);
    topoplot(power_theta, EEG_NonAutoCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 2);
    text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
    topoplot(power_alpha, EEG_NonAutoCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 3);
    text(-0.1, 0.7, 'Beta', 'FontSize', 18)
    topoplot(power_beta, EEG_NonAutoCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar; 
    subplot(2, 2, 4);
    text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
    topoplot(power_gamma, EEG_NonAutoCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    
    % Save the values onto a allSubjects variable.
    nonautocued_power_theta_allSubjects(:, subject) = power_theta;
    nonautocued_power_alpha_allSubjects(:, subject) = power_alpha;
    nonautocued_power_beta_allSubjects(:, subject) = power_beta;
    nonautocued_power_gamma_allSubjects(:, subject) = power_beta;
    
    disp(['These are the topoplots for subject ', char(sub), '.']);
    disp('Press any key to move onto the next subject.');
    pause;
    close all;
    
end

% Get the power spectrum density (PSD) averaged over all subjects.
% Auto Uncued.
autouncued_power_theta = mean(autouncued_power_theta_allSubjects, 2);
autouncued_power_alpha = mean(autouncued_power_alpha_allSubjects, 2);
autouncued_power_beta = mean(autouncued_power_beta_allSubjects, 2);
autouncued_power_gamma = mean(autouncued_power_gamma_allSubjects, 2);
% Non-Auto Uncued.
nonautouncued_power_theta = mean(nonautouncued_power_theta_allSubjects, 2);
nonautouncued_power_alpha = mean(nonautouncued_power_alpha_allSubjects, 2);
nonautouncued_power_beta = mean(nonautouncued_power_beta_allSubjects, 2);
nonautouncued_power_gamma = mean(nonautouncued_power_gamma_allSubjects, 2);
% Auto Cued.
autocued_power_theta = mean(autocued_power_theta_allSubjects, 2);
autocued_power_alpha = mean(autocued_power_alpha_allSubjects, 2);
autocued_power_beta = mean(autocued_power_beta_allSubjects, 2);
autocued_power_gamma = mean(autocued_power_gamma_allSubjects, 2);
% Non-Auto Cued.
nonautocued_power_theta = mean(nonautocued_power_theta_allSubjects, 2);
nonautocued_power_alpha = mean(nonautocued_power_alpha_allSubjects, 2);
nonautocued_power_beta = mean(nonautocued_power_beta_allSubjects, 2);
nonautocued_power_gamma = mean(nonautocued_power_gamma_allSubjects, 2);

% Topographic distribution of the frequency bands over the head
% (topoplot).
% Auto Uncued.
figure;
subplot(2, 2, 1);
text(-0.13, 0.7, 'Theta', 'FontSize', 18);
topoplot(autouncued_power_theta, EEG_AutoUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 2);
text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
topoplot(autouncued_power_alpha, EEG_AutoUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 3);
text(-0.1, 0.7, 'Beta', 'FontSize', 18)
topoplot(autouncued_power_beta, EEG_AutoUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar; 
subplot(2, 2, 4);
text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
topoplot(autouncued_power_gamma, EEG_AutoUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
% Non-Auto Uncued.
figure;
subplot(2, 2, 1);
text(-0.13, 0.7, 'Theta', 'FontSize', 18);
topoplot(nonautouncued_power_theta, EEG_NonAutoUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 2);
text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
topoplot(nonautouncued_power_alpha, EEG_NonAutoUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 3);
text(-0.1, 0.7, 'Beta', 'FontSize', 18)
topoplot(nonautouncued_power_beta, EEG_NonAutoUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar; 
subplot(2, 2, 4);
text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
topoplot(nonautouncued_power_gamma, EEG_NonAutoUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
% Auto Cued.
figure;
subplot(2, 2, 1);
text(-0.13, 0.7, 'Theta', 'FontSize', 18);
topoplot(autocued_power_theta, EEG_AutoCued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 2);
text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
topoplot(autocued_power_alpha, EEG_AutoCued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 3);
text(-0.1, 0.7, 'Beta', 'FontSize', 18)
topoplot(autocued_power_beta, EEG_AutoCued.chanlocs, 'electrodes', 'ptslabels');
colorbar; 
subplot(2, 2, 4);
text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
topoplot(autocued_power_gamma, EEG_AutoCued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
% Non-Auto Cued.
figure;
subplot(2, 2, 1);
text(-0.13, 0.7, 'Theta', 'FontSize', 18);
topoplot(nonautocued_power_theta, EEG_NonAutoCued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 2);
text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
topoplot(nonautocued_power_alpha, EEG_NonAutoCued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 3);
text(-0.1, 0.7, 'Beta', 'FontSize', 18)
topoplot(nonautocued_power_beta, EEG_NonAutoCued.chanlocs, 'electrodes', 'ptslabels');
colorbar; 
subplot(2, 2, 4);
text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
topoplot(nonautocued_power_gamma, EEG_NonAutoCued.chanlocs, 'electrodes', 'ptslabels');
colorbar;

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