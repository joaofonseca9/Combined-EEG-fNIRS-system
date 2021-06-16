clear; clc; close all;
addpath('C:\Users\maria\OneDrive\Documentos\GitHub\Combined-EEG-fNIRS-system\Analysis');

laptop = 'laptopMariana';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);
eeglab;
ft_defaults;

sub='04';
rec='01';

file = getFileNames(mainpath_out, sub, rec);

load(file.EEG_divided, 'EEG_divided');

EEG_AutoUncued = EEG_divided.EEG_AutoNoCue;
EEG_NonAutoUncued = EEG_divided.EEG_NonAutoNoCue;
EEG_AutoCued = EEG_divided.EEG_AutoCue;
EEG_NonAutoCued = EEG_divided.EEG_NonAutoCue;

%% Auto Uncued

[power_theta, power_alpha, power_beta, freq_theta, freq_alpha,...
    freq_beta] = calculateAveragePowerAllTrials(EEG_AutoUncued);

figure;
subplot(1, 3, 1);
text(-0.13, 0.7, 'Theta', 'FontSize', 18);
topoplot(power_theta, EEG_AutoUncued.chanlocs, 'electrodes', 'on');
colorbar;
subplot(1, 3, 2);
text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
topoplot(power_alpha, EEG_AutoUncued.chanlocs, 'electrodes', 'on');
colorbar;
subplot(1, 3, 3);
text(-0.1, 0.7, 'Beta', 'FontSize', 18)
topoplot(power_beta, EEG_AutoUncued.chanlocs, 'electrodes', 'on');
colorbar;   

%% Functions

% Loop through the power from the individual trials and average them.
function [power_theta, power_alpha, power_beta, freq_theta, freq_alpha,...
    freq_beta] = calculateAveragePowerAllTrials(EEG)

    event_samp  = [EEG.event.latency];
    startTask = find(strcmp({EEG.event.type}, 's1703')==1);
    endTask = find(strcmp({EEG.event.type}, 's1711')==1);
    
    for trial=1:length(startTask)
    
        title = char(strcat('Trial_', string(trial)));
        startTask_times = event_samp(startTask(trial));
        endTask_times = event_samp(endTask(trial));

        EEG_trial = pop_select(EEG, 'point', [startTask_times endTask_times]);
        trial_data = EEG_trial.data;

        [power_theta_oneTrial, power_alpha_oneTrial, power_beta_oneTrial,...
            freq_theta, freq_alpha, freq_beta] =...
            calculatePowerPerTrial(EEG_trial, trial_data);
            
        power_theta_allTrials(:, trial) = power_theta_oneTrial;
        power_alpha_allTrials(:, trial) = power_alpha_oneTrial;
        power_beta_allTrials(:, trial) = power_beta_oneTrial;
    
    end
    
    power_theta = mean(power_theta_allTrials, 2);
    power_alpha = mean(power_alpha_allTrials, 2);
    power_beta = mean(power_beta_allTrials, 2);

end

% From the trial data, calculate the power over each frequency band for all
% electrodes.
% Theta - 4 to 8 Hz; Alpha - 8 to 13 Hz; Beta - 13 to 32 Hz.
function [power_theta, power_alpha, power_beta, freq_theta, freq_alpha,...
    freq_beta] = calculatePowerPerTrial(EEG_trial, trial_data)
    
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
            else
                pow_theta(:, channel, window_id) = NaN;
                pow_alpha(:, channel, window_id) = NaN;
                pow_beta(:, channel, window_id) = NaN;
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
    % Average power per channel over windows.
    power_theta = mean(mean(pow_theta,3,'omitnan'));
    power_alpha = mean(mean(pow_alpha,3,'omitnan'));
    power_beta = mean(mean(pow_beta,3,'omitnan'));

end