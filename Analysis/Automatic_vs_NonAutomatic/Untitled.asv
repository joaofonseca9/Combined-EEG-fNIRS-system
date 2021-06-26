%% Analysis of the EEG signals - topoplots.

clear; clc; close all;
addpath('C:\Users\maria\OneDrive\Documentos\GitHub\Combined-EEG-fNIRS-system\Analysis');

laptop = 'laptopMariana';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);
eeglab;
ft_defaults;
results_path = 'C:\Users\maria\OneDrive\Ambiente de Trabalho\Automaticity Results\Topoplots';

subrec = ["02" "02"];

sub = subrec(1);
rec = subrec(2);

% Load the subject's EEG signals.
load([mainpath_in, '\pre-processed\sub-', char(sub), '\eeg\sub-',...
    char(sub), '_rec-', char(rec), '_eeg_divided.mat']);

% Separate into the four different tasks.
EEG_AutoUncued = EEG_divided.EEG_AutoNoCue;
EEG_NonAutoUncued = EEG_divided.EEG_NonAutoNoCue;
EEG_AutoCued = EEG_divided.EEG_AutoCue;
EEG_NonAutoCued = EEG_divided.EEG_NonAutoCue;

%% Auto Uncued
event_samp  = [EEG_AutoUncued.event.latency];
startBaseline = find(strcmp({EEG_AutoUncued.event.type}, 'boundary')==1);
%    startBaseline(12) = [];
%    startBaseline(1) = [];
startTask = find(strcmp({EEG_AutoUncued.event.type}, 's1703')==1);
endTask = find(strcmp({EEG_AutoUncued.event.type}, 's1711')==1);

for i=length(startBaseline):-1:1
    if ~any(startBaseline(i)+1==startTask(:))
        startBaseline(i) = [];
    end
end

%%
% Get the power spectrum density (PSD) averaged over all trials.
% For the baseline.
% [power_base_theta, power_base_alpha, power_base_beta,...
%     power_base_gamma, freq_base_theta, freq_base_alpha,...
%     freq_base_beta, freq_base_gamma] =...
%     calculateAveragePowerAllTrials(EEG_AutoUncued, event_samp,...
%     startBaseline, startTask, f);
% For the task.
% [power_theta, power_alpha, power_beta, power_gamma, freq_theta,...
%     freq_alpha, freq_beta, freq_gamma] =...
%     calculateAveragePowerAllTrials(EEG_AutoUncued, event_samp,...
%     startTask, endTask);

for trial=1:length(startBaseline)
    
    startTask_times = event_samp(startBaseline(trial));
    endTask_times = event_samp(startTask(trial));
    
    EEG_trial = pop_select(EEG_AutoUncued, 'point', [startTask_times endTask_times]);
    trial_data = EEG_trial.data;
    
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
                [P, ~] = periodogram(data_window(channel, :),...
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

%%

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