function [power_theta, power_alpha, power_beta, power_gamma, freq_theta,...
    freq_alpha, freq_beta, freq_gamma] =...
    calculateAveragePowerERDERSAllTrials(EEG, event_samp, startTask, endTask)
% Loop through the power from the individual trials and average them.
% From the trial data, calculate the power over each frequency band for all
% electrodes.
% Theta - 4 to 8 Hz; Alpha - 8 to 13 Hz; Beta - 13 to 32 Hz; Gamma - 32 to
% 48 Hz.

for trial=1:length(startTask)
    
startTask_times = event_samp(startTask(trial));
endTask_times = event_samp(endTask(trial));

EEG_trial = pop_select(EEG, 'point', [startTask_times endTask_times]);
trial_data = EEG_trial.data;

% Using a sliding Hann window
window_id = 1;
window = 1:1*EEG_trial.srate;
while window(end) <= size(trial_data, 2)
    % Select the data of this specific window [channel x time]
    data_window = trial_data(:, window);

    % Channel loop
    for channel = 1:size(data_window, 1)
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
    end
    % Increase indices and window (increase sliding window with
    % 0.5*fs)
    window_id = window_id + 1;
    window = window+0.5*EEG_trial.srate;
end

% Change frequency variable for frequencies of interest
freq_theta = f(f(:,1)>=4 & f(:,1)<=8);
freq_alpha = f(f(:,1)>=8 & f(:,1)<=13);
freq_beta = f(f(:,1)>=13 & f(:,1)<=32);
freq_gamma = f(f(:,1)>=32 & f(:,1)<=48);
% Average power per channel over windows and then average over the
% different channels
power_theta_oneTrial = mean(mean(pow_theta,3));
power_alpha_oneTrial = mean(mean(pow_alpha,3));
power_beta_oneTrial = mean(mean(pow_beta,3));
power_gamma_oneTrial = mean(mean(pow_gamma,3));

power_theta_allTrials(:, trial) = power_theta_oneTrial;
power_alpha_allTrials(:, trial) = power_alpha_oneTrial;
power_beta_allTrials(:, trial) = power_beta_oneTrial;
power_gamma_allTrials(:, trial) = power_gamma_oneTrial;

end

% Take the average of every trials
power_theta = mean(power_theta_allTrials, 2);
power_alpha = mean(power_alpha_allTrials, 2);
power_beta = mean(power_beta_allTrials, 2);
power_gamma = mean(power_gamma_allTrials, 2);

end