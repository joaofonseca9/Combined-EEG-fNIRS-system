function [power, freq] = calculateAveragePowerAllTrials(EEG, event_samp,...
    startTask, endTask)
% Loop through the power from the individual trials and average them.

for trial=1:length(startTask)

    title = char(strcat('Trial_', string(trial)));
    startTask_times = event_samp(startTask(trial));
    endTask_times = event_samp(endTask(trial));

    EEG_trial = pop_select(EEG, 'point', [startTask_times endTask_times]);
    trial_data = EEG_trial.data;

    [power_oneTrial, freq] = calculatePowerPerTrial(EEG_trial, trial_data);

    power_allTrials(:, :, trial) = power_oneTrial;

end

power = mean(power_allTrials, 3);

end