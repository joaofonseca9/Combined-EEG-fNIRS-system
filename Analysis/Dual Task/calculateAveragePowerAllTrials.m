function [power, freq] = calculateAveragePowerAllTrials(EEG, event_samp,...
    startTask, endTask, keypresses, th)
% Loop through the power from the individual trials and average them.
for trial=1:length(startTask)
    
if trial==1
    size_power_allEpochs = 1;
end

% Get the keypresses within that trial.
keypresses_trial = keypresses(keypresses > startTask(trial)...
    & keypresses < endTask(trial));
keypresses_times = event_samp(keypresses_trial);

% Epoch the data into the different keypresses.
for epoch = 1:length(keypresses_times)

    EEG_epoch = pop_select(EEG, 'point',...
        [keypresses_times(epoch)-floor(0.4*EEG.srate)...
        keypresses_times(epoch)+ceil(0.4*EEG.srate)]);
    epoch_data = EEG_epoch.data;

    [power_oneEpoch, freq] = calculatePowerPerEpoch(EEG_epoch,...
        epoch_data, th);

    power_allEpochs(:, :, size_power_allEpochs) = power_oneEpoch;
    size_power_allEpochs = size_power_allEpochs+1;

end

% Take the average of every epoch.
power = mean(power_allEpochs, 3, 'omitnan');

end
end