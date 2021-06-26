function [power, freq] = calculatePowerPerTrial(EEG_trial, trial_data)
% From the trial data, calculate the power over the frequencies of the signal
% for all electrodes.

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
            % Save the power for the frequencies of  the signal
            pow(:, channel, window_id) = P((f(:,1)>=1 & f(:,1)<=48),1);
        else
            pow(:, channel, window_id) = NaN;
        end
    end
    % Increase indices and window (increase sliding window with
    % 0.5*fs).
    window_id = window_id + 1;
    window = window+0.5*EEG_trial.srate;
end

% Change frequency variable for frequencies of the signal.
freq = f(f(:,1)>=1 & f(:,1)<=48);
% Average power per channel over windows.
power = mean(pow, 3, 'omitnan');

end