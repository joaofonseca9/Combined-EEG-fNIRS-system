function [power, freq] = calculatePowerPerEpoch(EEG_epoch, epoch_data, th)
% From the trial data, calculate the power over the frequencies of the signal
% for all electrodes.

% Using a Hann window.
window = 1:0.8*EEG_epoch.srate;
% Select the data of this specific window [channel x time].
data_window = epoch_data(:, window);

% Channel loop.
for channel = 1:size(data_window, 1)
    % Calculate PSD
    [P, f] = periodogram(data_window(channel, :),...
        hann(size(data_window, 2)),...
        2^(2 + nextpow2(size(data_window, 2))), EEG_epoch.srate);
    
    % Save the power for the frequencies of  the signal.
    pow(:, channel) = P((f(:,1)>=4 & f(:,1)<=48),1);
    
end

% Change frequency variable for frequencies of the signal.
freq = f(f(:,1)>=4 & f(:,1)<=48);

% For the bad channels, give NaN value.
power = pow;
for channel = 1:size(data_window, 1)
    if  pow(channel) > th(channel)
        power(:, channel) = NaN;
    end
end
end