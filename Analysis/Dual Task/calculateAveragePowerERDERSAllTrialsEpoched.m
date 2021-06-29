function [power_theta, power_alpha, power_beta, power_gamma, freq_theta,...
    freq_alpha, freq_beta, freq_gamma] =...
    calculateAveragePowerERDERSAllTrialsEpoched(EEG, event_samp, startTask, endTask, keypresses, th)
% Loop through the power from the individual trials and average them.
% From the trial data, calculate the power over each frequency band for all
% electrodes.
% Theta - 4 to 8 Hz; Alpha - 8 to 13 Hz; Beta - 13 to 32 Hz; Gamma - 32 to
% 48 Hz.

for trial=1:length(startTask)
    
if trial==1
    size_power_theta_allEpochs = 1;
    size_power_alpha_allEpochs = 1;
    size_power_beta_allEpochs = 1;
    size_power_gamma_allEpochs = 1;
end

% Get the beggining and the end of the trial
title = char(strcat('Trial_', string(trial)));
startTask_times = event_samp(startTask(trial));
endTask_times = event_samp(endTask(trial));

% Get the keypresses within that trial
keypresses_trial = keypresses(keypresses > startTask(trial)...
    & keypresses > endTask(trial));
keypresses_times = event_samp(keypresses_trial);

% Epoch the data into the different keypresses
for epoch = 1:length(keypresses_times)

    EEG_epoch = pop_select(EEG, 'point',...
        [keypresses_times(epoch)-floor(0.4*EEG.srate)...
        keypresses_times(epoch)+ceil(0.4*EEG.srate)]);
    epoch_data = EEG_epoch.data;

    % Using a sliding Hann window
    window = 1:0.8*EEG_epoch.srate;
    
    % Select the data of this specific window [channel x time].
    data_window = epoch_data(:, window);

    % Channel loop
    for channel = 1:size(data_window, 1)
        % Calculate PSD
        [P, f] = periodogram(data_window(channel, :),...
                hann(size(data_window, 2)),...
                2^(2 + nextpow2(size(data_window, 2))),...
                EEG_epoch.srate);

        % Save the power for the frequencies of interest in the
        % different pow variables (all windows will be saved
        % here)
         pow_theta(:, channel) =...
                P((f(:,1)>=4 & f(:,1)<=8),1);
         pow_alpha(:, channel) =...
                P((f(:,1)>=8 & f(:,1)<=13),1);
         pow_beta(:, channel) =...
                P((f(:,1)>=13 & f(:,1)<=32),1);
         pow_gamma(:, channel) =...
                P((f(:,1)>=32 & f(:,1)<=48),1);
         pow_all(:, channel) =...
                P((f(:,1)>=1 & f(:,1)<=48),1);
        end
end

% Change frequency variable for frequencies of interest
freq_theta = f(f(:,1)>=4 & f(:,1)<=8);
freq_alpha = f(f(:,1)>=8 & f(:,1)<=13);
freq_beta = f(f(:,1)>=13 & f(:,1)<=32);
freq_gamma = f(f(:,1)>=32 & f(:,1)<=48);

% Average power over the different frequencies
power_theta_oneEpoch = mean(pow_theta);
power_alpha_oneEpoch = mean(pow_alpha);
power_beta_oneEpoch = mean(pow_beta);
power_gamma_oneEpoch = mean(pow_gamma);

% For the bad channels, give NaN value
for channel = 1:size(data_window, 1)
    if  pow_all(channel) > th(channel)
       power_theta_oneEpoch(channel) = NaN;
       power_alpha_oneEpoch(channel) = NaN;
       power_beta_oneEpoch(channel) = NaN;
       power_gamma_oneEpoch(channel) = NaN;
    end
end

% Add to all epochs array.
power_theta_allEpochs(:, size_power_theta_allEpochs) =...
    power_theta_oneEpoch';
power_alpha_allEpochs(:, size_power_alpha_allEpochs) =...
    power_alpha_oneEpoch';
power_beta_allEpochs(:, size_power_beta_allEpochs) =...
    power_beta_oneEpoch';
power_gamma_allEpochs(:, size_power_gamma_allEpochs) =...
    power_gamma_oneEpoch';

size_power_theta_allEpochs = size(power_theta_allEpochs, 2)+1;
size_power_alpha_allEpochs = size(power_theta_allEpochs, 2)+1;
size_power_beta_allEpochs = size(power_theta_allEpochs, 2)+1;
size_power_gamma_allEpochs = size(power_theta_allEpochs, 2)+1;

end

% Take the average of every epoch
power_theta = mean(power_theta_allEpochs, 2, 'omitnan');
power_alpha = mean(power_alpha_allEpochs, 2, 'omitnan');
power_beta = mean(power_beta_allEpochs, 2, 'omitnan');
power_gamma = mean(power_gamma_allEpochs, 2, 'omitnan');
    
end