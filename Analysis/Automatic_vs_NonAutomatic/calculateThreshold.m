function th = calculateThreshold(EEG_divided)

EEG_AutoUncued = EEG_divided.EEG_AutoNoCue;
EEG_NonAutoUncued = EEG_divided.EEG_NonAutoNoCue;
EEG_AutoCued = EEG_divided.EEG_AutoCue;
EEG_NonAutoCued = EEG_divided.EEG_NonAutoCue;

size_power = 1;

%% Non-Auto Uncued.
event_samp  = [EEG_NonAutoUncued.event.latency];
startBaseline = find(strcmp({EEG_NonAutoUncued.event.type}, 'boundary')==1);
startTask = find(strcmp({EEG_NonAutoUncued.event.type}, 's1705')==1);
endTask = find(strcmp({EEG_NonAutoUncued.event.type}, 's1713')==1);
for i=length(startBaseline):-1:1
    if ~any(startBaseline(i)+1==startTask(:))
        startBaseline(i) = [];
    end
end

for trial=1:length(startTask)
    
    startBaseline_times = event_samp(startBaseline(trial));
    endTask_times = event_samp(endTask(trial));
    
    EEG_trial = pop_select(EEG_NonAutoUncued, 'point',...
        [startBaseline_times endTask_times]);
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
                [P, f] = periodogram(data_window(channel, :),...
                    hann(size(data_window, 2)),...
                    2^(2 + nextpow2(size(data_window, 2))), EEG_trial.srate);
                % Save the power for frequencies in between 1 and 48
                % Hz.
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
    
    % Average power per channel over windows.
    power = mean(pow, 3, 'omitnan');
    
    power_all(:, :, size_power) = power;
    size_power = size(power_all, 3)+1;
end

%% Auto Uncued.
event_samp  = [EEG_AutoUncued.event.latency];
startBaseline = find(strcmp({EEG_AutoUncued.event.type}, 'boundary')==1);
startTask = find(strcmp({EEG_AutoUncued.event.type}, 's1703')==1);
endTask = find(strcmp({EEG_AutoUncued.event.type}, 's1711')==1);
for i=length(startBaseline):-1:1
    if ~any(startBaseline(i)+1==startTask(:))
        startBaseline(i) = [];
    end
end

for trial=1:length(startTask)
    
    startBaseline_times = event_samp(startBaseline(trial));
    endTask_times = event_samp(endTask(trial));
    
    EEG_trial = pop_select(EEG_AutoUncued, 'point',...
        [startBaseline_times endTask_times]);
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
                [P, f] = periodogram(data_window(channel, :),...
                    hann(size(data_window, 2)),...
                    2^(2 + nextpow2(size(data_window, 2))), EEG_trial.srate);
                % Save the power for frequencies in between 1 and 48
                % Hz.
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
    
    % Average power per channel over windows.
    power = mean(pow, 3, 'omitnan');
    
    power_all(:, :, size_power) = power;
    size_power = size(power_all, 3)+1;
end

%% Non-Auto Cued.
event_samp  = [EEG_NonAutoCued.event.latency];
startBaseline = find(strcmp({EEG_NonAutoCued.event.type}, 'boundary')==1);
startTask = find(strcmp({EEG_NonAutoCued.event.type}, 's1704')==1);
endTask = find(strcmp({EEG_NonAutoCued.event.type}, 's1712')==1);
for i=length(startBaseline):-1:1
    if ~any(startBaseline(i)+1==startTask(:))
        startBaseline(i) = [];
    end
end

for trial=1:length(startTask)
    
    startBaseline_times = event_samp(startBaseline(trial));
    endTask_times = event_samp(endTask(trial));
    
    EEG_trial = pop_select(EEG_NonAutoCued, 'point',...
        [startBaseline_times endTask_times]);
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
                [P, f] = periodogram(data_window(channel, :),...
                    hann(size(data_window, 2)),...
                    2^(2 + nextpow2(size(data_window, 2))), EEG_trial.srate);
                % Save the power for frequencies in between 1 and 48
                % Hz.
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
    
    % Average power per channel over windows.
    power = mean(pow, 3, 'omitnan');
    
    power_all(:, :, size_power) = power;
    size_power = size(power_all, 3)+1;
end

%% Auto Cued.
event_samp  = [EEG_AutoCued.event.latency];
startBaseline = find(strcmp({EEG_AutoCued.event.type}, 'boundary')==1);
startTask = find(strcmp({EEG_AutoCued.event.type}, 's1702')==1);
endTask = find(strcmp({EEG_AutoCued.event.type}, 's1710')==1);
for i=length(startBaseline):-1:1
    if ~any(startBaseline(i)+1==startTask(:))
        startBaseline(i) = [];
    end
end

for trial=1:length(startTask)
    
    startBaseline_times = event_samp(startBaseline(trial));
    endTask_times = event_samp(endTask(trial));
    
    EEG_trial = pop_select(EEG_AutoCued, 'point',...
        [startBaseline_times endTask_times]);
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
                [P, f] = periodogram(data_window(channel, :),...
                    hann(size(data_window, 2)),...
                    2^(2 + nextpow2(size(data_window, 2))), EEG_trial.srate);
                % Save the power for frequencies in between 1 and 48
                % Hz.
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
    
    % Average power per channel over windows.
    power = mean(pow, 3, 'omitnan');
    
    power_all(:, :, size_power) = power;
    size_power = size(power_all, 3)+1;
end

%% Final value for threshold.

% Average over frequencies.
power_all = squeeze(mean(power_all, 1));

% Average over frequencies to get the mean of the power for every channel.
% Same thing for std.
mean_power = mean(power_all, 2);
std_power = std(power_all, 0, 2);

% Calculate the threshold.
th = mean_power + 3*std_power;

end