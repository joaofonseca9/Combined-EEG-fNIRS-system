function th = calculateThreshold(EEG_divided, sub)

EEG_DualUncued = EEG_divided.EEG_NonAutoDualNoCue;
EEG_SingleUncued = EEG_divided.EEG_NonAutoNoCue;
EEG_DualCued = EEG_divided.EEG_NonAutoDualCue;
EEG_SingleCued = EEG_divided.EEG_NonAutoCue;
    
size_power = 1;

%% Dual Uncued
event_samp  = [EEG_DualUncued.event.latency];
startBaseline = find(strcmp({EEG_DualUncued.event.type}, 'boundary')==1);
startTask = find(strcmp({EEG_DualUncued.event.type}, 's1709')==1);
endTask = find(strcmp({EEG_DualUncued.event.type}, 's1717')==1);

if sub=="64"
    startBaseline = [2 36 70];
else
    for i=length(startBaseline):-1:1
        if ~any(startBaseline(i)+1==startTask(:)) 
            startBaseline(i) = [];
        end
    end
end
    

for trial=1:length(startTask)
    
    startBaseline_times = event_samp(startBaseline(trial));
    endTask_times = event_samp(endTask(trial));
    
    EEG_trial = pop_select(EEG_DualUncued, 'point',...
        [startBaseline_times endTask_times]);
    trial_data = EEG_trial.data;
    
    % Using a sliding Hann window
    window_id = 1;
    window = 1:1*EEG_trial.srate;
    while window(end) <= size(trial_data, 2)
        % Select the data of this specific window [channel x time]
        data_window = trial_data(:, window);
        
        % Channel loop
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
        % 0.5*fs)
        window_id = window_id + 1;
        window = window+0.5*EEG_trial.srate;
    end
    
    % Average power per channel over windows
    power = mean(pow, 3, 'omitnan');
   
    power_all(:, :, size_power) = power;
    size_power = size(power_all, 3)+1;
end

%% Single Uncued
event_samp  = [EEG_SingleUncued.event.latency];
startBaseline = find(strcmp({EEG_SingleUncued.event.type}, 'boundary')==1);
startTask = find(strcmp({EEG_SingleUncued.event.type}, 's1705')==1);
endTask = find(strcmp({EEG_SingleUncued.event.type}, 's1713')==1);
if sub=="64"
    startBaseline = [27 48];
else
    for i=length(startBaseline):-1:1
        if ~any(startBaseline(i)+1==startTask(:)) 
            startBaseline(i) = [];
        end
    end
end

for trial=1:length(startTask)
    
    startBaseline_times = event_samp(startBaseline(trial));
    endTask_times = event_samp(endTask(trial));
    
    EEG_trial = pop_select(EEG_SingleUncued, 'point',...
        [startBaseline_times endTask_times]);
    trial_data = EEG_trial.data;
    
    % Using a sliding Hann window
    window_id = 1;
    window = 1:1*EEG_trial.srate;
    while window(end) <= size(trial_data, 2)
        % Select the data of this specific window [channel x time]
        data_window = trial_data(:, window);
        
        % Channel loop
        for channel = 1:size(data_window, 1)
            % If window is NOT removed because of badchannel (=NaN)
            if isempty(find(isnan(data_window(channel, :))))
                % Calculate PSD
                [P, f] = periodogram(data_window(channel, :),...
                    hann(size(data_window, 2)),...
                    2^(2 + nextpow2(size(data_window, 2))), EEG_trial.srate);
                % Save the power for frequencies in between 1 and 48
                % Hz
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

%% Dual Cued
event_samp  = [EEG_DualCued.event.latency];
startBaseline = find(strcmp({EEG_DualCued.event.type}, 'boundary')==1);
startTask = find(strcmp({EEG_DualCued.event.type}, 's1708')==1);
endTask = find(strcmp({EEG_DualCued.event.type}, 's1701')==1);
for i=length(startBaseline):-1:1
    if startBaseline(i)+1==startTask(:)
        startBaseline(i) = [];
    end
end

for trial=1:length(startTask)
    
    startBaseline_times = event_samp(startBaseline(trial));
    endTask_times = event_samp(endTask(trial));
    
    EEG_trial = pop_select(EEG_DualCued, 'point',...
        [startBaseline_times endTask_times]);
    trial_data = EEG_trial.data;
    
    % Using a sliding Hann window
    window_id = 1;
    window = 1:1*EEG_trial.srate;
    while window(end) <= size(trial_data, 2)
        % Select the data of this specific window [channel x time]
        data_window = trial_data(:, window);
        
        % Channel loop
        for channel = 1:size(data_window, 1)
            % If window is NOT removed because of badchannel (=NaN)
            if isempty(find(isnan(data_window(channel, :))))
                % Calculate PSD
                [P, f] = periodogram(data_window(channel, :),...
                    hann(size(data_window, 2)),...
                    2^(2 + nextpow2(size(data_window, 2))), EEG_trial.srate);
                % Save the power for frequencies in between 1 and 48
                % Hz
                pow(:, channel, window_id) = P((f(:,1)>=1 & f(:,1)<=48),1);
            else
                pow(:, channel, window_id) = NaN;
            end
        end
        % Increase indices and window (increase sliding window with
        % 0.5*fs)
        window_id = window_id + 1;
        window = window+0.5*EEG_trial.srate;
    end
    
    % Average power per channel over windows
    power = mean(pow, 3, 'omitnan');
    
    power_all(:, :, size_power) = power;
    size_power = size(power_all, 3)+1;
end

%% Single Cued
event_samp  = [EEG_SingleCued.event.latency];
startBaseline = find(strcmp({EEG_SingleCued.event.type}, 'boundary')==1);
startTask = find(strcmp({EEG_SingleCued.event.type}, 's1704')==1);
endTask = find(strcmp({EEG_SingleCued.event.type}, 's1701')==1);
if sub=="64"
    startBaseline = [162];
    startTask = [164];
else
    for i=length(startBaseline):-1:1
        if ~any(startBaseline(i)+1==startTask(:)) 
            startBaseline(i) = [];
        end
    end
end

for trial=1:length(startTask)
    
    startBaseline_times = event_samp(startBaseline(trial));
    endTask_times = event_samp(endTask(trial));
    
    EEG_trial = pop_select(EEG_SingleCued, 'point',...
        [startBaseline_times endTask_times]);
    trial_data = EEG_trial.data;
    
    % Using a sliding Hann window
    window_id = 1;
    window = 1:1*EEG_trial.srate;
    while window(end) <= size(trial_data, 2)
        % Select the data of this specific window [channel x time]
        data_window = trial_data(:, window);
        
        % Channel loop
        for channel = 1:size(data_window, 1)
            % If window is NOT removed because of badchannel (=NaN)
            if isempty(find(isnan(data_window(channel, :))))
                % Calculate PSD
                [P, f] = periodogram(data_window(channel, :),...
                    hann(size(data_window, 2)),...
                    2^(2 + nextpow2(size(data_window, 2))), EEG_trial.srate);
                % Save the power for frequencies in between 1 and 48
                % Hz
                pow(:, channel, window_id) = P((f(:,1)>=1 & f(:,1)<=48),1);
            else
                pow(:, channel, window_id) = NaN;
            end
        end
        % Increase indices and window (increase sliding window with
        % 0.5*fs)
        window_id = window_id + 1;
        window = window+0.5*EEG_trial.srate;
    end
    
    % Average power per channel over windows
    power = mean(pow, 3, 'omitnan');
    
    power_all(:, :, size_power) = power;
    size_power = size(power_all, 3)+1;
end

%% Final value for threshold

% Average over frequencies
power_all = squeeze(mean(power_all, 1));

% Average over frequencies to get the mean of the power for every channel
mean_power = mean(power_all, 2);

% Same thing for std
std_power = std(power_all, 0, 2);

% Calculate the threshold
th = mean_power + 3*std_power;

end