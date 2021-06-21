%% Analysis of the EEG signals - separate channels.

clear; clc; close all;
addpath('C:\Users\maria\OneDrive\Documentos\GitHub\Combined-EEG-fNIRS-system\Analysis');

laptop = 'laptopMariana';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);
eeglab;
ft_defaults;

subrec = ["04" "01"];
%%
% Loop through every subject.
for subject = 1:size(subrec, 1)
    sub = subrec(subject, 1);
    rec = subrec(subject, 2);
  
    % Load the subject's EEG signals.
    load([mainpath_in, '\pre-processed\sub-', char(sub), '\eeg\sub-',...
        char(sub), '_rec-', char(rec), '_eeg_divided.mat']);
    
    % Separate into the four different tasks.
    EEG_AutoUncued = EEG_divided.EEG_AutoNoCue;
    EEG_NonAutoUncued = EEG_divided.EEG_NonAutoNoCue;
    EEG_AutoCued = EEG_divided.EEG_AutoCue;
    EEG_NonAutoCued = EEG_divided.EEG_NonAutoCue;
    
    %% Auto Uncued.

    event_samp  = [EEG_AutoUncued.event.latency];
    startTask = find(strcmp({EEG_AutoUncued.event.type}, 's1703')==1);
    endTask = find(strcmp({EEG_AutoUncued.event.type}, 's1711')==1);

    % Get the power spectrum density (PSD) averaged over all trials.
    [power, freq] = calculateAveragePowerAllTrials(EEG_AutoUncued,...
        event_samp, startTask, endTask);
    
    % Save the values onto a allSubjects variable.
    autouncued_power_allSubjects(:, :, subject) = power;
    
    %% Non-Auto Uncued.

    event_samp  = [EEG_NonAutoUncued.event.latency];
    startTask = find(strcmp({EEG_NonAutoUncued.event.type}, 's1705')==1);
    endTask = find(strcmp({EEG_NonAutoUncued.event.type}, 's1713')==1);
    
    % Get the power spectrum density (PSD) averaged over all trials.
    [power, freq] = calculateAveragePowerAllTrials(EEG_NonAutoUncued,...
        event_samp, startTask, endTask);
    
    % Save the values onto a allSubjects variable.
    nonautouncued_power_allSubjects(:, :, subject) = power;
    
    %% Auto Cued.

    event_samp  = [EEG_AutoCued.event.latency];
    startTask = find(strcmp({EEG_AutoCued.event.type}, 's1702')==1);
    endTask = find(strcmp({EEG_AutoCued.event.type}, 's1710')==1);

    % Get the power spectrum density (PSD) averaged over all trials.
    [power, freq] = calculateAveragePowerAllTrials(EEG_AutoCued,...
        event_samp, startTask, endTask);
    
    % Save the values onto a allSubjects variable.
    autocued_power_allSubjects(:, :, subject) = power;
    
    %% Non-Auto Cued.

    event_samp  = [EEG_NonAutoCued.event.latency];
    startTask = find(strcmp({EEG_NonAutoCued.event.type}, 's1704')==1);
    endTask = find(strcmp({EEG_NonAutoCued.event.type}, 's1712')==1);

    % Get the power spectrum density (PSD) averaged over all trials.   
    [power, freq] = calculateAveragePowerAllTrials(EEG_NonAutoCued,...
        event_samp, startTask, endTask);
    
    % Save the values onto a allSubjects variable.
    nonautocued_power_allSubjects(:, :, subject) = power;

    disp(['These are the results for subject ', char(sub), '.']);
    disp('Press any key to move onto the next subject.');
    pause;
    close all;
    
end

% Get the power spectrum density (PSD) averaged over all subjects.
% Auto Uncued.
autouncued_power = mean(autouncued_power_allSubjects, 3);
% Non-Auto Uncued.
nonautouncued_power = mean(nonautouncued_power_allSubjects, 3);
% Auto Cued.
autocued_power = mean(autocued_power_allSubjects, 3);
% Non-Auto Cued.
nonautocued_power = mean(nonautocued_power_allSubjects, 3);

disp('This was the end of individual subjects.');
disp('These are the results for the average of all subjects.');

%% Functions

% Loop through the power from the individual trials and average them.
function [power, freq] = calculateAveragePowerAllTrials(EEG, event_samp,...
    startTask, endTask)

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

% From the trial data, calculate the power over the frequencies of the signal
% for all electrodes.
function [power, freq] =...
    calculatePowerPerTrial(EEG_trial, trial_data)
    
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