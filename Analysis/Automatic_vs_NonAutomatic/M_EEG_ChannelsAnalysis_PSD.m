%% Analysis of the EEG signals - separate channels.

clear; clc; close all;
addpath('C:\Users\maria\OneDrive\Documentos\GitHub\Combined-EEG-fNIRS-system\Analysis');

laptop = 'laptopMariana';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);
eeglab;
ft_defaults;
results_path = 'C:\Users\maria\OneDrive\Ambiente de Trabalho\Automaticity Results\Separate Channels';

subrec = ["28" "04"];

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

    %% Get the locations of the channels of interest.
    
    locs = {EEG_AutoCued.chanlocs.labels};
    F7_loc = find(contains(locs, 'F7'));
    F8_loc = find(contains(locs, 'F8'));
    AFFz_loc = find(contains(locs, 'AFFz'));
    FC1_loc = find(contains(locs, 'FC1'));
    
    %% Get the values of power for the specific channels and plot them.
    
    % F7.
    autouncued_power_F7 = autouncued_power_allSubjects(:, F7_loc, subject);
    autocued_power_F7 = autocued_power_allSubjects(:, F7_loc, subject);
    nonautouncued_power_F7 = nonautouncued_power_allSubjects(:, F7_loc, subject);
    nonautocued_power_F7 = nonautocued_power_allSubjects(:, F7_loc, subject);
    
    figure; title('F7');
    plot(freq, autouncued_power_F7); hold on;
    plot(freq, autocued_power_F7); hold on;
    plot(freq, nonautouncued_power_F7); hold on;
    plot(freq, nonautocued_power_F7); hold on;
    xline(4); hold on;
    xline(8); hold on;
    xline(13); hold on;
    xline(32); hold off;
    legend('Auto Uncued','Auto Cued', 'Non-Auto Uncued', 'Non-Auto Cued');
    
    % F8.
    autouncued_power_F8 = autouncued_power_allSubjects(:, F8_loc, subject);
    autocued_power_F8 = autocued_power_allSubjects(:, F8_loc, subject);
    nonautouncued_power_F8 = nonautouncued_power_allSubjects(:, F8_loc, subject);
    nonautocued_power_F8 = nonautocued_power_allSubjects(:, F8_loc, subject);
    
    figure; title('F8');
    plot(freq, autouncued_power_F8); hold on;
    plot(freq, autocued_power_F8); hold on;
    plot(freq, nonautouncued_power_F8); hold on;
    plot(freq, nonautocued_power_F8); hold on;
    xline(4); hold on;
    xline(8); hold on;
    xline(13); hold on;
    xline(32); hold off;
    legend('Auto Uncued','Auto Cued', 'Non-Auto Uncued', 'Non-Auto Cued');
    
    % AFFz.
    autouncued_power_AFFz = autouncued_power_allSubjects(:, AFFz_loc, subject);
    autocued_power_AFFz = autocued_power_allSubjects(:, AFFz_loc, subject);
    nonautouncued_power_AFFz = nonautouncued_power_allSubjects(:, AFFz_loc, subject);
    nonautocued_power_AFFz = nonautocued_power_allSubjects(:, AFFz_loc, subject);
    
    figure; title('AFFz');
    plot(freq, autouncued_power_AFFz); hold on;
    plot(freq, autocued_power_AFFz); hold on;
    plot(freq, nonautouncued_power_AFFz); hold on;
    plot(freq, nonautocued_power_AFFz); hold on;
    xline(4); hold on;
    xline(8); hold on;
    xline(13); hold on;
    xline(32); hold off;
    legend('Auto Uncued','Auto Cued', 'Non-Auto Uncued', 'Non-Auto Cued');
    
    % FC1.
    autouncued_power_FC1 = autouncued_power_allSubjects(:, FC1_loc, subject);
    autocued_power_FC1 = autocued_power_allSubjects(:, FC1_loc, subject);
    nonautouncued_power_FC1 = nonautouncued_power_allSubjects(:, FC1_loc, subject);
    nonautocued_power_FC1 = nonautocued_power_allSubjects(:, FC1_loc, subject);
    
    figure; title('FC1');
    plot(freq, autouncued_power_FC1); hold on;
    plot(freq, autocued_power_FC1); hold on;
    plot(freq, nonautouncued_power_FC1); hold on;
    plot(freq, nonautocued_power_FC1); hold on;
    xline(4); hold on;
    xline(8); hold on;
    xline(13); hold on;
    xline(32); hold off;
    legend('Auto Uncued','Auto Cued', 'Non-Auto Uncued', 'Non-Auto Cued');
    
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

%% Plot the PSD for specific channels.

% F7.
autouncued_power_F7 = autouncued_power(:, F7_loc);
autocued_power_F7 = autocued_power(:, F7_loc);
nonautouncued_power_F7 = nonautouncued_power(:, F7_loc);
nonautocued_power_F7 = nonautocued_power(:, F7_loc);

figure; title('F7');
plot(freq, autouncued_power_F7); hold on;
plot(freq, autocued_power_F7); hold on;
plot(freq, nonautouncued_power_F7); hold on;
plot(freq, nonautocued_power_F7); hold on;
xline(4); hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
legend('Auto Uncued','Auto Cued', 'Non-Auto Uncued', 'Non-Auto Cued');

saveas(gcf, fullfile(results_path, 'F7_PowervsFreq'),'png');

% F8.
autouncued_power_F8 = autouncued_power(:, F8_loc);
autocued_power_F8 = autocued_power(:, F8_loc);
nonautouncued_power_F8 = nonautouncued_power(:, F8_loc);
nonautocued_power_F8 = nonautocued_power(:, F8_loc);

figure; title('F8');
plot(freq, autouncued_power_F8); hold on;
plot(freq, autocued_power_F8); hold on;
plot(freq, nonautouncued_power_F8); hold on;
plot(freq, nonautocued_power_F8); hold on;
xline(4); hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
legend('Auto Uncued','Auto Cued', 'Non-Auto Uncued', 'Non-Auto Cued');

saveas(gcf, fullfile(results_path, 'F8_PowervsFreq'),'png');

% AFFz.
autouncued_power_AFFz = autouncued_power(:, AFFz_loc);
autocued_power_AFFz = autocued_power(:, AFFz_loc);
nonautouncued_power_AFFz = nonautouncued_power(:, AFFz_loc);
nonautocued_power_AFFz = nonautocued_power(:, AFFz_loc);

figure; title('AFFz');
plot(freq, autouncued_power_AFFz); hold on;
plot(freq, autocued_power_AFFz); hold on;
plot(freq, nonautouncued_power_AFFz); hold on;
plot(freq, nonautocued_power_AFFz); hold on;
xline(4); hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
legend('Auto Uncued','Auto Cued', 'Non-Auto Uncued', 'Non-Auto Cued');

saveas(gcf, fullfile(results_path, 'AFFz_PowervsFreq'),'png');

% FC1.
autouncued_power_FC1 = autouncued_power(:, FC1_loc);
autocued_power_FC1 = autocued_power(:, FC1_loc);
nonautouncued_power_FC1 = nonautouncued_power(:, FC1_loc);
nonautocued_power_FC1 = nonautocued_power(:, FC1_loc);

figure; title('FC1');
plot(freq, autouncued_power_FC1); hold on;
plot(freq, autocued_power_FC1); hold on;
plot(freq, nonautouncued_power_FC1); hold on;
plot(freq, nonautocued_power_FC1); hold on;
xline(4); hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
legend('Auto Uncued','Auto Cued', 'Non-Auto Uncued', 'Non-Auto Cued');

saveas(gcf, fullfile(results_path, 'FC1_PowervsFreq'),'png');

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