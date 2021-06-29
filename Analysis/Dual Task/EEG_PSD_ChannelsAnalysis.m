%% Analysis of the EEG signal (PSD: separate channels)
clear; clc; close all;

%% Initialize data
addpath('C:\Users\catar\OneDrive - Universidade do Porto\Twente\Combined-EEG-fNIRS-system\Analysis');
addpath('C:\Users\catar\OneDrive - Universidade do Porto\Twente\Combined-EEG-fNIRS-system\Analysis\Dual Task');
laptop = 'laptopCatarina';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);
results_path = 'C:\Users\catar\OneDrive - Universidade do Porto\Twente\Data Analysis\eeg\psd';

eeglab;

subrec = ["28" "04"];

% List of the 30 channels in the cap
list_channels = ["Fp1"; "Fpz"; "Fp2"; "F7"; "F3"; "AFFz"; "F4"; "F8";...
    "FC5"; "FC1"; "FC2"; "FC6"; "T7"; "C3"; "Cz"; "C4"; "T8"; "CP5";...
    "CP1"; "CP2"; "CP6"; "P7"; "P3"; "Pz"; "P4"; "P8"; "POz"; "O1";...
    "Oz"; "O2"];

%% Load data + processing per subject
% Go through all subjects
for subject = 1:size(subrec, 1)
    sub = subrec(subject, 1);
    rec = subrec(subject, 2);
  
    % Load EEG preprocessed signal
    load([mainpath_in, '\pre-processed\sub-', char(sub), '\eeg\sub-',...
        char(sub), '_rec-', char(rec), '_eeg_divided.mat'], 'EEG_divided');
    
    % Separate into the four different tasks
    EEG_DualUncued = EEG_divided.EEG_NonAutoDualNoCue;
    EEG_SingleUncued = EEG_divided.EEG_NonAutoNoCue;
    EEG_DualCued = EEG_divided.EEG_NonAutoDualCue;
    EEG_SingleCued = EEG_divided.EEG_NonAutoCue;
    
    %% Dual Uncued: PSD
    event_samp  = [EEG_DualUncued.event.latency];
    startTask = find(strcmp({EEG_DualUncued.event.type}, 's1798')==1);
    endTask = find(strcmp({EEG_DualUncued.event.type}, 's1717')==1);

    % Get the PSD averaged over all trials
    [power, freq] = calculateAveragePowerAllTrials(EEG_DualUncued,...
        event_samp, startTask, endTask);
    
    % Compensate for removed channels
    power = compensateRemovedChannels(power, EEG_DualUncued, list_channels);
    
    % Save the values onto a allSubjects variable
    dualuncued_power_allSubjects(:, :, subject) = power;
    
    %% Single Uncued: PSD
    event_samp  = [EEG_SingleUncued.event.latency];
    startTask = find(strcmp({EEG_SingleUncued.event.type}, 's1705')==1);
    endTask = find(strcmp({EEG_SingleUncued.event.type}, 's1713')==1);
    
    % Get the PSD averaged over all trials
    [power, freq] = calculateAveragePowerAllTrials(EEG_SingleUncued,...
        event_samp, startTask, endTask);
    
    % Compensate for removed channels
    power = compensateRemovedChannels(power, EEG_SingleUncued, list_channels);
    
    % Save the values onto a allSubjects variable
    singleuncued_power_allSubjects(:, :, subject) = power;
    
    %% Dual Cued: PSD
    event_samp  = [EEG_DualCued.event.latency];
    startTask = find(strcmp({EEG_DualCued.event.type}, 's1798')==1);
    endTask = find(strcmp({EEG_DualCued.event.type}, 's1701')==1);

    % Get the PSD averaged over all trials
    [power, freq] = calculateAveragePowerAllTrials(EEG_DualCued,...
        event_samp, startTask, endTask);
    
    % Compensate for removed channels
    power = compensateRemovedChannels(power, EEG_DualCued, list_channels);
    
    % Save the values onto a allSubjects variable
    dualcued_power_allSubjects(:, :, subject) = power;
    
    %% Single Cued: PSD
    event_samp  = [EEG_SingleCued.event.latency];
    startTask = find(strcmp({EEG_SingleCued.event.type}, 's1704')==1);
    endTask = find(strcmp({EEG_SingleCued.event.type}, 's1701')==1);

    % Get the power PSD averaged over all trials
    [power, freq] = calculateAveragePowerAllTrials(EEG_SingleCued,...
        event_samp, startTask, endTask);
    
    % Compensate for removed channels
    power = compensateRemovedChannels(power, EEG_SingleCued, list_channels);
    
    % Save the values onto a allSubjects variable
    singlecued_power_allSubjects(:, :, subject) = power;

    %% Get the locations of the channels of interest
    locs = {EEG_DualCued.chanlocs.labels};
    % PFC:
    F7_loc = find(contains(locs, 'F7'));
    F8_loc = find(contains(locs, 'F8'));
    % SMA:
    FC1_loc = find(contains(locs, 'FC1'));
    FC2_loc = find(contains(locs, 'FC2'));
    % M1:
    C3_loc = find(contains(locs, 'C3'));
    CP1_loc = find(contains(locs, 'CP1'));
    % PPC:
    P3_loc = find(contains(locs, 'P3'));
    P4_loc = find(contains(locs, 'P4'));
    
    %% Get the values of power for the specific channels and plot them
    % PFC: F7
    dualuncued_power_F7 = dualuncued_power_allSubjects(:, F7_loc, subject);
    dualcued_power_F7 = dualcued_power_allSubjects(:, F7_loc, subject);
    singleuncued_power_F7 = singleuncued_power_allSubjects(:, F7_loc, subject);
    singlecued_power_F7 = singlecued_power_allSubjects(:, F7_loc, subject);
    
    figure; title('F7');
    plot(freq, dualuncued_power_F7); hold on;
    plot(freq, dualcued_power_F7); hold on;
    plot(freq, singleuncued_power_F7); hold on;
    plot(freq, singlecued_power_F7); hold on;
    xline(4); hold on;
    xline(8); hold on;
    xline(13); hold on;
    xline(32); hold off;
    legend('Dual Uncued','Dual Cued', 'Single Uncued', 'Single Cued');
    
    % PFC: F8
    dualuncued_power_F8 = dualuncued_power_allSubjects(:, F8_loc, subject);
    dualcued_power_F8 = dualcued_power_allSubjects(:, F8_loc, subject);
    singleuncued_power_F8 = singleuncued_power_allSubjects(:, F8_loc, subject);
    singlecued_power_F8 = singlecued_power_allSubjects(:, F8_loc, subject);
    
    figure; title('F8');
    plot(freq, dualuncued_power_F8); hold on;
    plot(freq, dualcued_power_F8); hold on;
    plot(freq, singleuncued_power_F8); hold on;
    plot(freq, singlecued_power_F8); hold on;
    xline(4); hold on;
    xline(8); hold on;
    xline(13); hold on;
    xline(32); hold off;
    legend('Dual Uncued','Dual Cued', 'Single Uncued', 'Single Cued');
    
    % SMA: FC1
    dualuncued_power_FC1 = dualuncued_power_allSubjects(:, FC1_loc, subject);
    dualcued_power_FC1 = dualcued_power_allSubjects(:, FC1_loc, subject);
    singleuncued_power_FC1 = singleuncued_power_allSubjects(:, FC1_loc, subject);
    singlecued_power_FC1 = singlecued_power_allSubjects(:, FC1_loc, subject);
    
    figure; title('FC1');
    plot(freq, dualuncued_power_FC1); hold on;
    plot(freq, dualcued_power_FC1); hold on;
    plot(freq, singleuncued_power_FC1); hold on;
    plot(freq, singlecued_power_FC1); hold on;
    xline(4); hold on;
    xline(8); hold on;
    xline(13); hold on;
    xline(32); hold off;
    legend('Dual Uncued','Dual Cued', 'Single Uncued', 'Single Cued');
    
    % SMA: FC2
    dualuncued_power_FC2 = dualuncued_power_allSubjects(:, FC2_loc, subject);
    dualcued_power_FC2 = dualcued_power_allSubjects(:, FC2_loc, subject);
    singleuncued_power_FC2 = singleuncued_power_allSubjects(:, FC2_loc, subject);
    singlecued_power_FC2 = singlecued_power_allSubjects(:, FC2_loc, subject);
    
    figure; title('FC2');
    plot(freq, dualuncued_power_FC2); hold on;
    plot(freq, dualcued_power_FC2); hold on;
    plot(freq, singleuncued_power_FC2); hold on;
    plot(freq, singlecued_power_FC2); hold on;
    xline(4); hold on;
    xline(8); hold on;
    xline(13); hold on;
    xline(32); hold off;
    legend('Dual Uncued','Dual Cued', 'Single Uncued', 'Single Cued');
    
    % M1: C3
    dualuncued_power_C3 = dualuncued_power_allSubjects(:, C3_loc, subject);
    dualcued_power_C3 = dualcued_power_allSubjects(:, C3_loc, subject);
    singleuncued_power_C3 = singleuncued_power_allSubjects(:, C3_loc, subject);
    singlecued_power_C3 = singlecued_power_allSubjects(:, C3_loc, subject);
    
    figure; title('C3');
    plot(freq, dualuncued_power_C3); hold on;
    plot(freq, dualcued_power_C3); hold on;
    plot(freq, singleuncued_power_C3); hold on;
    plot(freq, singlecued_power_C3); hold on;
    xline(4); hold on;
    xline(8); hold on;
    xline(13); hold on;
    xline(32); hold off;
    legend('Dual Uncued','Dual Cued', 'Single Uncued', 'Single Cued');
    
    % M1: CP1
    dualuncued_power_CP1 = dualuncued_power_allSubjects(:, CP1_loc, subject);
    dualcued_power_CP1 = dualcued_power_allSubjects(:, CP1_loc, subject);
    singleuncued_power_CP1 = singleuncued_power_allSubjects(:, CP1_loc, subject);
    singlecued_power_CP1 = singlecued_power_allSubjects(:, CP1_loc, subject);
    
    figure; title('CP1');
    plot(freq, dualuncued_power_CP1); hold on;
    plot(freq, dualcued_power_CP1); hold on;
    plot(freq, singleuncued_power_CP1); hold on;
    plot(freq, singlecued_power_CP1); hold on;
    xline(4); hold on;
    xline(8); hold on;
    xline(13); hold on;
    xline(32); hold off;
    legend('Dual Uncued','Dual Cued', 'Single Uncued', 'Single Cued');
    
    % PPC: P3
    dualuncued_power_P3 = dualuncued_power_allSubjects(:, P3_loc, subject);
    dualcued_power_P3 = dualcued_power_allSubjects(:, P3_loc, subject);
    singleuncued_power_P3 = singleuncued_power_allSubjects(:, P3_loc, subject);
    singlecued_power_P3 = singlecued_power_allSubjects(:, P3_loc, subject);
    
    figure; title('P3');
    plot(freq, dualuncued_power_P3); hold on;
    plot(freq, dualcued_power_P3); hold on;
    plot(freq, singleuncued_power_P3); hold on;
    plot(freq, singlecued_power_P3); hold on;
    xline(4); hold on;
    xline(8); hold on;
    xline(13); hold on;
    xline(32); hold off;
    legend('Dual Uncued','Dual Cued', 'Single Uncued', 'Single Cued');
    
    % PPC: P4
    dualuncued_power_P4 = dualuncued_power_allSubjects(:, P4_loc, subject);
    dualcued_power_P4 = dualcued_power_allSubjects(:, P4_loc, subject);
    singleuncued_power_P4 = singleuncued_power_allSubjects(:, P4_loc, subject);
    singlecued_power_P4 = singlecued_power_allSubjects(:, P4_loc, subject);
    
    figure; title('P4');
    plot(freq, dualuncued_power_P4); hold on;
    plot(freq, dualcued_power_P4); hold on;
    plot(freq, singleuncued_power_P4); hold on;
    plot(freq, singlecued_power_P4); hold on;
    xline(4); hold on;
    xline(8); hold on;
    xline(13); hold on;
    xline(32); hold off;
    legend('Dual Uncued','Dual Cued', 'Single Uncued', 'Single Cued');
    
    disp(['These are the results for subject ', char(sub), '.']);
    disp('Press any key to move onto the next subject.');
    pause;
    close all;
    
end

%% Get the PSD averaged over all subjects
% Dual Cued
dualcued_power = mean(dualcued_power_allSubjects, 3, 'omitnan');
% Dual Uncued
dualuncued_power = mean(dualuncued_power_allSubjects, 3, 'omitnan');
% Single Cued
singlecued_power = mean(singlecued_power_allSubjects, 3, 'omitnan');
% Single Uncued
singleuncued_power = mean(singleuncued_power_allSubjects, 3, 'omitnan');

%% Plot the PSD for specific channels
% PFC: F7
dualuncued_power_F7 = dualuncued_power(:, F7_loc);
dualcued_power_F7 = dualcued_power(:, F7_loc);
singleuncued_power_F7 = singleuncued_power(:, F7_loc);
singlecued_power_F7 = singlecued_power(:, F7_loc);

figure; title('F7');
plot(freq, dualuncued_power_F7); hold on;
plot(freq, dualcued_power_F7); hold on;
plot(freq, singleuncued_power_F7); hold on;
plot(freq, singlecued_power_F7); hold on;
xline(4); hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
legend('Dual Uncued','Dual Cued', 'Single Uncued', 'Single Cued');

saveas(gcf, fullfile(results_path, 'F7_PowervsFreq'),'png');

% PFC: F8
dualuncued_power_F8 = dualuncued_power(:, F8_loc);
dualcued_power_F8 = dualcued_power(:, F8_loc);
singleuncued_power_F8 = singleuncued_power(:, F8_loc);
singlecued_power_F8 = singlecued_power(:, F8_loc);

figure; title('F8');
plot(freq, dualuncued_power_F8); hold on;
plot(freq, dualcued_power_F8); hold on;
plot(freq, singleuncued_power_F8); hold on;
plot(freq, singlecued_power_F8); hold on;
xline(4); hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
legend('Dual Uncued','Dual Cued', 'Single Uncued', 'Single Cued');

saveas(gcf, fullfile(results_path, 'F8_PowervsFreq'),'png');

% SMA: FC1
dualuncued_power_FC1 = dualuncued_power(:, FC1_loc);
dualcued_power_FC1 = dualcued_power(:, FC1_loc);
singleuncued_power_FC1 = singleuncued_power(:, FC1_loc);
singlecued_power_FC1 = singlecued_power(:, FC1_loc);

figure; title('FC1');
plot(freq, dualuncued_power_FC1); hold on;
plot(freq, dualcued_power_FC1); hold on;
plot(freq, singleuncued_power_FC1); hold on;
plot(freq, singlecued_power_FC1); hold on;
xline(4); hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
legend('Dual Uncued','Dual Cued', 'Single Uncued', 'Single Cued');

saveas(gcf, fullfile(results_path, 'FC1_PowervsFreq'),'png');

% SMA: FC2
dualuncued_power_FC2 = dualuncued_power(:, FC2_loc);
dualcued_power_FC2 = dualcued_power(:, FC2_loc);
singleuncued_power_FC2 = singleuncued_power(:, FC2_loc);
singlecued_power_FC2 = singlecued_power(:, FC2_loc);

figure; title('FC2');
plot(freq, dualuncued_power_FC2); hold on;
plot(freq, dualcued_power_FC2); hold on;
plot(freq, singleuncued_power_FC2); hold on;
plot(freq, singlecued_power_FC2); hold on;
xline(4); hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
legend('Dual Uncued','Dual Cued', 'Single Uncued', 'Single Cued');

saveas(gcf, fullfile(results_path, 'FC2_PowervsFreq'),'png');

% M1: C3
dualuncued_power_C3 = dualuncued_power(:, C3_loc);
dualcued_power_C3 = dualcued_power(:, C3_loc);
singleuncued_power_C3 = singleuncued_power(:, C3_loc);
singlecued_power_C3 = singlecued_power(:, C3_loc);

figure; title('C3');
plot(freq, dualuncued_power_C3); hold on;
plot(freq, dualcued_power_C3); hold on;
plot(freq, singleuncued_power_C3); hold on;
plot(freq, singlecued_power_C3); hold on;
xline(4); hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
legend('Dual Uncued','Dual Cued', 'Single Uncued', 'Single Cued');

saveas(gcf, fullfile(results_path, 'C3_PowervsFreq'),'png');

% M1: CP1
dualuncued_power_CP1 = dualuncued_power(:, CP1_loc);
dualcued_power_CP1 = dualcued_power(:, CP1_loc);
singleuncued_power_CP1 = singleuncued_power(:, CP1_loc);
singlecued_power_CP1 = singlecued_power(:, CP1_loc);

figure; title('CP1');
plot(freq, dualuncued_power_CP1); hold on;
plot(freq, dualcued_power_CP1); hold on;
plot(freq, singleuncued_power_CP1); hold on;
plot(freq, singlecued_power_CP1); hold on;
xline(4); hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
legend('Dual Uncued','Dual Cued', 'Single Uncued', 'Single Cued');

saveas(gcf, fullfile(results_path, 'CP1_PowervsFreq'),'png');

% PPC: P3
dualuncued_power_P3 = dualuncued_power(:, P3_loc);
dualcued_power_P3 = dualcued_power(:, P3_loc);
singleuncued_power_P3 = singleuncued_power(:, P3_loc);
singlecued_power_P3 = singlecued_power(:, P3_loc);

figure; title('P3');
plot(freq, dualuncued_power_P3); hold on;
plot(freq, dualcued_power_P3); hold on;
plot(freq, singleuncued_power_P3); hold on;
plot(freq, singlecued_power_P3); hold on;
xline(4); hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
legend('Dual Uncued','Dual Cued', 'Single Uncued', 'Single Cued');

saveas(gcf, fullfile(results_path, 'P3_PowervsFreq'),'png');

% PPC: P4
dualuncued_power_P4 = dualuncued_power(:, P4_loc);
dualcued_power_P4 = dualcued_power(:, P4_loc);
singleuncued_power_P4 = singleuncued_power(:, P4_loc);
singlecued_power_P4 = singlecued_power(:, P4_loc);

figure; title('P4');
plot(freq, dualuncued_power_P4); hold on;
plot(freq, dualcued_power_P4); hold on;
plot(freq, singleuncued_power_P4); hold on;
plot(freq, singlecued_power_P4); hold on;
xline(4); hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
legend('Dual Uncued','Dual Cued', 'Single Uncued', 'Single Cued');

saveas(gcf, fullfile(results_path, 'P4_PowervsFreq'),'png');

disp('This was the end of individual subjects.');
disp('These are the results for the average of all subjects.');
