%% Analysis of the EEG signal (PSD: separate channels)
clear; clc; close all;

%% Initialize data
addpath('C:\Users\catar\OneDrive - Universidade do Porto\Twente\Combined-EEG-fNIRS-system\Analysis');
addpath('C:\Users\catar\OneDrive - Universidade do Porto\Twente\Combined-EEG-fNIRS-system\Analysis\Dual Task');
laptop = 'laptopCatarina';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);
results_path = 'C:\Users\catar\OneDrive - Universidade do Porto\Twente\Data Analysis\eeg\psd';

eeglab;

subrec = ["02" "02";"64" "01";"28" "04"];

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
    
    % Calculate threshold to eliminate noisy epochs.
    th = calculateThreshold(EEG_divided,sub);
    
    %% Dual Uncued: PSD
    event_samp  = [EEG_DualUncued.event.latency];
    startTask = find(strcmp({EEG_DualUncued.event.type}, 's1709')==1);
    endTask = find(strcmp({EEG_DualUncued.event.type}, 's1717')==1);
    keypresses = find((strcmp({EEG_DualUncued.event.type}, 's1777') | strcmp({EEG_DualUncued.event.type}, 's1797'))==1);  
   
    % Get the PSD averaged over all trials
    [power, freq] = calculateAveragePowerAllTrials(EEG_DualUncued,...
        event_samp, startTask, endTask, keypresses, th);
    
    % Compensate for removed channels
    power = compensateRemovedChannelsPSD(power, EEG_DualUncued, list_channels, sub);
    
    if sub == "64"
       power(:,31)=[]; 
    end
    
    % Save the values onto a allSubjects variable
    dualuncued_power_allSubjects(:, :, subject) = power;
    
    %% Single Uncued: PSD
    event_samp  = [EEG_SingleUncued.event.latency];
    startTask = find(strcmp({EEG_SingleUncued.event.type}, 's1705')==1);
    endTask = find(strcmp({EEG_SingleUncued.event.type}, 's1713')==1);
    keypresses = find((strcmp({EEG_SingleUncued.event.type}, 's1777') | strcmp({EEG_SingleUncued.event.type}, 's1797'))==1);  
   
    % Get the PSD averaged over all trials
    [power, freq] = calculateAveragePowerAllTrials(EEG_SingleUncued,...
        event_samp, startTask, endTask, keypresses, th);
    
    % Compensate for removed channels
    power = compensateRemovedChannelsPSD(power, EEG_SingleUncued, list_channels, sub);
    
    if sub == "64"
       power(:,31)=[]; 
    end
    
    % Save the values onto a allSubjects variable
    singleuncued_power_allSubjects(:, :, subject) = power;
    
    %% Dual Cued: PSD
    event_samp  = [EEG_DualCued.event.latency];
    startTask = find(strcmp({EEG_DualCued.event.type}, 's1708')==1);
    endTask = find(strcmp({EEG_DualCued.event.type}, 's1701')==1);
    keypresses = find((strcmp({EEG_DualCued.event.type}, 's1777') | strcmp({EEG_DualCued.event.type}, 's1797'))==1);  
   
    % Get the PSD averaged over all trials
    [power, freq] = calculateAveragePowerAllTrials(EEG_DualCued,...
        event_samp, startTask, endTask, keypresses, th);
    
    % Compensate for removed channels
    power = compensateRemovedChannelsPSD(power, EEG_DualCued, list_channels, sub);
    
    if sub == "64"
       power(:,31)=[]; 
    end
    
    % Save the values onto a allSubjects variable
    dualcued_power_allSubjects(:, :, subject) = power;
    
    %% Single Cued: PSD
    event_samp  = [EEG_SingleCued.event.latency];
    startTask = find(strcmp({EEG_SingleCued.event.type}, 's1704')==1);
    endTask = find(strcmp({EEG_SingleCued.event.type}, 's1701')==1);
    keypresses = find((strcmp({EEG_SingleCued.event.type}, 's1777') | strcmp({EEG_SingleCued.event.type}, 's1797'))==1);  
   
    % Get the power PSD averaged over all trials
    [power, freq] = calculateAveragePowerAllTrials(EEG_SingleCued,...
        event_samp, startTask, endTask, keypresses, th);
    
    % Compensate for removed channels
    power = compensateRemovedChannelsPSD(power, EEG_SingleCued, list_channels, sub);
    
    if sub == "64"
       power(:,31)=[]; 
    end
    
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
    plot(freq, dualuncued_power_F7, '-b'); hold on;
    plot(freq, dualcued_power_F7, '-y'); hold on;
    plot(freq, singleuncued_power_F7, '-g'); hold on;
    plot(freq, singlecued_power_F7, '-r'); hold on;
    xline(4); hold on;
    xline(8); hold on;
    xline(13); hold on;
    xline(32); hold off;
    legend('Dual Uncued','Dual Cued', 'Single Uncued', 'Single Cued');
    
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, fullfile(results_path, ['sub-', char(sub)], 'F7'),'png');
    
    % Save the values onto a subject struct
    s.dualuncued_power_F7 = dualuncued_power_F7;
    s.singleuncued_power_F7 = singleuncued_power_F7;
    s.dualcued_power_F7 = dualcued_power_F7;
    s.singlecued_power_F7 = singlecued_power_F7;
    
    % PFC: F8
    dualuncued_power_F8 = dualuncued_power_allSubjects(:, F8_loc, subject);
    dualcued_power_F8 = dualcued_power_allSubjects(:, F8_loc, subject);
    singleuncued_power_F8 = singleuncued_power_allSubjects(:, F8_loc, subject);
    singlecued_power_F8 = singlecued_power_allSubjects(:, F8_loc, subject);
    
    figure; title('F8');
    plot(freq, dualuncued_power_F8, '-b'); hold on;
    plot(freq, dualcued_power_F8, '-y'); hold on;
    plot(freq, singleuncued_power_F8, '-g'); hold on;
    plot(freq, singlecued_power_F8, '-r'); hold on;
    xline(4); hold on;
    xline(8); hold on;
    xline(13); hold on;
    xline(32); hold off;
    legend('Dual Uncued','Dual Cued', 'Single Uncued', 'Single Cued');
    
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, fullfile(results_path, ['sub-', char(sub)], 'F8'),'png');
    
    % Save the values onto a subject struct
    s.dualuncued_power_F8 = dualuncued_power_F8;
    s.singleuncued_power_F8 = singleuncued_power_F8;
    s.dualcued_power_F8 = dualcued_power_F8;
    s.singlecued_power_F8 = singlecued_power_F8;
    
    % SMA: FC1
    dualuncued_power_FC1 = dualuncued_power_allSubjects(:, FC1_loc, subject);
    dualcued_power_FC1 = dualcued_power_allSubjects(:, FC1_loc, subject);
    singleuncued_power_FC1 = singleuncued_power_allSubjects(:, FC1_loc, subject);
    singlecued_power_FC1 = singlecued_power_allSubjects(:, FC1_loc, subject);
    
    figure; title('FC1');
    plot(freq, dualuncued_power_FC1, '-b'); hold on;
    plot(freq, dualcued_power_FC1, '-y'); hold on;
    plot(freq, singleuncued_power_FC1, '-g'); hold on;
    plot(freq, singlecued_power_FC1, '-r'); hold on;
    xline(4); hold on;
    xline(8); hold on;
    xline(13); hold on;
    xline(32); hold off;
    legend('Dual Uncued','Dual Cued', 'Single Uncued', 'Single Cued');
    
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, fullfile(results_path, ['sub-', char(sub)], 'FC1'),'png');
    
    % Save the values onto a subject struct
    s.dualuncued_power_FC1 = dualuncued_power_FC1;
    s.singleuncued_power_FC1 = singleuncued_power_FC1;
    s.dualcued_power_FC1 = dualcued_power_FC1;
    s.singlecued_power_FC1 = singlecued_power_FC1;
    
    % SMA: FC2
    dualuncued_power_FC2 = dualuncued_power_allSubjects(:, FC2_loc, subject);
    dualcued_power_FC2 = dualcued_power_allSubjects(:, FC2_loc, subject);
    singleuncued_power_FC2 = singleuncued_power_allSubjects(:, FC2_loc, subject);
    singlecued_power_FC2 = singlecued_power_allSubjects(:, FC2_loc, subject);
    
    figure; title('FC2');
    plot(freq, dualuncued_power_FC2, '-b'); hold on;
    plot(freq, dualcued_power_FC2, '-y'); hold on;
    plot(freq, singleuncued_power_FC2, '-g'); hold on;
    plot(freq, singlecued_power_FC2, '-r'); hold on;
    xline(4); hold on;
    xline(8); hold on;
    xline(13); hold on;
    xline(32); hold off;
    legend('Dual Uncued','Dual Cued', 'Single Uncued', 'Single Cued');
    
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, fullfile(results_path, ['sub-', char(sub)], 'FC2'),'png');
    
    % Save the values onto a subject struct
    s.dualuncued_power_FC2 = dualuncued_power_FC2;
    s.singleuncued_power_FC2 = singleuncued_power_FC2;
    s.dualcued_power_FC2 = dualcued_power_FC2;
    s.singlecued_power_FC2 = singlecued_power_FC2;
    
    % M1: C3
    dualuncued_power_C3 = dualuncued_power_allSubjects(:, C3_loc, subject);
    dualcued_power_C3 = dualcued_power_allSubjects(:, C3_loc, subject);
    singleuncued_power_C3 = singleuncued_power_allSubjects(:, C3_loc, subject);
    singlecued_power_C3 = singlecued_power_allSubjects(:, C3_loc, subject);
    
    figure; title('C3');
    plot(freq, dualuncued_power_C3, '-b'); hold on;
    plot(freq, dualcued_power_C3, '-y'); hold on;
    plot(freq, singleuncued_power_C3, '-g'); hold on;
    plot(freq, singlecued_power_C3, '-r'); hold on;
    xline(4); hold on;
    xline(8); hold on;
    xline(13); hold on;
    xline(32); hold off;
    legend('Dual Uncued','Dual Cued', 'Single Uncued', 'Single Cued');
    
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, fullfile(results_path, ['sub-', char(sub)], 'C3'),'png');
    
    % Save the values onto a subject struct
    s.dualuncued_power_C3 = dualuncued_power_C3;
    s.singleuncued_power_C3 = singleuncued_power_C3;
    s.dualcued_power_C3 = dualcued_power_C3;
    s.singlecued_power_C3 = singlecued_power_C3;
    
    % M1: CP1
    if CP1_loc ~= 0
    dualuncued_power_CP1 = dualuncued_power_allSubjects(:, CP1_loc, subject);
    dualcued_power_CP1 = dualcued_power_allSubjects(:, CP1_loc, subject);
    singleuncued_power_CP1 = singleuncued_power_allSubjects(:, CP1_loc, subject);
    singlecued_power_CP1 = singlecued_power_allSubjects(:, CP1_loc, subject);
    
    figure; title('CP1');
    plot(freq, dualuncued_power_CP1, '-b'); hold on;
    plot(freq, dualcued_power_CP1, '-y'); hold on;
    plot(freq, singleuncued_power_CP1, '-g'); hold on;
    plot(freq, singlecued_power_CP1, '-r'); hold on;
    xline(4); hold on;
    xline(8); hold on;
    xline(13); hold on;
    xline(32); hold off;
    legend('Dual Uncued','Dual Cued', 'Single Uncued', 'Single Cued');
    
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, fullfile(results_path, ['sub-', char(sub)], 'CP1'),'png');
    
    % Save the values onto a subject struct
    s.dualuncued_power_CP1 = dualuncued_power_CP1;
    s.singleuncued_power_CP1 = singleuncued_power_CP1;
    s.dualcued_power_CP1 = dualcued_power_CP1;
    s.singlecued_power_CP1 = singlecued_power_CP1;
    
    end
    
    % PPC: P3
    dualuncued_power_P3 = dualuncued_power_allSubjects(:, P3_loc, subject);
    dualcued_power_P3 = dualcued_power_allSubjects(:, P3_loc, subject);
    singleuncued_power_P3 = singleuncued_power_allSubjects(:, P3_loc, subject);
    singlecued_power_P3 = singlecued_power_allSubjects(:, P3_loc, subject);
    
    figure; title('P3');
    plot(freq, dualuncued_power_P3, '-b'); hold on;
    plot(freq, dualcued_power_P3, '-y'); hold on;
    plot(freq, singleuncued_power_P3, '-g'); hold on;
    plot(freq, singlecued_power_P3, '-r'); hold on;
    xline(4); hold on;
    xline(8); hold on;
    xline(13); hold on;
    xline(32); hold off;
    legend('Dual Uncued','Dual Cued', 'Single Uncued', 'Single Cued');
    
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, fullfile(results_path, ['sub-', char(sub)], 'P3'),'png');
    
    % Save the values onto a subject struct
    s.dualuncued_power_P3 = dualuncued_power_P3;
    s.singleuncued_power_P3 = singleuncued_power_P3;
    s.dualcued_power_P3 = dualcued_power_P3;
    s.singlecued_power_P3 = singlecued_power_P3;
    
    % PPC: P4
    dualuncued_power_P4 = dualuncued_power_allSubjects(:, P4_loc, subject);
    dualcued_power_P4 = dualcued_power_allSubjects(:, P4_loc, subject);
    singleuncued_power_P4 = singleuncued_power_allSubjects(:, P4_loc, subject);
    singlecued_power_P4 = singlecued_power_allSubjects(:, P4_loc, subject);
    
    figure; title('P4');
    plot(freq, dualuncued_power_P4, '-b'); hold on;
    plot(freq, dualcued_power_P4, '-y'); hold on;
    plot(freq, singleuncued_power_P4, '-g'); hold on;
    plot(freq, singlecued_power_P4, '-r'); hold on;
    xline(4); hold on;
    xline(8); hold on;
    xline(13); hold on;
    xline(32); hold off;
    legend('Dual Uncued','Dual Cued', 'Single Uncued', 'Single Cued');
    
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, fullfile(results_path, ['sub-', char(sub)], 'P4'),'png');
    
    % Save the values onto a subject struct
    s.dualuncued_power_P4 = dualuncued_power_P4;
    s.singleuncued_power_P4 = singleuncued_power_P4;
    s.dualcued_power_P4 = dualcued_power_P4;
    s.singlecued_power_P4 = singlecued_power_P4;
    
    % Add struct of current subject to all subjects struct.
    allsubs.(genvarname(strcat('sub', char(sub)))) = s;
    
    disp(['These are the results for subject ', char(sub), '.']);
    disp('Press any key to move onto the next subject.');
    %pause;
    close all;
    
end
disp('This was the end of individual subjects.');

%% Get the PSD averaged over all subjects
% Dual Cued
dualcued_power = mean(dualcued_power_allSubjects, 3, 'omitnan');
% Dual Uncued
dualuncued_power = mean(dualuncued_power_allSubjects, 3, 'omitnan');
% Single Cued
singlecued_power = mean(singlecued_power_allSubjects, 3, 'omitnan');
% Single Uncued
singleuncued_power = mean(singleuncued_power_allSubjects, 3, 'omitnan');

%% Get the standard deviation over all subjects
% Dual Cued
std_dualcued_power = std(dualcued_power_allSubjects, 0, 3, 'omitnan');
% Dual Uncued
std_dualuncued_power = std(dualuncued_power_allSubjects, 0, 3, 'omitnan');
% Single Cued
std_singlecued_power = std(singlecued_power_allSubjects, 0, 3, 'omitnan');
% Single Uncued
std_singleuncued_power = std(singleuncued_power_allSubjects, 0, 3, 'omitnan');

%% Plot the PSD for specific channels
% PFC: F7
dualuncued_power_F7 = dualuncued_power(:, F7_loc);
dualcued_power_F7 = dualcued_power(:, F7_loc);
singleuncued_power_F7 = singleuncued_power(:, F7_loc);
singlecued_power_F7 = singlecued_power(:, F7_loc);
std_dualuncued_power_F7 = std_dualuncued_power(:, F7_loc);
std_dualcued_power_F7 = std_dualcued_power(:, F7_loc);
std_singleuncued_power_F7 = std_singleuncued_power(:, F7_loc);
std_singlecued_power_F7 = std_singlecued_power(:, F7_loc);

figure; title('F7');
tiledlayout(1, 2);
% Dual Task
ax1 = nexttile;
plot(freq, dualuncued_power_F7, '-b', 'LineWidth', 1.5); 
hold on;
h1 = patch([freq; flipud(freq)],...
    [(dualuncued_power_F7 - std_dualuncued_power_F7);...
    flipud(dualuncued_power_F7 + std_dualuncued_power_F7)], 'b',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, dualcued_power_F7, '-y', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(dualcued_power_F7 - std_dualcued_power_F7);...
    flipud(dualcued_power_F7 + std_dualcued_power_F7)], 'y',...
    'LineStyle', 'none'); set(h2,'FaceAlpha',0.2);
hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
xlim([4 48]);
legend('Dual Uncued', '', 'Dual Cued', '');
xlabel('Frequency (Hz)');
ylabel('PSD');

% Single Task
ax2 = nexttile;
plot(freq, singleuncued_power_F7, '-g', 'LineWidth', 1.5);
hold on;
h1 = patch([freq; flipud(freq)],...
    [(singleuncued_power_F7 - std_singleuncued_power_F7);...
    flipud(singleuncued_power_F7 + std_singleuncued_power_F7)], 'g',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, singlecued_power_F7, '-r', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(singlecued_power_F7 - std_singlecued_power_F7);...
    flipud(singlecued_power_F7 + std_singlecued_power_F7)], 'r',...
    'LineStyle', 'none'); set(h2,'FaceAlpha',0.2);
hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
xlim([4 48]);
legend('Single Uncued', '', 'Single Cued', '');
xlabel('Frequency (Hz)');
ylabel('PSD');
linkaxes([ax1 ax2],'y')

set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, fullfile(results_path, 'F7_PowervsFreq'),'png');

% PFC: F8
dualuncued_power_F8 = dualuncued_power(:, F8_loc);
dualcued_power_F8 = dualcued_power(:, F8_loc);
singleuncued_power_F8 = singleuncued_power(:, F8_loc);
singlecued_power_F8 = singlecued_power(:, F8_loc);
std_dualuncued_power_F8 = std_dualuncued_power(:, F8_loc);
std_dualcued_power_F8 = std_dualcued_power(:, F8_loc);
std_singleuncued_power_F8 = std_singleuncued_power(:, F8_loc);
std_singlecued_power_F8 = std_singlecued_power(:, F8_loc);

figure; title('F8');
tiledlayout(1, 2);
% Dual Task
ax1 = nexttile;
plot(freq, dualuncued_power_F8, '-b', 'LineWidth', 1.5); 
hold on;
h1 = patch([freq; flipud(freq)],...
    [(dualuncued_power_F8 - std_dualuncued_power_F8);...
    flipud(dualuncued_power_F8 + std_dualuncued_power_F8)], 'b',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, dualcued_power_F8, '-y', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(dualcued_power_F8 - std_dualcued_power_F8);...
    flipud(dualcued_power_F8 + std_dualcued_power_F8)], 'y',...
    'LineStyle', 'none'); set(h2,'FaceAlpha',0.2);
hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
xlim([4 48]);
legend('Dual Uncued', '', 'Dual Cued', '');
xlabel('Frequency (Hz)');
ylabel('PSD');

% Single Task
ax2 = nexttile;
plot(freq, singleuncued_power_F8, '-g', 'LineWidth', 1.5);
hold on;
h1 = patch([freq; flipud(freq)],...
    [(singleuncued_power_F8 - std_singleuncued_power_F8);...
    flipud(singleuncued_power_F8 + std_singleuncued_power_F8)], 'g',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, singlecued_power_F8, '-r', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(singlecued_power_F8 - std_singlecued_power_F8);...
    flipud(singlecued_power_F8 + std_singlecued_power_F8)], 'r',...
    'LineStyle', 'none'); set(h2,'FaceAlpha',0.2);
hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
xlim([4 48]);
legend('Single Uncued', '', 'Single Cued', '');
xlabel('Frequency (Hz)');
ylabel('PSD');
linkaxes([ax1 ax2],'y')

set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, fullfile(results_path, 'F8_PowervsFreq'),'png');

% SMA: FC1
dualuncued_power_FC1 = dualuncued_power(:, FC1_loc);
dualcued_power_FC1 = dualcued_power(:, FC1_loc);
singleuncued_power_FC1 = singleuncued_power(:, FC1_loc);
singlecued_power_FC1 = singlecued_power(:, FC1_loc);
std_dualuncued_power_FC1 = std_dualuncued_power(:, FC1_loc);
std_dualcued_power_FC1 = std_dualcued_power(:, FC1_loc);
std_singleuncued_power_FC1 = std_singleuncued_power(:, FC1_loc);
std_singlecued_power_FC1 = std_singlecued_power(:, FC1_loc);

figure; title('FC1');
tiledlayout(1, 2);
% Dual Task
ax1 = nexttile;
plot(freq, dualuncued_power_FC1, '-b', 'LineWidth', 1.5); 
hold on;
h1 = patch([freq; flipud(freq)],...
    [(dualuncued_power_FC1 - std_dualuncued_power_FC1);...
    flipud(dualuncued_power_FC1 + std_dualuncued_power_FC1)], 'b',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, dualcued_power_FC1, '-y', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(dualcued_power_FC1 - std_dualcued_power_FC1);...
    flipud(dualcued_power_FC1 + std_dualcued_power_FC1)], 'y',...
    'LineStyle', 'none'); set(h2,'FaceAlpha',0.2);
hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
xlim([4 48]);
legend('Dual Uncued', '', 'Dual Cued', '');
xlabel('Frequency (Hz)');
ylabel('PSD');

% Single Task
ax2 = nexttile;
plot(freq, singleuncued_power_FC1, '-g', 'LineWidth', 1.5);
hold on;
h1 = patch([freq; flipud(freq)],...
    [(singleuncued_power_FC1 - std_singleuncued_power_FC1);...
    flipud(singleuncued_power_FC1 + std_singleuncued_power_FC1)], 'g',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, singlecued_power_FC1, '-r', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(singlecued_power_FC1 - std_singlecued_power_FC1);...
    flipud(singlecued_power_FC1 + std_singlecued_power_FC1)], 'r',...
    'LineStyle', 'none'); set(h2,'FaceAlpha',0.2);
hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
xlim([4 48]);
legend('Single Uncued', '', 'Single Cued', '');
xlabel('Frequency (Hz)');
ylabel('PSD');
linkaxes([ax1 ax2],'y')

set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, fullfile(results_path, 'FC1_PowervsFreq'),'png');

% SMA: FC2
dualuncued_power_FC2 = dualuncued_power(:, FC2_loc);
dualcued_power_FC2 = dualcued_power(:, FC2_loc);
singleuncued_power_FC2 = singleuncued_power(:, FC2_loc);
singlecued_power_FC2 = singlecued_power(:, FC2_loc);
std_dualuncued_power_FC2 = std_dualuncued_power(:, FC2_loc);
std_dualcued_power_FC2 = std_dualcued_power(:, FC2_loc);
std_singleuncued_power_FC2 = std_singleuncued_power(:, FC2_loc);
std_singlecued_power_FC2 = std_singlecued_power(:, FC2_loc);

figure; title('FC2');
tiledlayout(1, 2);
% Dual Task
ax1 = nexttile;
plot(freq, dualuncued_power_FC2, '-b', 'LineWidth', 1.5); 
hold on;
h1 = patch([freq; flipud(freq)],...
    [(dualuncued_power_FC2 - std_dualuncued_power_FC2);...
    flipud(dualuncued_power_FC2 + std_dualuncued_power_FC2)], 'b',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, dualcued_power_FC2, '-y', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(dualcued_power_FC2 - std_dualcued_power_FC2);...
    flipud(dualcued_power_FC2 + std_dualcued_power_FC2)], 'y',...
    'LineStyle', 'none'); set(h2,'FaceAlpha',0.2);
hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
xlim([4 48]);
legend('Dual Uncued', '', 'Dual Cued', '');
xlabel('Frequency (Hz)');
ylabel('PSD');

% Single Task
ax2 = nexttile;
plot(freq, singleuncued_power_FC2, '-g', 'LineWidth', 1.5);
hold on;
h1 = patch([freq; flipud(freq)],...
    [(singleuncued_power_FC2 - std_singleuncued_power_FC2);...
    flipud(singleuncued_power_FC2 + std_singleuncued_power_FC2)], 'g',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, singlecued_power_FC2, '-r', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(singlecued_power_FC2 - std_singlecued_power_FC2);...
    flipud(singlecued_power_FC2 + std_singlecued_power_FC2)], 'r',...
    'LineStyle', 'none'); set(h2,'FaceAlpha',0.2);
hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
xlim([4 48]);
legend('Single Uncued', '', 'Single Cued', '');
xlabel('Frequency (Hz)');
ylabel('PSD');
linkaxes([ax1 ax2],'y')

set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, fullfile(results_path, 'FC2_PowervsFreq'),'png');

% M1: C3
dualuncued_power_C3 = dualuncued_power(:, C3_loc);
dualcued_power_C3 = dualcued_power(:, C3_loc);
singleuncued_power_C3 = singleuncued_power(:, C3_loc);
singlecued_power_C3 = singlecued_power(:, C3_loc);
std_dualuncued_power_C3 = std_dualuncued_power(:, C3_loc);
std_dualcued_power_C3 = std_dualcued_power(:, C3_loc);
std_singleuncued_power_C3 = std_singleuncued_power(:, C3_loc);
std_singlecued_power_C3 = std_singlecued_power(:, C3_loc);

figure; title('C3');
tiledlayout(1, 2);
% Dual Task
ax1 = nexttile;
plot(freq, dualuncued_power_C3, '-b', 'LineWidth', 1.5); 
hold on;
h1 = patch([freq; flipud(freq)],...
    [(dualuncued_power_C3 - std_dualuncued_power_C3);...
    flipud(dualuncued_power_C3 + std_dualuncued_power_C3)], 'b',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, dualcued_power_C3, '-y', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(dualcued_power_C3 - std_dualcued_power_C3);...
    flipud(dualcued_power_C3 + std_dualcued_power_C3)], 'y',...
    'LineStyle', 'none'); set(h2,'FaceAlpha',0.2);
hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
xlim([4 48]);
legend('Dual Uncued', '', 'Dual Cued', '');
xlabel('Frequency (Hz)');
ylabel('PSD');

% Single Task
ax2 = nexttile;
plot(freq, singleuncued_power_C3, '-g', 'LineWidth', 1.5);
hold on;
h1 = patch([freq; flipud(freq)],...
    [(singleuncued_power_C3 - std_singleuncued_power_C3);...
    flipud(singleuncued_power_C3 + std_singleuncued_power_C3)], 'g',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, singlecued_power_C3, '-r', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(singlecued_power_C3 - std_singlecued_power_C3);...
    flipud(singlecued_power_C3 + std_singlecued_power_C3)], 'r',...
    'LineStyle', 'none'); set(h2,'FaceAlpha',0.2);
hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
xlim([4 48]);
legend('Single Uncued', '', 'Single Cued', '');
xlabel('Frequency (Hz)');
ylabel('PSD');
linkaxes([ax1 ax2],'y')

set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, fullfile(results_path, 'C3_PowervsFreq'),'png');

% M1: CP1
dualuncued_power_CP1 = dualuncued_power(:, CP1_loc);
dualcued_power_CP1 = dualcued_power(:, CP1_loc);
singleuncued_power_CP1 = singleuncued_power(:, CP1_loc);
singlecued_power_CP1 = singlecued_power(:, CP1_loc);
std_dualuncued_power_CP1 = std_dualuncued_power(:, CP1_loc);
std_dualcued_power_CP1 = std_dualcued_power(:, CP1_loc);
std_singleuncued_power_CP1 = std_singleuncued_power(:, CP1_loc);
std_singlecued_power_CP1 = std_singlecued_power(:, CP1_loc);

figure; title('CP1');
tiledlayout(1, 2);
% Dual Task
ax1 = nexttile;
plot(freq, dualuncued_power_CP1, '-b', 'LineWidth', 1.5); 
hold on;
h1 = patch([freq; flipud(freq)],...
    [(dualuncued_power_CP1 - std_dualuncued_power_CP1);...
    flipud(dualuncued_power_CP1 + std_dualuncued_power_CP1)], 'b',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, dualcued_power_CP1, '-y', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(dualcued_power_CP1 - std_dualcued_power_CP1);...
    flipud(dualcued_power_CP1 + std_dualcued_power_CP1)], 'y',...
    'LineStyle', 'none'); set(h2,'FaceAlpha',0.2);
hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
xlim([4 48]);
legend('Dual Uncued', '', 'Dual Cued', '');
xlabel('Frequency (Hz)');
ylabel('PSD');

% Single Task
ax2 = nexttile;
plot(freq, singleuncued_power_CP1, '-g', 'LineWidth', 1.5);
hold on;
h1 = patch([freq; flipud(freq)],...
    [(singleuncued_power_CP1 - std_singleuncued_power_CP1);...
    flipud(singleuncued_power_CP1 + std_singleuncued_power_CP1)], 'g',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, singlecued_power_CP1, '-r', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(singlecued_power_CP1 - std_singlecued_power_CP1);...
    flipud(singlecued_power_CP1 + std_singlecued_power_CP1)], 'r',...
    'LineStyle', 'none'); set(h2,'FaceAlpha',0.2);
hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
xlim([4 48]);
legend('Single Uncued', '', 'Single Cued', '');
xlabel('Frequency (Hz)');
ylabel('PSD');
linkaxes([ax1 ax2],'y')

set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, fullfile(results_path, 'CP1_PowervsFreq'),'png');

% PPC: P3
dualuncued_power_P3 = dualuncued_power(:, P3_loc);
dualcued_power_P3 = dualcued_power(:, P3_loc);
singleuncued_power_P3 = singleuncued_power(:, P3_loc);
singlecued_power_P3 = singlecued_power(:, P3_loc);
std_dualuncued_power_P3 = std_dualuncued_power(:, P3_loc);
std_dualcued_power_P3 = std_dualcued_power(:, P3_loc);
std_singleuncued_power_P3 = std_singleuncued_power(:, P3_loc);
std_singlecued_power_P3 = std_singlecued_power(:, P3_loc);

figure; title('P3');
tiledlayout(1, 2);
% Dual Task
ax1 = nexttile;
plot(freq, dualuncued_power_P3, '-b', 'LineWidth', 1.5); 
hold on;
h1 = patch([freq; flipud(freq)],...
    [(dualuncued_power_P3 - std_dualuncued_power_P3);...
    flipud(dualuncued_power_P3 + std_dualuncued_power_P3)], 'b',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, dualcued_power_P3, '-y', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(dualcued_power_P3 - std_dualcued_power_P3);...
    flipud(dualcued_power_P3 + std_dualcued_power_P3)], 'y',...
    'LineStyle', 'none'); set(h2,'FaceAlpha',0.2);
hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
xlim([4 48]);
legend('Dual Uncued', '', 'Dual Cued', '');
xlabel('Frequency (Hz)');
ylabel('PSD');

% Single Task
ax2 = nexttile;
plot(freq, singleuncued_power_P3, '-g', 'LineWidth', 1.5);
hold on;
h1 = patch([freq; flipud(freq)],...
    [(singleuncued_power_P3 - std_singleuncued_power_P3);...
    flipud(singleuncued_power_P3 + std_singleuncued_power_P3)], 'g',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, singlecued_power_P3, '-r', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(singlecued_power_P3 - std_singlecued_power_P3);...
    flipud(singlecued_power_P3 + std_singlecued_power_P3)], 'r',...
    'LineStyle', 'none'); set(h2,'FaceAlpha',0.2);
hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
xlim([4 48]);
legend('Single Uncued', '', 'Single Cued', '');
xlabel('Frequency (Hz)');
ylabel('PSD');
linkaxes([ax1 ax2],'y')

set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, fullfile(results_path, 'P3_PowervsFreq'),'png');

% PPC: P4
dualuncued_power_P4 = dualuncued_power(:, P4_loc);
dualcued_power_P4 = dualcued_power(:, P4_loc);
singleuncued_power_P4 = singleuncued_power(:, P4_loc);
singlecued_power_P4 = singlecued_power(:, P4_loc);
std_dualuncued_power_P4 = std_dualuncued_power(:, P4_loc);
std_dualcued_power_P4 = std_dualcued_power(:, P4_loc);
std_singleuncued_power_P4 = std_singleuncued_power(:, P4_loc);
std_singlecued_power_P4 = std_singlecued_power(:, P4_loc);

figure; title('P4');
tiledlayout(1, 2);
% Dual Task
ax1 = nexttile;
plot(freq, dualuncued_power_P4, '-b', 'LineWidth', 1.5); 
hold on;
h1 = patch([freq; flipud(freq)],...
    [(dualuncued_power_P4 - std_dualuncued_power_P4);...
    flipud(dualuncued_power_P4 + std_dualuncued_power_P4)], 'b',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, dualcued_power_P4, '-y', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(dualcued_power_P4 - std_dualcued_power_P4);...
    flipud(dualcued_power_P4 + std_dualcued_power_P4)], 'y',...
    'LineStyle', 'none'); set(h2,'FaceAlpha',0.2);
hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
xlim([4 48]);
legend('Dual Uncued', '', 'Dual Cued', '');
xlabel('Frequency (Hz)');
ylabel('PSD');

% Single Task
ax2 = nexttile;
plot(freq, singleuncued_power_P4, '-g', 'LineWidth', 1.5);
hold on;
h1 = patch([freq; flipud(freq)],...
    [(singleuncued_power_P4 - std_singleuncued_power_P4);...
    flipud(singleuncued_power_P4 + std_singleuncued_power_P4)], 'g',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, singlecued_power_P4, '-r', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(singlecued_power_P4 - std_singlecued_power_P4);...
    flipud(singlecued_power_P4 + std_singlecued_power_P4)], 'r',...
    'LineStyle', 'none'); set(h2,'FaceAlpha',0.2);
hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
xlim([4 48]);
legend('Single Uncued', '', 'Single Cued', '');
xlabel('Frequency (Hz)');
ylabel('PSD');
linkaxes([ax1 ax2],'y')

set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, fullfile(results_path, 'P4_PowervsFreq'),'png');

%% Save the values onto all subs struct
avg.dualuncued_power_F7 = dualuncued_power_F7;
avg.singleuncued_power_F7 = singleuncued_power_F7;
avg.dualcued_power_F7 = dualcued_power_F7;
avg.singlecued_power_F7 = singlecued_power_F7;
avg.dualuncued_power_F8 = dualuncued_power_F8;
avg.singleuncued_power_F8 = singleuncued_power_F8;
avg.dualcued_power_F8 = dualcued_power_F8;
avg.singlecued_power_F8 = singlecued_power_F8;
avg.dualuncued_power_FC1 = dualuncued_power_FC1;
avg.singleuncued_power_FC1 = singleuncued_power_FC1;
avg.dualcued_power_FC1 = dualcued_power_FC1;
avg.singlecued_power_FC1 = singlecued_power_FC1;
avg.dualuncued_power_FC2 = dualuncued_power_FC2;
avg.singleuncued_power_FC2 = singleuncued_power_FC2;
avg.dualcued_power_FC2 = dualcued_power_FC2;
avg.singlecued_power_FC2 = singlecued_power_FC2;
avg.dualuncued_power_C3 = dualuncued_power_C3;
avg.singleuncued_power_C3 = singleuncued_power_C3;
avg.dualcued_power_C3 = dualcued_power_C3;
avg.singlecued_power_C3 = singlecued_power_C3;
avg.dualuncued_power_CP1 = dualuncued_power_CP1;
avg.singleuncued_power_CP1 = singleuncued_power_CP1;
avg.dualcued_power_CP1 = dualcued_power_CP1;
avg.singlecued_power_CP1 = singlecued_power_CP1;
avg.dualuncued_power_P3 = dualuncued_power_P3;
avg.singleuncued_power_P3 = singleuncued_power_P3;
avg.dualcued_power_P3 = dualcued_power_P3;
avg.singlecued_power_P3 = singlecued_power_P3;
avg.dualuncued_power_P4 = dualuncued_power_P4;
avg.singleuncued_power_P4 = singleuncued_power_P4;
avg.dualcued_power_P4 = dualcued_power_P4;
avg.singlecued_power_P4 = singlecued_power_P4;
allsubs.avg = avg;

std.dualuncued_power_F7 = std_dualuncued_power_F7;
std.singleuncued_power_F7 = std_singleuncued_power_F7;
std.dualcued_power_F7 = std_dualcued_power_F7;
std.singlecued_power_F7 = std_singlecued_power_F7;
std.dualuncued_power_F8 = std_dualuncued_power_F8;
std.singleuncued_power_F8 = std_singleuncued_power_F8;
std.dualcued_power_F8 = std_dualcued_power_F8;
std.singlecued_power_F8 = std_singlecued_power_F8;
std.dualuncued_power_FC1 = std_dualuncued_power_FC1;
std.singleuncued_power_FC1 = std_singleuncued_power_FC1;
std.dualcued_power_FC1 = std_dualcued_power_FC1;
std.singlecued_power_FC1 = std_singlecued_power_FC1;
std.dualuncued_power_FC2 = std_dualuncued_power_FC2;
std.singleuncued_power_FC2 = std_singleuncued_power_FC2;
std.dualcued_power_FC2 = std_dualcued_power_FC2;
std.singlecued_power_FC2 = std_singlecued_power_FC2;
std.dualuncued_power_C3 = std_dualuncued_power_C3;
std.singleuncued_power_C3 = std_singleuncued_power_C3;
std.dualcued_power_C3 = std_dualcued_power_C3;
std.singlecued_power_C3 = std_singlecued_power_C3;
std.dualuncued_power_CP1 = std_dualuncued_power_CP1;
std.singleuncued_power_CP1 = std_singleuncued_power_CP1;
std.dualcued_power_CP1 = std_dualcued_power_CP1;
std.singlecued_power_CP1 = std_singlecued_power_CP1;
std.dualuncued_power_P3 = std_dualuncued_power_P3;
std.singleuncued_power_P3 = std_singleuncued_power_P3;
std.dualcued_power_P3 = std_dualcued_power_P3;
std.singlecued_power_P3 = std_singlecued_power_P3;
std.dualuncued_power_P4 = std_dualuncued_power_P4;
std.singleuncued_power_P4 = std_singleuncued_power_P4;
std.dualcued_power_P4 = std_dualcued_power_P4;
std.singlecued_power_P4 = std_singlecued_power_P4;
allsubs.std = std;

% Save the struct from all subs
save(strcat(results_path, '\psdchannel_allsubs.mat'), 'allsubs')

disp('These are the results for the average of all subjects.');
