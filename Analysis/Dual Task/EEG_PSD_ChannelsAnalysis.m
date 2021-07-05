%% Analysis of the EEG signal (PSD: separate channels)
clear; clc; close all;

%% Initialize data
addpath('C:\Users\catar\OneDrive - Universidade do Porto\Twente\Combined-EEG-fNIRS-system\Analysis');
addpath('C:\Users\catar\OneDrive - Universidade do Porto\Twente\Combined-EEG-fNIRS-system\Analysis\Dual Task');
laptop = 'laptopCatarina';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);
results_path = 'C:\Users\catar\OneDrive - Universidade do Porto\Twente\Data Analysis\eeg\psd';

eeglab;

subrec = ["76" "01";"64" "01";"02" "02";"28" "04"];

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
   
    if sub=="64"
        startTask = [37 50];
        endTask = [47 62];
    end
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
   
     if sub=="64"
        startTask = [164];
        endTask = [177];
     end
     
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
    % DLPFC
    dualuncued_power_F7 = dualuncued_power_allSubjects(:, F7_loc, subject);
    dualcued_power_F7 = dualcued_power_allSubjects(:, F7_loc, subject);
    singleuncued_power_F7 = singleuncued_power_allSubjects(:, F7_loc, subject);
    singlecued_power_F7 = singlecued_power_allSubjects(:, F7_loc, subject);
    
    dualuncued_power_F8 = dualuncued_power_allSubjects(:, F8_loc, subject);
    dualcued_power_F8 = dualcued_power_allSubjects(:, F8_loc, subject);
    singleuncued_power_F8 = singleuncued_power_allSubjects(:, F8_loc, subject);
    singlecued_power_F8 = singlecued_power_allSubjects(:, F8_loc, subject);
    
    dualuncued_power_DLPFC = mean([dualuncued_power_F7 dualuncued_power_F8], 2);
    dualcued_power_DLPFC = mean([dualcued_power_F7 dualcued_power_F8], 2);
    singleuncued_power_DLPFC = mean([singleuncued_power_F7 singleuncued_power_F8], 2);
    singlecued_power_DLPFC = mean([singlecued_power_F7 singlecued_power_F8], 2);
    
    figure; title('DLPFC');
    plot(freq, dualuncued_power_DLPFC, '-b'); hold on;
    plot(freq, dualcued_power_DLPFC, '-m'); hold on;
    plot(freq, singleuncued_power_DLPFC, '-g'); hold on;
    plot(freq, singlecued_power_DLPFC, '-r'); hold on;
    xline(4); hold on;
    xline(8); hold on;
    xline(13); hold on;
    xline(32); hold off;
    legend('Dual Uncued','Dual Cued', 'Single Uncued', 'Single Cued');
    
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, fullfile(results_path, ['sub-', char(sub)], 'DLPFC'),'png');
    
    % Save the values onto a subject struct.
    s.dualuncued_power_DLPFC = dualuncued_power_DLPFC;
    s.dualcued_power_DLPFC = dualcued_power_DLPFC;
    s.singleuncued_power_DLPFC = singleuncued_power_DLPFC;
    s.singlecued_power_DLPFC= singlecued_power_DLPFC;
    
    % SMA
    dualuncued_power_FC1 = dualuncued_power_allSubjects(:, FC1_loc, subject);
    dualcued_power_FC1 = dualcued_power_allSubjects(:, FC1_loc, subject);
    singleuncued_power_FC1 = singleuncued_power_allSubjects(:, FC1_loc, subject);
    singlecued_power_FC1 = singlecued_power_allSubjects(:, FC1_loc, subject);
    
    dualuncued_power_FC2 = dualuncued_power_allSubjects(:, FC2_loc, subject);
    dualcued_power_FC2 = dualcued_power_allSubjects(:, FC2_loc, subject);
    singleuncued_power_FC2 = singleuncued_power_allSubjects(:, FC2_loc, subject);
    singlecued_power_FC2 = singlecued_power_allSubjects(:, FC2_loc, subject);
    
    dualuncued_power_SMA = mean([dualuncued_power_FC1 dualuncued_power_FC2], 2);
    dualcued_power_SMA = mean([dualcued_power_FC1 dualcued_power_FC2], 2);
    singleuncued_power_SMA = mean([singleuncued_power_FC1 singleuncued_power_FC2], 2);
    singlecued_power_SMA = mean([singlecued_power_FC1 singlecued_power_FC2], 2);
    
    figure; title('SMA');
    plot(freq, dualuncued_power_SMA, '-b'); hold on;
    plot(freq, dualcued_power_SMA, '-m'); hold on;
    plot(freq, singleuncued_power_SMA, '-g'); hold on;
    plot(freq, singlecued_power_SMA, '-r'); hold on;
    xline(4); hold on;
    xline(8); hold on;
    xline(13); hold on;
    xline(32); hold off;
    legend('Dual Uncued','Dual Cued', 'Single Uncued', 'Single Cued');
    
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, fullfile(results_path, ['sub-', char(sub)], 'SMA'),'png');
    
    % Save the values onto a subject struct.
    s.dualuncued_power_SMA = dualuncued_power_SMA;
    s.dualcued_power_SMA = dualcued_power_SMA;
    s.singleuncued_power_SMA = singleuncued_power_SMA;
    s.singlecued_power_SMA= singlecued_power_SMA;
    
    % M1
    dualuncued_power_C3 = dualuncued_power_allSubjects(:, C3_loc, subject);
    dualcued_power_C3 = dualcued_power_allSubjects(:, C3_loc, subject);
    singleuncued_power_C3 = singleuncued_power_allSubjects(:, C3_loc, subject);
    singlecued_power_C3 = singlecued_power_allSubjects(:, C3_loc, subject);
    
    if CP1_loc ~= 0
    dualuncued_power_CP1 = dualuncued_power_allSubjects(:, CP1_loc, subject);
    dualcued_power_CP1 = dualcued_power_allSubjects(:, CP1_loc, subject);
    singleuncued_power_CP1 = singleuncued_power_allSubjects(:, CP1_loc, subject);
    singlecued_power_CP1 = singlecued_power_allSubjects(:, CP1_loc, subject);
    
    dualuncued_power_M1 = mean([dualuncued_power_C3 dualuncued_power_CP1], 2);
    dualcued_power_M1 = mean([dualcued_power_C3 dualcued_power_CP1], 2);
    singleuncued_power_M1 = mean([singleuncued_power_C3 singleuncued_power_CP1], 2);
    singlecued_power_M1 = mean([singlecued_power_C3 singlecued_power_CP1], 2);
    
    figure; title('M1');
    plot(freq, dualuncued_power_M1, '-b'); hold on;
    plot(freq, dualcued_power_M1, '-m'); hold on;
    plot(freq, singleuncued_power_M1, '-g'); hold on;
    plot(freq, singlecued_power_M1, '-r'); hold on;
    xline(4); hold on;
    xline(8); hold on;
    xline(13); hold on;
    xline(32); hold off;
    legend('Dual Uncued','Dual Cued', 'Single Uncued', 'Single Cued');
    
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, fullfile(results_path, ['sub-', char(sub)], 'M1'),'png');
    
    % Save the values onto a subject struct.
    s.dualuncued_power_M1 = dualuncued_power_M1;
    s.dualcued_power_M1 = dualcued_power_M1;
    s.singleuncued_power_M1 = singleuncued_power_M1;
    s.singlecued_power_M1 = singlecued_power_M1;
    
    else
        
    dualuncued_power_M1 = dualuncued_power_C3;
    dualcued_power_M1 = dualcued_power_C3;
    singleuncued_power_M1 = singleuncued_power_C3;
    singlecued_power_M1 = singlecued_power_C3;
    
    figure; title('M1');
    plot(freq, dualuncued_power_M1, '-b'); hold on;
    plot(freq, dualcued_power_M1, '-m'); hold on;
    plot(freq, singleuncued_power_M1, '-g'); hold on;
    plot(freq, singlecued_power_M1, '-r'); hold on;
    xline(4); hold on;
    xline(8); hold on;
    xline(13); hold on;
    xline(32); hold off;
    legend('Dual Uncued','Dual Cued', 'Single Uncued', 'Single Cued');
    
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, fullfile(results_path, ['sub-', char(sub)], 'M1'),'png');
    
    % Save the values onto a subject struct.
    s.dualuncued_power_M1 = dualuncued_power_M1;
    s.dualcued_power_M1 = dualcued_power_M1;
    s.singleuncued_power_M1 = singleuncued_power_M1;
    s.singlecued_power_M1 = singlecued_power_M1;
        
    end
    
    % PPC
    dualuncued_power_P3 = dualuncued_power_allSubjects(:, P3_loc, subject);
    dualcued_power_P3 = dualcued_power_allSubjects(:, P3_loc, subject);
    singleuncued_power_P3 = singleuncued_power_allSubjects(:, P3_loc, subject);
    singlecued_power_P3 = singlecued_power_allSubjects(:, P3_loc, subject);
    
    dualuncued_power_P4 = dualuncued_power_allSubjects(:, P4_loc, subject);
    dualcued_power_P4 = dualcued_power_allSubjects(:, P4_loc, subject);
    singleuncued_power_P4 = singleuncued_power_allSubjects(:, P4_loc, subject);
    singlecued_power_P4 = singlecued_power_allSubjects(:, P4_loc, subject);
    
    dualuncued_power_PPC = mean([dualuncued_power_P3 dualuncued_power_P4], 2);
    dualcued_power_PPC = mean([dualcued_power_P3 dualcued_power_P4], 2);
    singleuncued_power_PPC = mean([singleuncued_power_P3 singleuncued_power_P4], 2);
    singlecued_power_PPC = mean([singlecued_power_P3 singlecued_power_P4], 2);
    
    figure; title('PPC');
    plot(freq, dualuncued_power_PPC, '-b'); hold on;
    plot(freq, dualcued_power_PPC, '-m'); hold on;
    plot(freq, singleuncued_power_PPC, '-g'); hold on;
    plot(freq, singlecued_power_PPC, '-r'); hold on;
    xline(4); hold on;
    xline(8); hold on;
    xline(13); hold on;
    xline(32); hold off;
    legend('Dual Uncued','Dual Cued', 'Single Uncued', 'Single Cued');
    
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, fullfile(results_path, ['sub-', char(sub)], 'PPC'),'png');
    
    % Save the values onto a subject struct.
    s.dualuncued_power_PPC = dualuncued_power_PPC;
    s.dualcued_power_PPC = dualcued_power_PPC;
    s.singleuncued_power_PPC = singleuncued_power_PPC;
    s.singlecued_power_PPC = singlecued_power_PPC;
    
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

%% Plot the PSD for specific regions of interest
%% DLPFC
dualuncued_power_F7 = dualuncued_power(:, F7_loc);
dualcued_power_F7 = dualcued_power(:, F7_loc);
singleuncued_power_F7 = singleuncued_power(:, F7_loc);
singlecued_power_F7 = singlecued_power(:, F7_loc);
std_dualuncued_power_F7 = std_dualuncued_power(:, F7_loc);
std_dualcued_power_F7 = std_dualcued_power(:, F7_loc);
std_singleuncued_power_F7 = std_singleuncued_power(:, F7_loc);
std_singlecued_power_F7 = std_singlecued_power(:, F7_loc);

dualuncued_power_F8 = dualuncued_power(:, F8_loc);
dualcued_power_F8 = dualcued_power(:, F8_loc);
singleuncued_power_F8 = singleuncued_power(:, F8_loc);
singlecued_power_F8 = singlecued_power(:, F8_loc);
std_dualuncued_power_F8 = std_dualuncued_power(:, F8_loc);
std_dualcued_power_F8 = std_dualcued_power(:, F8_loc);
std_singleuncued_power_F8 = std_singleuncued_power(:, F8_loc);
std_singlecued_power_F8 = std_singlecued_power(:, F8_loc);

dualuncued_power_DLPFC = mean([dualuncued_power_F7 dualuncued_power_F8], 2);
dualcued_power_DLPFC = mean([dualcued_power_F7 dualcued_power_F8], 2);
singleuncued_power_DLPFC = mean([singleuncued_power_F7 singleuncued_power_F8], 2);
singlecued_power_DLPFC = mean([singlecued_power_F7 singlecued_power_F8], 2);
std_dualuncued_power_DLPFC = mean([std_dualuncued_power_F7 std_dualuncued_power_F8], 2);
std_dualcued_power_DLPFC = mean([std_dualcued_power_F7 std_dualcued_power_F8], 2);
std_singleuncued_power_DLPFC = mean([std_singleuncued_power_F7 std_singleuncued_power_F8], 2);
std_singlecued_power_DLPFC = mean([std_singlecued_power_F7 std_singlecued_power_F8], 2);

figure; title('DLPFC');
tiledlayout(1, 2);
% Dual Task
ax1 = nexttile;
plot(freq, dualuncued_power_DLPFC, '-b', 'LineWidth', 1.5); 
hold on;
h1 = patch([freq; flipud(freq)],...
    [(dualuncued_power_DLPFC - std_dualuncued_power_DLPFC);...
    flipud(dualuncued_power_DLPFC + std_dualuncued_power_DLPFC)], 'b',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, dualcued_power_DLPFC, '-m', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(dualcued_power_DLPFC - std_dualcued_power_DLPFC);...
    flipud(dualcued_power_DLPFC + std_dualcued_power_DLPFC)], 'm',...
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
plot(freq, singleuncued_power_DLPFC, '-g', 'LineWidth', 1.5);
hold on;
h1 = patch([freq; flipud(freq)],...
    [(singleuncued_power_DLPFC - std_singleuncued_power_DLPFC);...
    flipud(singleuncued_power_DLPFC + std_singleuncued_power_DLPFC)], 'g',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, singlecued_power_DLPFC, '-r', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(singlecued_power_DLPFC - std_singlecued_power_DLPFC);...
    flipud(singlecued_power_DLPFC + std_singlecued_power_DLPFC)], 'r',...
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
saveas(gcf, fullfile(results_path, 'DLPFC_PowervsFreq'),'png');

%% SMA
dualuncued_power_FC1 = dualuncued_power(:, FC1_loc);
dualcued_power_FC1 = dualcued_power(:, FC1_loc);
singleuncued_power_FC1 = singleuncued_power(:, FC1_loc);
singlecued_power_FC1 = singlecued_power(:, FC1_loc);
std_dualuncued_power_FC1 = std_dualuncued_power(:, FC1_loc);
std_dualcued_power_FC1 = std_dualcued_power(:, FC1_loc);
std_singleuncued_power_FC1 = std_singleuncued_power(:, FC1_loc);
std_singlecued_power_FC1 = std_singlecued_power(:, FC1_loc);

dualuncued_power_FC2 = dualuncued_power(:, FC2_loc);
dualcued_power_FC2 = dualcued_power(:, FC2_loc);
singleuncued_power_FC2 = singleuncued_power(:, FC2_loc);
singlecued_power_FC2 = singlecued_power(:, FC2_loc);
std_dualuncued_power_FC2 = std_dualuncued_power(:, FC2_loc);
std_dualcued_power_FC2 = std_dualcued_power(:, FC2_loc);
std_singleuncued_power_FC2 = std_singleuncued_power(:, FC2_loc);
std_singlecued_power_FC2 = std_singlecued_power(:, FC2_loc);

dualuncued_power_SMA = mean([dualuncued_power_FC1 dualuncued_power_FC2], 2);
dualcued_power_SMA = mean([dualcued_power_FC1 dualcued_power_FC2], 2);
singleuncued_power_SMA = mean([singleuncued_power_FC1 singleuncued_power_FC2], 2);
singlecued_power_SMA = mean([singlecued_power_FC1 singlecued_power_FC2], 2);
std_dualuncued_power_SMA = mean([std_dualuncued_power_FC1 std_dualuncued_power_FC2], 2);
std_dualcued_power_SMA = mean([std_dualcued_power_FC1 std_dualcued_power_FC2], 2);
std_singleuncued_power_SMA = mean([std_singleuncued_power_FC1 std_singleuncued_power_FC2], 2);
std_singlecued_power_SMA = mean([std_singlecued_power_FC1 std_singlecued_power_FC2], 2);

figure; title('SMA');
tiledlayout(1, 2);
% Dual Task
ax1 = nexttile;
plot(freq, dualuncued_power_SMA, '-b', 'LineWidth', 1.5); 
hold on;
h1 = patch([freq; flipud(freq)],...
    [(dualuncued_power_SMA - std_dualuncued_power_SMA);...
    flipud(dualuncued_power_SMA + std_dualuncued_power_SMA)], 'b',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, dualcued_power_SMA, '-m', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(dualcued_power_SMA - std_dualcued_power_SMA);...
    flipud(dualcued_power_SMA + std_dualcued_power_SMA)], 'm',...
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
plot(freq, singleuncued_power_SMA, '-g', 'LineWidth', 1.5);
hold on;
h1 = patch([freq; flipud(freq)],...
    [(singleuncued_power_SMA - std_singleuncued_power_SMA);...
    flipud(singleuncued_power_SMA + std_singleuncued_power_SMA)], 'g',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, singlecued_power_SMA, '-r', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(singlecued_power_SMA - std_singlecued_power_SMA);...
    flipud(singlecued_power_SMA + std_singlecued_power_SMA)], 'r',...
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
saveas(gcf, fullfile(results_path, 'SMA_PowervsFreq'),'png');

%% M1
dualuncued_power_C3 = dualuncued_power(:, C3_loc);
dualcued_power_C3 = dualcued_power(:, C3_loc);
singleuncued_power_C3 = singleuncued_power(:, C3_loc);
singlecued_power_C3 = singlecued_power(:, C3_loc);
std_dualuncued_power_C3 = std_dualuncued_power(:, C3_loc);
std_dualcued_power_C3 = std_dualcued_power(:, C3_loc);
std_singleuncued_power_C3 = std_singleuncued_power(:, C3_loc);
std_singlecued_power_C3 = std_singlecued_power(:, C3_loc);

dualuncued_power_CP1 = dualuncued_power(:, CP1_loc);
dualcued_power_CP1 = dualcued_power(:, CP1_loc);
singleuncued_power_CP1 = singleuncued_power(:, CP1_loc);
singlecued_power_CP1 = singlecued_power(:, CP1_loc);
std_dualuncued_power_CP1 = std_dualuncued_power(:, CP1_loc);
std_dualcued_power_CP1 = std_dualcued_power(:, CP1_loc);
std_singleuncued_power_CP1 = std_singleuncued_power(:, CP1_loc);
std_singlecued_power_CP1 = std_singlecued_power(:, CP1_loc);

dualuncued_power_M1 = mean([dualuncued_power_C3 dualuncued_power_CP1], 2);
dualcued_power_M1 = mean([dualcued_power_C3 dualcued_power_CP1], 2);
singleuncued_power_M1 = mean([singleuncued_power_C3 singleuncued_power_CP1], 2);
singlecued_power_M1 = mean([singlecued_power_C3 singlecued_power_CP1], 2);
std_dualuncued_power_M1 = mean([std_dualuncued_power_C3 std_dualuncued_power_CP1], 2);
std_dualcued_power_M1 = mean([std_dualcued_power_C3 std_dualcued_power_CP1], 2);
std_singleuncued_power_M1 = mean([std_singleuncued_power_C3 std_singleuncued_power_CP1], 2);
std_singlecued_power_M1 = mean([std_singlecued_power_C3 std_singlecued_power_CP1], 2);

figure; title('M1');
tiledlayout(1, 2);
% Dual Task
ax1 = nexttile;
plot(freq, dualuncued_power_M1, '-b', 'LineWidth', 1.5); 
hold on;
h1 = patch([freq; flipud(freq)],...
    [(dualuncued_power_M1 - std_dualuncued_power_M1);...
    flipud(dualuncued_power_M1 + std_dualuncued_power_M1)], 'b',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, dualcued_power_M1, '-m', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(dualcued_power_M1 - std_dualcued_power_M1);...
    flipud(dualcued_power_M1 + std_dualcued_power_M1)], 'm',...
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
plot(freq, singleuncued_power_M1, '-g', 'LineWidth', 1.5);
hold on;
h1 = patch([freq; flipud(freq)],...
    [(singleuncued_power_M1 - std_singleuncued_power_M1);...
    flipud(singleuncued_power_M1 + std_singleuncued_power_M1)], 'g',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, singlecued_power_M1, '-r', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(singlecued_power_M1 - std_singlecued_power_M1);...
    flipud(singlecued_power_M1 + std_singlecued_power_M1)], 'r',...
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
saveas(gcf, fullfile(results_path, 'M1_PowervsFreq'),'png');

%% PPC
dualuncued_power_P3 = dualuncued_power(:, P3_loc);
dualcued_power_P3 = dualcued_power(:, P3_loc);
singleuncued_power_P3 = singleuncued_power(:, P3_loc);
singlecued_power_P3 = singlecued_power(:, P3_loc);
std_dualuncued_power_P3 = std_dualuncued_power(:, P3_loc);
std_dualcued_power_P3 = std_dualcued_power(:, P3_loc);
std_singleuncued_power_P3 = std_singleuncued_power(:, P3_loc);
std_singlecued_power_P3 = std_singlecued_power(:, P3_loc);

dualuncued_power_P4 = dualuncued_power(:, P4_loc);
dualcued_power_P4 = dualcued_power(:, P4_loc);
singleuncued_power_P4 = singleuncued_power(:, P4_loc);
singlecued_power_P4 = singlecued_power(:, P4_loc);
std_dualuncued_power_P4 = std_dualuncued_power(:, P4_loc);
std_dualcued_power_P4 = std_dualcued_power(:, P4_loc);
std_singleuncued_power_P4 = std_singleuncued_power(:, P4_loc);
std_singlecued_power_P4 = std_singlecued_power(:, P4_loc);

dualuncued_power_PPC = mean([dualuncued_power_P3 dualuncued_power_P4], 2);
dualcued_power_PPC = mean([dualcued_power_P3 dualcued_power_P4], 2);
singleuncued_power_PPC = mean([singleuncued_power_P3 singleuncued_power_P4], 2);
singlecued_power_PPC = mean([singlecued_power_P3 singlecued_power_P4], 2);
std_dualuncued_power_PPC = mean([std_dualuncued_power_P3 std_dualuncued_power_P4], 2);
std_dualcued_power_PPC = mean([std_dualcued_power_P3 std_dualcued_power_P4], 2);
std_singleuncued_power_PPC = mean([std_singleuncued_power_P3 std_singleuncued_power_P4], 2);
std_singlecued_power_PPC = mean([std_singlecued_power_P3 std_singlecued_power_P4], 2);

figure; title('PPC');
tiledlayout(1, 2);
% Dual Task
ax1 = nexttile;
plot(freq, dualuncued_power_PPC, '-b', 'LineWidth', 1.5); 
hold on;
h1 = patch([freq; flipud(freq)],...
    [(dualuncued_power_PPC - std_dualuncued_power_PPC);...
    flipud(dualuncued_power_PPC + std_dualuncued_power_PPC)], 'b',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, dualcued_power_PPC, '-m', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(dualcued_power_PPC - std_dualcued_power_PPC);...
    flipud(dualcued_power_PPC + std_dualcued_power_PPC)], 'm',...
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
plot(freq, singleuncued_power_PPC, '-g', 'LineWidth', 1.5);
hold on;
h1 = patch([freq; flipud(freq)],...
    [(singleuncued_power_PPC - std_singleuncued_power_PPC);...
    flipud(singleuncued_power_PPC + std_singleuncued_power_PPC)], 'g',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, singlecued_power_PPC, '-r', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(singlecued_power_PPC - std_singlecued_power_PPC);...
    flipud(singlecued_power_PPC + std_singlecued_power_PPC)], 'r',...
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
saveas(gcf, fullfile(results_path, 'PPC_PowervsFreq'),'png');

%% Save the values onto all subs struct
avg.dualuncued_power_DLPFC = dualuncued_power_DLPFC;
avg.singleuncued_power_DLPFC = singleuncued_power_DLPFC;
avg.dualcued_power_DLPFC = dualcued_power_DLPFC;
avg.singlecued_power_DLPFC = singlecued_power_DLPFC;
avg.dualuncued_power_SMA = dualuncued_power_SMA;
avg.singleuncued_power_SMA = singleuncued_power_SMA;
avg.dualcued_power_SMA = dualcued_power_SMA;
avg.singlecued_power_SMA = singlecued_power_SMA;
avg.dualuncued_power_M1 = dualuncued_power_M1;
avg.singleuncued_power_M1 = singleuncued_power_M1;
avg.dualcued_power_M1 = dualcued_power_M1;
avg.singlecued_power_M1 = singlecued_power_M1;
avg.dualuncued_power_PPC = dualuncued_power_PPC;
avg.singleuncued_power_PPC = singleuncued_power_PPC;
avg.dualcued_power_PPC = dualcued_power_PPC;
avg.singlecued_power_PPC = singlecued_power_PPC;
allsubs.avg = avg;

std.dualuncued_power_DLPFC = std_dualuncued_power_DLPFC;
std.singleuncued_power_DLPFC = std_singleuncued_power_DLPFC;
std.dualcued_power_DLPFC = std_dualcued_power_DLPFC;
std.singlecued_power_DLPFC = std_singlecued_power_DLPFC;
std.dualuncued_power_SMA = std_dualuncued_power_SMA;
std.singleuncued_power_SMA = std_singleuncued_power_SMA;
std.dualcued_power_SMA = std_dualcued_power_SMA;
std.singlecued_power_SMA = std_singlecued_power_SMA;
std.dualuncued_power_M1 = std_dualuncued_power_M1;
std.singleuncued_power_M1 = std_singleuncued_power_M1;
std.dualcued_power_M1 = std_dualcued_power_M1;
std.singlecued_power_M1 = std_singlecued_power_M1;
std.dualuncued_power_PPC = std_dualuncued_power_PPC;
std.singleuncued_power_PPC = std_singleuncued_power_PPC;
std.dualcued_power_PPC = std_dualcued_power_PPC;
std.singlecued_power_PPC = std_singlecued_power_PPC;
allsubs.std = std;

% Save the struct from all subs
save(strcat(results_path, '\psdchannel_allsubs.mat'), 'allsubs')

disp('These are the results for the average of all subjects.');
