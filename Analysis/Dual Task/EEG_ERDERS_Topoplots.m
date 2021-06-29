%% Analysis of the EEG signal (ERD/ERS: topoplots)
clear; clc; close all;

%% Initialize data
addpath('C:\Users\catar\OneDrive - Universidade do Porto\Twente\Combined-EEG-fNIRS-system\Analysis');
addpath('C:\Users\catar\OneDrive - Universidade do Porto\Twente\Combined-EEG-fNIRS-system\Analysis\Dual Task');
laptop = 'laptopCatarina';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);
results_path = 'C:\Users\catar\OneDrive - Universidade do Porto\Twente\Data Analysis\eeg\erders';

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
    
    %% Dual Uncued: ERD/ERS
    event_samp  = [EEG_DualUncued.event.latency];
    startBaseline = find(strcmp({EEG_DualUncued.event.type}, 'boundary')==1);
    startTask = find(strcmp({EEG_DualUncued.event.type}, 's1798')==1);
    endTask = find(strcmp({EEG_DualUncued.event.type}, 's1717')==1);
   
    % Remove any extra boundaries
    startBaseline = removeExtraBoundaries(startBaseline, startTask);

    % Get the PSD averaged over all trials
    % For the baseline
    [power_base_theta, power_base_alpha, power_base_beta,...
        power_base_gamma, freq_base_theta, freq_base_alpha,...
        freq_base_beta, freq_base_gamma] =...
        calculateAveragePowerERDERSAllTrials(EEG_DualUncued, event_samp,...
        startBaseline, startTask);
    % For the task
    [power_theta, power_alpha, power_beta, power_gamma, freq_theta,...
        freq_alpha, freq_beta, freq_gamma] =...
        calculateAveragePowerERDERSAllTrials(EEG_DualUncued, event_samp,...
        startTask, endTask);
    
    % Calculate the ERD/ERS for each of the frequency bands above
    ERD_ERS_theta = (power_theta - power_base_theta)./power_base_theta; 
    ERD_ERS_alpha = (power_alpha - power_base_alpha)./power_base_alpha; 
    ERD_ERS_beta = (power_beta - power_base_beta)./power_base_beta; 
    ERD_ERS_gamma = (power_gamma - power_base_gamma)./power_base_gamma; 

    % Topographic distribution of the frequency bands over the head (topoplot)
    figure;
    subplot(2, 2, 1);
    text(-0.13, 0.7, 'Theta', 'FontSize', 18);
    topoplot(ERD_ERS_theta, EEG_DualUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 2);
    text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
    topoplot(ERD_ERS_alpha, EEG_DualUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 3);
    text(-0.1, 0.7, 'Beta', 'FontSize', 18)
    topoplot(ERD_ERS_beta, EEG_DualUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 4);
    text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
    topoplot(ERD_ERS_gamma, EEG_DualUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    
    % Compensate for removed channels
    ERD_ERS_theta = compensateRemovedChannels(ERD_ERS_theta, EEG_DualUncued, list_channels);
    ERD_ERS_alpha = compensateRemovedChannels(ERD_ERS_alpha, EEG_DualUncued, list_channels);
    ERD_ERS_beta = compensateRemovedChannels(ERD_ERS_beta, EEG_DualUncued, list_channels);
    ERD_ERS_gamma = compensateRemovedChannels(ERD_ERS_gamma, EEG_DualUncued, list_channels);
    
    % Save the values onto a allSubjects variable
    dualuncued_ERD_ERS_theta_allSubjects(:, subject) = ERD_ERS_theta;
    dualuncued_ERD_ERS_alpha_allSubjects(:, subject) = ERD_ERS_alpha;
    dualuncued_ERD_ERS_beta_allSubjects(:, subject) = ERD_ERS_beta;
    dualuncued_ERD_ERS_gamma_allSubjects(:, subject) = ERD_ERS_gamma;
    
    %% Single Uncued: ERD/ERS
    event_samp  = [EEG_SingleUncued.event.latency];
    startBaseline = find(strcmp({EEG_SingleUncued.event.type}, 'boundary')==1);
    startTask = find(strcmp({EEG_SingleUncued.event.type}, 's1705')==1);
    endTask = find(strcmp({EEG_SingleUncued.event.type}, 's1713')==1);
    
    % Remove any extra boundaries
    startBaseline = removeExtraBoundaries(startBaseline, startTask);
    
    % Get the PSD averaged over all trials
    % For the baseline
    [power_base_theta, power_base_alpha, power_base_beta,...
        power_base_gamma, freq_base_theta, freq_base_alpha,...
        freq_base_beta, freq_base_gamma] =...
        calculateAveragePowerERDERSAllTrials(EEG_SingleUncued, event_samp,...
        startBaseline, startTask);
    % For the task
    [power_theta, power_alpha, power_beta, power_gamma, freq_theta,...
        freq_alpha, freq_beta, freq_gamma] =...
        calculateAveragePowerERDERSAllTrials(EEG_SingleUncued, event_samp,...
        startTask, endTask);
    
    % Calculate the ERD/ERS for each of the frequency bands above
    ERD_ERS_theta = (power_theta - power_base_theta)./power_base_theta; 
    ERD_ERS_alpha = (power_alpha - power_base_alpha)./power_base_alpha; 
    ERD_ERS_beta = (power_beta - power_base_beta)./power_base_beta; 
    ERD_ERS_gamma = (power_gamma - power_base_gamma)./power_base_gamma; 
    
    % Topographic distribution of the frequency bands over the head (topoplot)
    figure;
    subplot(2, 2, 1);
    text(-0.13, 0.7, 'Theta', 'FontSize', 18);
    topoplot(ERD_ERS_theta, EEG_SingleUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 2);
    text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
    topoplot(ERD_ERS_alpha, EEG_SingleUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 3);
    text(-0.1, 0.7, 'Beta', 'FontSize', 18)
    topoplot(ERD_ERS_beta, EEG_SingleUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 4);
    text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
    topoplot(ERD_ERS_gamma, EEG_SingleUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    
    % Compensate for removed channels
    ERD_ERS_theta = compensateRemovedChannels(ERD_ERS_theta, EEG_SingleUncued, list_channels);
    ERD_ERS_alpha = compensateRemovedChannels(ERD_ERS_alpha, EEG_SingleUncued, list_channels);
    ERD_ERS_beta = compensateRemovedChannels(ERD_ERS_beta, EEG_SingleUncued, list_channels);
    ERD_ERS_gamma = compensateRemovedChannels(ERD_ERS_gamma, EEG_SingleUncued, list_channels);
    
    % Save the values onto a allSubjects variable
    singleuncued_ERD_ERS_theta_allSubjects(:, subject) = ERD_ERS_theta;
    singleuncued_ERD_ERS_alpha_allSubjects(:, subject) = ERD_ERS_alpha;
    singleuncued_ERD_ERS_beta_allSubjects(:, subject) = ERD_ERS_beta;
    singleuncued_ERD_ERS_gamma_allSubjects(:, subject) = ERD_ERS_gamma;
       
    %% Dual Cued: ERD/ERS
    event_samp  = [EEG_DualCued.event.latency];
    startBaseline = find(strcmp({EEG_DualCued.event.type}, 'boundary')==1);
    startTask = find(strcmp({EEG_DualCued.event.type}, 's1798')==1);
    endTask = find(strcmp({EEG_DualCued.event.type}, 's1701')==1);

    % Remove any extra boundaries
    startBaseline = removeExtraBoundaries(startBaseline, startTask);

    % Get the PSD averaged over all trials
    % For the baseline
    [power_base_theta, power_base_alpha, power_base_beta,...
        power_base_gamma, freq_base_theta, freq_base_alpha,...
        freq_base_beta, freq_base_gamma] =...
        calculateAveragePowerERDERSAllTrials(EEG_DualCued, event_samp,...
        startBaseline, startTask);
    % For the task
    [power_theta, power_alpha, power_beta, power_gamma, freq_theta,...
        freq_alpha, freq_beta, freq_gamma] =...
        calculateAveragePowerERDERSAllTrials(EEG_DualCued, event_samp,...
        startTask, endTask);
    
    % Calculate the ERD/ERS for each of the frequency bands above
    ERD_ERS_theta = (power_theta - power_base_theta)./power_base_theta; 
    ERD_ERS_alpha = (power_alpha - power_base_alpha)./power_base_alpha; 
    ERD_ERS_beta = (power_beta - power_base_beta)./power_base_beta; 
    ERD_ERS_gamma = (power_gamma - power_base_gamma)./power_base_gamma; 
    
    % Topographic distribution of the frequency bands over the head (topoplot)
    figure;
    subplot(2, 2, 1);
    text(-0.13, 0.7, 'Theta', 'FontSize', 18);
    topoplot(ERD_ERS_theta, EEG_DualCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 2);
    text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
    topoplot(ERD_ERS_alpha, EEG_DualCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 3);
    text(-0.1, 0.7, 'Beta', 'FontSize', 18)
    topoplot(ERD_ERS_beta, EEG_DualCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;   
    subplot(2, 2, 4);
    text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
    topoplot(ERD_ERS_gamma, EEG_DualCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    
    % Compensate for removed channels
    ERD_ERS_theta = compensateRemovedChannels(ERD_ERS_theta, EEG_DualCued, list_channels);
    ERD_ERS_alpha = compensateRemovedChannels(ERD_ERS_alpha, EEG_DualCued, list_channels);
    ERD_ERS_beta = compensateRemovedChannels(ERD_ERS_beta, EEG_DualCued, list_channels);
    ERD_ERS_gamma = compensateRemovedChannels(ERD_ERS_gamma, EEG_DualCued, list_channels);
    
    % Save the values onto a allSubjects variable
    dualcued_ERD_ERS_theta_allSubjects(:, subject) = ERD_ERS_theta;
    dualcued_ERD_ERS_alpha_allSubjects(:, subject) = ERD_ERS_alpha;
    dualcued_ERD_ERS_beta_allSubjects(:, subject) = ERD_ERS_beta;
    dualcued_ERD_ERS_gamma_allSubjects(:, subject) = ERD_ERS_gamma;
    
    %% Single Cued: ERD/ERS
    event_samp  = [EEG_SingleCued.event.latency];
    startBaseline = find(strcmp({EEG_SingleCued.event.type}, 'boundary')==1);
    startTask = find(strcmp({EEG_SingleCued.event.type}, 's1704')==1);
    endTask = find(strcmp({EEG_SingleCued.event.type}, 's1701')==1);

    % Remove any extra boundaries
    startBaseline = removeExtraBoundaries(startBaseline, startTask);

    % Get the PSD averaged over all trials 
    % For the baseline
    [power_base_theta, power_base_alpha, power_base_beta,...
        power_base_gamma, freq_base_theta, freq_base_alpha,...
        freq_base_beta, freq_base_gamma] =...
        calculateAveragePowerERDERSAllTrials(EEG_SingleCued, event_samp,...
        startBaseline, startTask);
    % For the task
    [power_theta, power_alpha, power_beta, power_gamma, freq_theta,...
        freq_alpha, freq_beta, freq_gamma] =...
        calculateAveragePowerERDERSAllTrials(EEG_SingleCued, event_samp,...
        startTask, endTask);
    
    % Calculate the ERD/ERS for each of the frequency bands above
    ERD_ERS_theta = (power_theta - power_base_theta)./power_base_theta; 
    ERD_ERS_alpha = (power_alpha - power_base_alpha)./power_base_alpha; 
    ERD_ERS_beta = (power_beta - power_base_beta)./power_base_beta; 
    ERD_ERS_gamma = (power_gamma - power_base_gamma)./power_base_gamma; 

    % Topographic distribution of the frequency bands over the head (topoplot)
    figure;
    subplot(2, 2, 1);
    text(-0.13, 0.7, 'Theta', 'FontSize', 18);
    topoplot(ERD_ERS_theta, EEG_SingleCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 2);
    text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
    topoplot(ERD_ERS_alpha, EEG_SingleCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 3);
    text(-0.1, 0.7, 'Beta', 'FontSize', 18)
    topoplot(ERD_ERS_beta, EEG_SingleCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar; 
    subplot(2, 2, 4);
    text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
    topoplot(ERD_ERS_gamma, EEG_SingleCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    
    % Compensate for removed channels
    ERD_ERS_theta = compensateRemovedChannels(ERD_ERS_theta, EEG_SingleCued, list_channels);
    ERD_ERS_alpha = compensateRemovedChannels(ERD_ERS_alpha, EEG_SingleCued, list_channels);
    ERD_ERS_beta = compensateRemovedChannels(ERD_ERS_beta, EEG_SingleCued, list_channels);
    ERD_ERS_gamma = compensateRemovedChannels(ERD_ERS_gamma, EEG_SingleCued, list_channels);
    
    % Save the values onto a allSubjects variable
    singlecued_ERD_ERS_theta_allSubjects(:, subject) = ERD_ERS_theta;
    singlecued_ERD_ERS_alpha_allSubjects(:, subject) = ERD_ERS_alpha;
    singlecued_ERD_ERS_beta_allSubjects(:, subject) = ERD_ERS_beta;
    singlecued_ERD_ERS_gamma_allSubjects(:, subject) = ERD_ERS_beta;
    
    disp(['These are the topoplots for subject ', char(sub), '.']);
    disp('Press any key to move onto the next subject.');
    pause;
    close all;
    
end

%% Get the PSD averaged over all subjects
% Dual Uncued
dualuncued_ERD_ERS_theta = mean(dualuncued_ERD_ERS_theta_allSubjects, 2, 'omitnan');
dualuncued_ERD_ERS_alpha = mean(dualuncued_ERD_ERS_alpha_allSubjects, 2, 'omitnan');
dualuncued_ERD_ERS_beta = mean(dualuncued_ERD_ERS_beta_allSubjects, 2, 'omitnan');
dualuncued_ERD_ERS_gamma = mean(dualuncued_ERD_ERS_gamma_allSubjects, 2, 'omitnan');

% Single Uncued
singleuncued_ERD_ERS_theta = mean(singleuncued_ERD_ERS_theta_allSubjects, 2, 'omitnan');
singleuncued_ERD_ERS_alpha = mean(singleuncued_ERD_ERS_alpha_allSubjects, 2, 'omitnan');
singleuncued_ERD_ERS_beta = mean(singleuncued_ERD_ERS_beta_allSubjects, 2, 'omitnan');
singleuncued_ERD_ERS_gamma = mean(singleuncued_ERD_ERS_gamma_allSubjects, 2, 'omitnan');

% Dual Cued
dualcued_ERD_ERS_theta = mean(dualcued_ERD_ERS_theta_allSubjects, 2, 'omitnan');
dualcued_ERD_ERS_alpha = mean(dualcued_ERD_ERS_alpha_allSubjects, 2, 'omitnan');
dualcued_ERD_ERS_beta = mean(dualcued_ERD_ERS_beta_allSubjects, 2, 'omitnan');
dualcued_ERD_ERS_gamma = mean(dualcued_ERD_ERS_gamma_allSubjects, 2, 'omitnan');

% Single Cued
singlecued_ERD_ERS_theta = mean(singlecued_ERD_ERS_theta_allSubjects, 2, 'omitnan');
singlecued_ERD_ERS_alpha = mean(singlecued_ERD_ERS_alpha_allSubjects, 2, 'omitnan');
singlecued_ERD_ERS_beta = mean(singlecued_ERD_ERS_beta_allSubjects, 2, 'omitnan');
singlecued_ERD_ERS_gamma = mean(singlecued_ERD_ERS_gamma_allSubjects, 2, 'omitnan');

%% Averaged topographic distribution of the frequency bands over the head (topoplot)
% Dual Uncued
figure;
subplot(2, 2, 1);
text(-0.13, 0.7, 'Theta', 'FontSize', 18);
topoplot(dualuncued_ERD_ERS_theta, EEG_DualUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 2);
text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
topoplot(dualuncued_ERD_ERS_alpha, EEG_DualUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 3);
text(-0.1, 0.7, 'Beta', 'FontSize', 18)
topoplot(dualuncued_ERD_ERS_beta, EEG_DualUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar; 
subplot(2, 2, 4);
text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
topoplot(dualuncued_ERD_ERS_gamma, EEG_DualUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
% Save figure
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, fullfile(results_path, 'dualuncued_erders'),'png');

% Single Uncued
figure;
subplot(2, 2, 1);
text(-0.13, 0.7, 'Theta', 'FontSize', 18);
topoplot(singleuncued_ERD_ERS_theta, EEG_SingleUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 2);
text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
topoplot(singleuncued_ERD_ERS_alpha, EEG_SingleUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 3);
text(-0.1, 0.7, 'Beta', 'FontSize', 18)
topoplot(singleuncued_ERD_ERS_beta, EEG_SingleUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar; 
subplot(2, 2, 4);
text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
topoplot(singleuncued_ERD_ERS_gamma, EEG_SingleUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
% Save figure
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, fullfile(results_path, 'singleuncued_erders'),'png');

% Dual Cued
figure;
subplot(2, 2, 1);
text(-0.13, 0.7, 'Theta', 'FontSize', 18);
topoplot(dualcued_ERD_ERS_theta, EEG_DualCued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 2);
text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
topoplot(dualcued_ERD_ERS_alpha, EEG_DualCued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 3);
text(-0.1, 0.7, 'Beta', 'FontSize', 18)
topoplot(dualcued_ERD_ERS_beta, EEG_DualCued.chanlocs, 'electrodes', 'ptslabels');
colorbar; 
subplot(2, 2, 4);
text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
topoplot(dualcued_ERD_ERS_gamma, EEG_DualCued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
% Save figure
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, fullfile(results_path, 'dualcued_erders'),'png');

% Single Cued
figure;
subplot(2, 2, 1);
text(-0.13, 0.7, 'Theta', 'FontSize', 18);
topoplot(singlecued_ERD_ERS_theta, EEG_SingleCued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 2);
text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
topoplot(singlecued_ERD_ERS_alpha, EEG_SingleCued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 3);
text(-0.1, 0.7, 'Beta', 'FontSize', 18)
topoplot(singlecued_ERD_ERS_beta, EEG_SingleCued.chanlocs, 'electrodes', 'ptslabels');
colorbar; 
subplot(2, 2, 4);
text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
topoplot(singlecued_ERD_ERS_gamma, EEG_SingleCued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
% Save figure
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, fullfile(results_path, 'singlecued_erders'),'png');

disp('This was the end of individual subjects.');
disp('These are the topoplots for the average of all subjects.');
