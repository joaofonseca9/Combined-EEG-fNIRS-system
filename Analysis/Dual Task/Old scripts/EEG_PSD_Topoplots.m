%% Analysis of the EEG signal (PSD: topoplots)
clear; clc; close all;

%% Initialize data
addpath('C:\Users\catar\OneDrive - Universidade do Porto\Twente\Combined-EEG-fNIRS-system\Analysis');
addpath('C:\Users\catar\OneDrive - Universidade do Porto\Twente\Combined-EEG-fNIRS-system\Analysis\Dual Task');
laptop = 'laptopCatarina';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);
results_path = 'C:\Users\catar\OneDrive - Universidade do Porto\Twente\Data Analysis\eeg\psd';

eeglab;

subrec = ["64" "01";"28" "04";"02" "02"];

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
    startTask = find(strcmp({EEG_DualUncued.event.type}, 's1709')==1);
    endTask = find(strcmp({EEG_DualUncued.event.type}, 's1717')==1);

    % Get the PSD averaged over all trials
    [power_theta, power_alpha, power_beta, power_gamma, freq_theta,...
        freq_alpha, freq_beta, freq_gamma] =...
        calculateAveragePowerBandsAllTrials(EEG_DualUncued, event_samp,...
        startTask, endTask);

    % Topographic distribution of the frequency bands over the head (topoplot)
    figure;
    subplot(2, 2, 1);
    text(-0.13, 0.7, 'Theta', 'FontSize', 18);
    topoplot(power_theta, EEG_DualUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 2);
    text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
    topoplot(power_alpha, EEG_DualUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 3);
    text(-0.1, 0.7, 'Beta', 'FontSize', 18)
    topoplot(power_beta, EEG_DualUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 4);
    text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
    topoplot(power_gamma, EEG_DualUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    
    if sub == "64"
       power_theta(31,:)=[]; 
       power_alpha(31,:)=[]; 
       power_beta(31,:)=[]; 
       power_gamma(31,:)=[]; 
    end
    if sub == "02"
       power_theta(30,:)=0; 
       power_alpha(30,:)=0; 
       power_beta(30,:)=0; 
       power_gamma(30,:)=0; 
       power_theta(29,:)=0; 
       power_alpha(29,:)=0; 
       power_beta(29,:)=0; 
       power_gamma(29,:)=0; 
    end
    
    % Save the values onto a allSubjects variable
    dualuncued_power_theta_allSubjects(:, subject) = power_theta;
    dualuncued_power_alpha_allSubjects(:, subject) = power_alpha;
    dualuncued_power_beta_allSubjects(:, subject) = power_beta;
    dualuncued_power_gamma_allSubjects(:, subject) = power_gamma;
    
    %% Single Uncued: PSD
    event_samp  = [EEG_SingleUncued.event.latency];
    startTask = find(strcmp({EEG_SingleUncued.event.type}, 's1705')==1);
    endTask = find(strcmp({EEG_SingleUncued.event.type}, 's1713')==1);
    
    % Get the PSD averaged over all trials
    [power_theta, power_alpha, power_beta, power_gamma, ~,...
        freq_alpha, freq_beta, freq_gamma] =...
        calculateAveragePowerBandsAllTrials(EEG_SingleUncued, event_samp,...
        startTask, endTask);
    
    % Topographic distribution of the frequency bands over the head (topoplot)
    figure;
    subplot(2, 2, 1);
    text(-0.13, 0.7, 'Theta', 'FontSize', 18);
    topoplot(power_theta, EEG_SingleUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 2);
    text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
    topoplot(power_alpha, EEG_SingleUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 3);
    text(-0.1, 0.7, 'Beta', 'FontSize', 18)
    topoplot(power_beta, EEG_SingleUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 4);
    text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
    topoplot(power_gamma, EEG_SingleUncued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    
    if sub == "64"
       power_theta(31,:)=[]; 
       power_alpha(31,:)=[]; 
       power_beta(31,:)=[]; 
       power_gamma(31,:)=[]; 
    end
    
    if sub == "02"
       power_theta(30,:)=0; 
       power_alpha(30,:)=0; 
       power_beta(30,:)=0; 
       power_gamma(30,:)=0; 
       power_theta(29,:)=0; 
       power_alpha(29,:)=0; 
       power_beta(29,:)=0; 
       power_gamma(29,:)=0; 
    end
    
    % Save the values onto a allSubjects variable
    singleuncued_power_theta_allSubjects(:, subject) = power_theta;
    singleuncued_power_alpha_allSubjects(:, subject) = power_alpha;
    singleuncued_power_beta_allSubjects(:, subject) = power_beta;
    singleuncued_power_gamma_allSubjects(:, subject) = power_gamma;
    
    %% Dual Cued: PSD
    event_samp  = [EEG_DualCued.event.latency];
    startTask = find(strcmp({EEG_DualCued.event.type}, 's1708')==1);
    endTask = find(strcmp({EEG_DualCued.event.type}, 's1701')==1);

    % Get the PSD averaged over all trials
    [power_theta, power_alpha, power_beta, power_gamma, freq_theta,...
        freq_alpha, freq_beta, freq_gamma] =...
        calculateAveragePowerBandsAllTrials(EEG_DualCued, event_samp,...
        startTask, endTask);

    % Topographic distribution of the frequency bands over the head (topoplot)
    figure;
    subplot(2, 2, 1);
    text(-0.13, 0.7, 'Theta', 'FontSize', 18);
    topoplot(power_theta, EEG_DualCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 2);
    text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
    topoplot(power_alpha, EEG_DualCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 3);
    text(-0.1, 0.7, 'Beta', 'FontSize', 18)
    topoplot(power_beta, EEG_DualCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;   
    subplot(2, 2, 4);
    text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
    topoplot(power_gamma, EEG_DualCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    
    if sub == "64"
       power_theta(31,:)=[]; 
       power_alpha(31,:)=[]; 
       power_beta(31,:)=[]; 
       power_gamma(31,:)=[]; 
    end
    
    if sub == "02"
       power_theta(30,:)=0; 
       power_alpha(30,:)=0; 
       power_beta(30,:)=0; 
       power_gamma(30,:)=0; 
       power_theta(29,:)=0; 
       power_alpha(29,:)=0; 
       power_beta(29,:)=0; 
       power_gamma(29,:)=0; 
    end
    
    % Save the values onto a allSubjects variable
    dualcued_power_theta_allSubjects(:, subject) = power_theta;
    dualcued_power_alpha_allSubjects(:, subject) = power_alpha;
    dualcued_power_beta_allSubjects(:, subject) = power_beta;
    dualcued_power_gamma_allSubjects(:, subject) = power_gamma;
    
    %% Single Cued: PSD
    event_samp  = [EEG_SingleCued.event.latency];
    startTask = find(strcmp({EEG_SingleCued.event.type}, 's1704')==1);
    endTask = find(strcmp({EEG_SingleCued.event.type}, 's1701')==1);

    % Get the PSD averaged over all trials 
    [power_theta, power_alpha, power_beta, power_gamma, freq_theta,...
        freq_alpha, freq_beta, freq_gamma] =...
        calculateAveragePowerBandsAllTrials(EEG_SingleCued, event_samp,...
        startTask, endTask);

    % Topographic distribution of the frequency bands over the head (topoplot)
    figure;
    subplot(2, 2, 1);
    text(-0.13, 0.7, 'Theta', 'FontSize', 18);
    topoplot(power_theta, EEG_SingleCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 2);
    text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
    topoplot(power_alpha, EEG_SingleCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    subplot(2, 2, 3);
    text(-0.1, 0.7, 'Beta', 'FontSize', 18)
    topoplot(power_beta, EEG_SingleCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar; 
    subplot(2, 2, 4);
    text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
    topoplot(power_gamma, EEG_SingleCued.chanlocs, 'electrodes', 'ptslabels');
    colorbar;
    
    if sub == "64"
       power_theta(31,:)=[]; 
       power_alpha(31,:)=[]; 
       power_beta(31,:)=[]; 
       power_gamma(31,:)=[]; 
    end
    
    if sub == "02"
       power_theta(30,:)=0; 
       power_alpha(30,:)=0; 
       power_beta(30,:)=0; 
       power_gamma(30,:)=0; 
       power_theta(29,:)=0; 
       power_alpha(29,:)=0; 
       power_beta(29,:)=0; 
       power_gamma(29,:)=0; 
    end
    
    % Save the values onto a allSubjects variable
    singlecued_power_theta_allSubjects(:, subject) = power_theta;
    singlecued_power_alpha_allSubjects(:, subject) = power_alpha;
    singlecued_power_beta_allSubjects(:, subject) = power_beta;
    singlecued_power_gamma_allSubjects(:, subject) = power_beta;
    
    disp(['These are the topoplots for subject ', char(sub), '.']);
    disp('Press any key to move onto the next subject.');
    %pause;
    close all;
    
end
disp('This was the end of individual subjects.');

%% Get the PSD averaged over all subjects
% Dual Uncued
dualuncued_power_theta = mean(dualuncued_power_theta_allSubjects, 2);
dualuncued_power_alpha = mean(dualuncued_power_alpha_allSubjects, 2);
dualuncued_power_beta = mean(dualuncued_power_beta_allSubjects, 2);
dualuncued_power_gamma = mean(dualuncued_power_gamma_allSubjects, 2);

% Single Uncued
singleuncued_power_theta = mean(singleuncued_power_theta_allSubjects, 2);
singleuncued_power_alpha = mean(singleuncued_power_alpha_allSubjects, 2);
singleuncued_power_beta = mean(singleuncued_power_beta_allSubjects, 2);
singleuncued_power_gamma = mean(singleuncued_power_gamma_allSubjects, 2);

% Dual Cued
dualcued_power_theta = mean(dualcued_power_theta_allSubjects, 2);
dualcued_power_alpha = mean(dualcued_power_alpha_allSubjects, 2);
dualcued_power_beta = mean(dualcued_power_beta_allSubjects, 2);
dualcued_power_gamma = mean(dualcued_power_gamma_allSubjects, 2);

% Single Cued
singlecued_power_theta = mean(singlecued_power_theta_allSubjects, 2);
singlecued_power_alpha = mean(singlecued_power_alpha_allSubjects, 2);
singlecued_power_beta = mean(singlecued_power_beta_allSubjects, 2);
singlecued_power_gamma = mean(singlecued_power_gamma_allSubjects, 2);

%% Averaged topographic distribution of the frequency bands over the head (topoplot)
% Dual Uncued
figure;
subplot(2, 2, 1);
text(-0.13, 0.7, 'Theta', 'FontSize', 18);
topoplot(dualuncued_power_theta, EEG_DualUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 2);
text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
topoplot(dualuncued_power_alpha, EEG_DualUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 3);
text(-0.1, 0.7, 'Beta', 'FontSize', 18)
topoplot(dualuncued_power_beta, EEG_DualUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar; 
subplot(2, 2, 4);
text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
topoplot(dualuncued_power_gamma, EEG_DualUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar;

set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, fullfile(results_path, 'dualuncued_psdtopoplot'),'png');

% Single Uncued
figure;
subplot(2, 2, 1);
text(-0.13, 0.7, 'Theta', 'FontSize', 18);
topoplot(singleuncued_power_theta, EEG_SingleUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 2);
text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
topoplot(singleuncued_power_alpha, EEG_SingleUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 3);
text(-0.1, 0.7, 'Beta', 'FontSize', 18)
topoplot(singleuncued_power_beta, EEG_SingleUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar; 
subplot(2, 2, 4);
text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
topoplot(singleuncued_power_gamma, EEG_SingleUncued.chanlocs, 'electrodes', 'ptslabels');
colorbar;

set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, fullfile(results_path, 'singleuncued_psdtopoplot'),'png');

% Dual Cued
figure;
subplot(2, 2, 1);
text(-0.13, 0.7, 'Theta', 'FontSize', 18);
topoplot(dualcued_power_theta, EEG_DualCued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 2);
text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
topoplot(dualcued_power_alpha, EEG_DualCued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 3);
text(-0.1, 0.7, 'Beta', 'FontSize', 18)
topoplot(dualcued_power_beta, EEG_DualCued.chanlocs, 'electrodes', 'ptslabels');
colorbar; 
subplot(2, 2, 4);
text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
topoplot(dualcued_power_gamma, EEG_DualCued.chanlocs, 'electrodes', 'ptslabels');
colorbar;

set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, fullfile(results_path, 'dualcued_psdtopoplot'),'png');

% Single Cued
figure;
subplot(2, 2, 1);
text(-0.13, 0.7, 'Theta', 'FontSize', 18);
topoplot(singlecued_power_theta, EEG_SingleCued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 2);
text(-0.13, 0.7, 'Alpha', 'FontSize', 18)
topoplot(singlecued_power_alpha, EEG_SingleCued.chanlocs, 'electrodes', 'ptslabels');
colorbar;
subplot(2, 2, 3);
text(-0.1, 0.7, 'Beta', 'FontSize', 18)
topoplot(singlecued_power_beta, EEG_SingleCued.chanlocs, 'electrodes', 'ptslabels');
colorbar; 
subplot(2, 2, 4);
text(-0.2, 0.7, 'Gamma', 'FontSize', 18)
topoplot(singlecued_power_gamma, EEG_SingleCued.chanlocs, 'electrodes', 'ptslabels');
colorbar;

set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, fullfile(results_path, 'singlecued_psdtopoplot'),'png');

disp('These are the topoplots for the average of all subjects.');
