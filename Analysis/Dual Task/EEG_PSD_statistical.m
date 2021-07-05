%% Statistical Analysis: EEG (PSD - Channel Analysis)
clear; clc; close all;

%% Initialize data
addpath('C:\Users\catar\OneDrive - Universidade do Porto\Twente\Combined-EEG-fNIRS-system\Analysis');
addpath('C:\Users\catar\OneDrive - Universidade do Porto\Twente\Combined-EEG-fNIRS-system\Analysis\Dual Task');
laptop = 'laptopCatarina';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);
results_path = 'C:\Users\catar\OneDrive - Universidade do Porto\Twente\Data Analysis\eeg\psd';
statistics_path = 'C:\Users\catar\OneDrive - Universidade do Porto\Twente\Data Analysis\statistics\eeg\psd (channel analysis)';

load(fullfile(results_path, 'psdchannel_allsubs.mat'), 'allsubs');

subrec = ["28" "04"; "02" "02"; "76" "01"; "64" "01"];

for subject=1:size(subrec, 1)
    sub = subrec(subject, 1);
    rec = subrec(subject, 2);
    
    sub_struct = allsubs.(genvarname(strcat('sub', char(sub))));
    
    %% Get average value over frequency bands for each region - dual
    % DLPFC
    dualuncued_DLPFC = sub_struct.dualuncued_power_DLPFC;
    dualcued_DLPFC = sub_struct.dualcued_power_DLPFC;
    
    % Theta band - 4:0.25:8
    dual_DLPFC_theta(subject, 1) = mean(dualuncued_DLPFC(1:length(4:0.25:8)));
    dual_DLPFC_theta(subject, 2) = mean(dualcued_DLPFC(1:length(4:0.25:8)));
    % Alpha band - 8:0.25:13
    dual_DLPFC_alpha(subject, 1) = mean(dualuncued_DLPFC(length(4:0.25:8):length(4:0.25:13)));
    dual_DLPFC_alpha(subject, 2) = mean(dualcued_DLPFC(length(4:0.25:8):length(4:0.25:13)));
    % Beta band - 13:0.25:32
    dual_DLPFC_beta(subject, 1) = mean(dualuncued_DLPFC(length(4:0.25:13):length(4:0.25:32)));
    dual_DLPFC_beta(subject, 2) = mean(dualcued_DLPFC(length(4:0.25:13):length(4:0.25:32)));
    % Gamma band - 32:0.25:48
    dual_DLPFC_gamma(subject, 1) = mean(dualuncued_DLPFC(length(4:0.25:32):length(4:0.25:48)));
    dual_DLPFC_gamma(subject, 2) = mean(dualcued_DLPFC(length(4:0.25:32):length(4:0.25:48)));
    
    % SMA
    dualuncued_SMA = sub_struct.dualuncued_power_SMA;
    dualcued_SMA = sub_struct.dualcued_power_SMA;
    
    % Theta band - 4:0.25:8
    dual_SMA_theta(subject, 1) = mean(dualuncued_SMA(1:length(4:0.25:8)));
    dual_SMA_theta(subject, 2) = mean(dualcued_SMA(1:length(4:0.25:8)));
    % Alpha band - 8:0.25:13
    dual_SMA_alpha(subject, 1) = mean(dualuncued_SMA(length(4:0.25:8):length(4:0.25:13)));
    dual_SMA_alpha(subject, 2) = mean(dualcued_SMA(length(4:0.25:8):length(4:0.25:13)));
    % Beta band - 13:0.25:32
    dual_SMA_beta(subject, 1) = mean(dualuncued_SMA(length(4:0.25:13):length(4:0.25:32)));
    dual_SMA_beta(subject, 2) = mean(dualcued_SMA(length(4:0.25:13):length(4:0.25:32)));
    % Gamma band - 32:0.25:48
    dual_SMA_gamma(subject, 1) = mean(dualuncued_SMA(length(4:0.25:32):length(4:0.25:48)));
    dual_SMA_gamma(subject, 2) = mean(dualcued_SMA(length(4:0.25:32):length(4:0.25:48)));
    
    % M1
    dualuncued_M1 = sub_struct.dualuncued_power_M1;
    dualcued_M1 = sub_struct.dualcued_power_M1;
    
    % Theta band - 4:0.25:8
    dual_M1_theta(subject, 1) = mean(dualuncued_M1(1:length(4:0.25:8)));
    dual_M1_theta(subject, 2) = mean(dualcued_M1(1:length(4:0.25:8)));
    % Alpha band - 8:0.25:13
    dual_M1_alpha(subject, 1) = mean(dualuncued_M1(length(4:0.25:8):length(4:0.25:13)));
    dual_M1_alpha(subject, 2) = mean(dualcued_M1(length(4:0.25:8):length(4:0.25:13)));
    % Beta band - 13:0.25:32
    dual_M1_beta(subject, 1) = mean(dualuncued_M1(length(4:0.25:13):length(4:0.25:32)));
    dual_M1_beta(subject, 2) = mean(dualcued_M1(length(4:0.25:13):length(4:0.25:32)));
    % Gamma band - 32:0.25:48
    dual_M1_gamma(subject, 1) = mean(dualuncued_M1(length(4:0.25:32):length(4:0.25:48)));
    dual_M1_gamma(subject, 2) = mean(dualcued_M1(length(4:0.25:32):length(4:0.25:48)));
    
    % PPC
    dualuncued_PPC = sub_struct.dualuncued_power_PPC;
    dualcued_PPC = sub_struct.dualcued_power_PPC;
    
    % Theta band - 4:0.25:8
    dual_PPC_theta(subject, 1) = mean(dualuncued_PPC(1:length(4:0.25:8)));
    dual_PPC_theta(subject, 2) = mean(dualcued_PPC(1:length(4:0.25:8)));
    % Alpha band - 8:0.25:13
    dual_PPC_alpha(subject, 1) = mean(dualuncued_PPC(length(4:0.25:8):length(4:0.25:13)));
    dual_PPC_alpha(subject, 2) = mean(dualcued_PPC(length(4:0.25:8):length(4:0.25:13)));
    % Beta band - 13:0.25:32
    dual_PPC_beta(subject, 1) = mean(dualuncued_PPC(length(4:0.25:13):length(4:0.25:32)));
    dual_PPC_beta(subject, 2) = mean(dualcued_PPC(length(4:0.25:13):length(4:0.25:32)));
    % Gamma band - 32:0.25:48
    dual_PPC_gamma(subject, 1) = mean(dualuncued_PPC(length(4:0.25:32):length(4:0.25:48)));
    dual_PPC_gamma(subject, 2) = mean(dualcued_PPC(length(4:0.25:32):length(4:0.25:48)));
    
    %% Get average value over frequency bands for each region - uncued
    % DLPFC
    singlecued_DLPFC = sub_struct.singleuncued_power_DLPFC;
    
    % Theta band - 4:0.25:8
    uncued_DLPFC_theta(subject, 1) = mean(singlecued_DLPFC(1:length(4:0.25:8)));
    uncued_DLPFC_theta(subject, 2) = mean(dualuncued_DLPFC(1:length(4:0.25:8)));
    % Alpha band - 8:0.25:13
    uncued_DLPFC_alpha(subject, 1) = mean(singlecued_DLPFC(length(4:0.25:8):length(4:0.25:13)));
    uncued_DLPFC_alpha(subject, 2) = mean(dualuncued_DLPFC(length(4:0.25:8):length(4:0.25:13)));
    % Beta band - 13:0.25:32
    uncued_DLPFC_beta(subject, 1) = mean(singlecued_DLPFC(length(4:0.25:13):length(4:0.25:32)));
    uncued_DLPFC_beta(subject, 2) = mean(dualuncued_DLPFC(length(4:0.25:13):length(4:0.25:32)));
    % Gamma band - 32:0.25:48
    uncued_DLPFC_gamma(subject, 1) = mean(singlecued_DLPFC(length(4:0.25:32):length(4:0.25:48)));
    uncued_DLPFC_gamma(subject, 2) = mean(dualuncued_DLPFC(length(4:0.25:32):length(4:0.25:48)));
    
    % SMA
    singlecued_SMA = sub_struct.singleuncued_power_SMA;
    
    % Theta band - 4:0.25:8
    uncued_SMA_theta(subject, 1) = mean(singlecued_SMA(1:length(4:0.25:8)));
    uncued_SMA_theta(subject, 2) = mean(dualuncued_SMA(1:length(4:0.25:8)));
    % Alpha band - 8:0.25:13
    uncued_SMA_alpha(subject, 1) = mean(singlecued_SMA(length(4:0.25:8):length(4:0.25:13)));
    uncued_SMA_alpha(subject, 2) = mean(dualuncued_SMA(length(4:0.25:8):length(4:0.25:13)));
    % Beta band - 13:0.25:32
    uncued_SMA_beta(subject, 1) = mean(singlecued_SMA(length(4:0.25:13):length(4:0.25:32)));
    uncued_SMA_beta(subject, 2) = mean(dualuncued_SMA(length(4:0.25:13):length(4:0.25:32)));
    % Gamma band - 32:0.25:48
    uncued_SMA_gamma(subject, 1) = mean(singlecued_SMA(length(4:0.25:32):length(4:0.25:48)));
    uncued_SMA_gamma(subject, 2) = mean(dualuncued_SMA(length(4:0.25:32):length(4:0.25:48)));
    
    % M1
    singlecued_M1 = sub_struct.singleuncued_power_M1;
    
    % Theta band - 4:0.25:8
    uncued_M1_theta(subject, 1) = mean(singlecued_M1(1:length(4:0.25:8)));
    uncued_M1_theta(subject, 2) = mean(dualuncued_M1(1:length(4:0.25:8)));
    % Alpha band - 8:0.25:13
    uncued_M1_alpha(subject, 1) = mean(singlecued_M1(length(4:0.25:8):length(4:0.25:13)));
    uncued_M1_alpha(subject, 2) = mean(dualuncued_M1(length(4:0.25:8):length(4:0.25:13)));
    % Beta band - 13:0.25:32
    uncued_M1_beta(subject, 1) = mean(singlecued_M1(length(4:0.25:13):length(4:0.25:32)));
    uncued_M1_beta(subject, 2) = mean(dualuncued_M1(length(4:0.25:13):length(4:0.25:32)));
    % Gamma band - 32:0.25:48
    uncued_M1_gamma(subject, 1) = mean(singlecued_M1(length(4:0.25:32):length(4:0.25:48)));
    uncued_M1_gamma(subject, 2) = mean(dualuncued_M1(length(4:0.25:32):length(4:0.25:48)));
    
    % PPC
    singlecued_PPC = sub_struct.singleuncued_power_PPC;
    
    % Theta band - 4:0.25:8
    uncued_PPC_theta(subject, 1) = mean(singlecued_PPC(1:length(4:0.25:8)));
    uncued_PPC_theta(subject, 2) = mean(dualuncued_PPC(1:length(4:0.25:8)));
    % Alpha band - 8:0.25:13
    uncued_PPC_alpha(subject, 1) = mean(singlecued_PPC(length(4:0.25:8):length(4:0.25:13)));
    uncued_PPC_alpha(subject, 2) = mean(dualuncued_PPC(length(4:0.25:8):length(4:0.25:13)));
    % Beta band - 13:0.25:32
    uncued_PPC_beta(subject, 1) = mean(singlecued_PPC(length(4:0.25:13):length(4:0.25:32)));
    uncued_PPC_beta(subject, 2) = mean(dualuncued_PPC(length(4:0.25:13):length(4:0.25:32)));
    % Gamma band - 32:0.25:48
    uncued_PPC_gamma(subject, 1) = mean(singlecued_PPC(length(4:0.25:32):length(4:0.25:48)));
    uncued_PPC_gamma(subject, 2) = mean(dualuncued_PPC(length(4:0.25:32):length(4:0.25:48)));
    
    %% Get average value over frequency bands for each region - cued
    % DLPFC
    singlecued_DLPFC = sub_struct.singlecued_power_DLPFC;
    
    % Theta band - 4:0.25:8
    cued_DLPFC_theta(subject, 1) = mean(singlecued_DLPFC(1:length(4:0.25:8)));
    cued_DLPFC_theta(subject, 2) = mean(dualcued_DLPFC(1:length(4:0.25:8)));
    % Alpha band - 8:0.25:13
    cued_DLPFC_alpha(subject, 1) = mean(singlecued_DLPFC(length(4:0.25:8):length(4:0.25:13)));
    cued_DLPFC_alpha(subject, 2) = mean(dualcued_DLPFC(length(4:0.25:8):length(4:0.25:13)));
    % Beta band - 13:0.25:32
    cued_DLPFC_beta(subject, 1) = mean(singlecued_DLPFC(length(4:0.25:13):length(4:0.25:32)));
    cued_DLPFC_beta(subject, 2) = mean(dualcued_DLPFC(length(4:0.25:13):length(4:0.25:32)));
    % Gamma band - 32:0.25:48
    cued_DLPFC_gamma(subject, 1) = mean(singlecued_DLPFC(length(4:0.25:32):length(4:0.25:48)));
    cued_DLPFC_gamma(subject, 2) = mean(dualcued_DLPFC(length(4:0.25:32):length(4:0.25:48)));
    
    % SMA
    singlecued_SMA = sub_struct.singlecued_power_SMA;
    
    % Theta band - 4:0.25:8
    cued_SMA_theta(subject, 1) = mean(singlecued_SMA(1:length(4:0.25:8)));
    cued_SMA_theta(subject, 2) = mean(dualcued_SMA(1:length(4:0.25:8)));
    % Alpha band - 8:0.25:13
    cued_SMA_alpha(subject, 1) = mean(singlecued_SMA(length(4:0.25:8):length(4:0.25:13)));
    cued_SMA_alpha(subject, 2) = mean(dualcued_SMA(length(4:0.25:8):length(4:0.25:13)));
    % Beta band - 13:0.25:32
    cued_SMA_beta(subject, 1) = mean(singlecued_SMA(length(4:0.25:13):length(4:0.25:32)));
    cued_SMA_beta(subject, 2) = mean(dualcued_SMA(length(4:0.25:13):length(4:0.25:32)));
    % Gamma band - 32:0.25:48
    cued_SMA_gamma(subject, 1) = mean(singlecued_SMA(length(4:0.25:32):length(4:0.25:48)));
    cued_SMA_gamma(subject, 2) = mean(dualcued_SMA(length(4:0.25:32):length(4:0.25:48)));
    
    % M1
    singlecued_M1 = sub_struct.singlecued_power_M1;
    
    % Theta band - 4:0.25:8
    cued_M1_theta(subject, 1) = mean(singlecued_M1(1:length(4:0.25:8)));
    cued_M1_theta(subject, 2) = mean(dualcued_M1(1:length(4:0.25:8)));
    % Alpha band - 8:0.25:13
    cued_M1_alpha(subject, 1) = mean(singlecued_M1(length(4:0.25:8):length(4:0.25:13)));
    cued_M1_alpha(subject, 2) = mean(dualcued_M1(length(4:0.25:8):length(4:0.25:13)));
    % Beta band - 13:0.25:32
    cued_M1_beta(subject, 1) = mean(singlecued_M1(length(4:0.25:13):length(4:0.25:32)));
    cued_M1_beta(subject, 2) = mean(dualcued_M1(length(4:0.25:13):length(4:0.25:32)));
    % Gamma band - 32:0.25:48
    cued_M1_gamma(subject, 1) = mean(singlecued_M1(length(4:0.25:32):length(4:0.25:48)));
    cued_M1_gamma(subject, 2) = mean(dualcued_M1(length(4:0.25:32):length(4:0.25:48)));
    
    % PPC
    singlecued_PPC = sub_struct.singlecued_power_PPC;
    
    % Theta band - 4:0.25:8
    cued_PPC_theta(subject, 1) = mean(singlecued_PPC(1:length(4:0.25:8)));
    cued_PPC_theta(subject, 2) = mean(dualcued_PPC(1:length(4:0.25:8)));
    % Alpha band - 8:0.25:13
    cued_PPC_alpha(subject, 1) = mean(singlecued_PPC(length(4:0.25:8):length(4:0.25:13)));
    cued_PPC_alpha(subject, 2) = mean(dualcued_PPC(length(4:0.25:8):length(4:0.25:13)));
    % Beta band - 13:0.25:32
    cued_PPC_beta(subject, 1) = mean(singlecued_PPC(length(4:0.25:13):length(4:0.25:32)));
    cued_PPC_beta(subject, 2) = mean(dualcued_PPC(length(4:0.25:13):length(4:0.25:32)));
    % Gamma band - 32:0.25:48
    cued_PPC_gamma(subject, 1) = mean(singlecued_PPC(length(4:0.25:32):length(4:0.25:48)));
    cued_PPC_gamma(subject, 2) = mean(dualcued_PPC(length(4:0.25:32):length(4:0.25:48)));
    
end

%% Wilcoxon signed-rank test - dual
% DLPFC
p = signrank(dual_DLPFC_theta(:, 1), dual_DLPFC_theta(:, 2));
stats_dual.dual_DLPFC_theta = p;
p = signrank(dual_DLPFC_alpha(:, 1), dual_DLPFC_alpha(:, 2));
stats_dual.dual_DLPFC_alpha = p;
p = signrank(dual_DLPFC_beta(:, 1), dual_DLPFC_beta(:, 2));
stats_dual.dual_DLPFC_beta = p;
p = signrank(dual_DLPFC_gamma(:, 1), dual_DLPFC_gamma(:, 2));
stats_dual.dual_DLPFC_gamma = p;

% SMA
p = signrank(dual_SMA_theta(:, 1), dual_SMA_theta(:, 2));
stats_dual.dual_SMA_theta = p;
p = signrank(dual_SMA_alpha(:, 1), dual_SMA_alpha(:, 2));
stats_dual.dual_SMA_alpha = p;
p = signrank(dual_SMA_beta(:, 1), dual_SMA_beta(:, 2));
stats_dual.dual_SMA_beta = p;
p = signrank(dual_SMA_gamma(:, 1), dual_SMA_gamma(:, 2));
stats_dual.dual_SMA_gamma = p;

% M1
p = signrank(dual_M1_theta(:, 1), dual_M1_theta(:, 2));
stats_dual.dual_M1_theta = p;
p = signrank(dual_M1_alpha(:, 1), dual_M1_alpha(:, 2));
stats_dual.dual_M1_alpha = p;
p = signrank(dual_M1_beta(:, 1), dual_M1_beta(:, 2));
stats_dual.dual_M1_beta = p;
p = signrank(dual_M1_gamma(:, 1), dual_M1_gamma(:, 2));
stats_dual.dual_M1_gamma = p;

% PPC
p = signrank(dual_PPC_theta(:, 1), dual_PPC_theta(:, 2));
stats_dual.dual_PPC_theta = p;
p = signrank(dual_PPC_alpha(:, 1), dual_PPC_alpha(:, 2));
stats_dual.dual_PPC_alpha = p;
p = signrank(dual_PPC_beta(:, 1), dual_PPC_beta(:, 2));
stats_dual.dual_PPC_beta = p;
p = signrank(dual_PPC_gamma(:, 1), dual_PPC_gamma(:, 2));
stats_dual.dual_PPC_gamma = p;

stats.stats_dual = stats_dual;

%% Wilcoxon signed-rank test - cued
% DLPFC
p = signrank(uncued_DLPFC_theta(:, 1), uncued_DLPFC_theta(:, 2));
stats_uncued.uncued_DLPFC_theta = p;
p = signrank(uncued_DLPFC_alpha(:, 1), uncued_DLPFC_alpha(:, 2));
stats_uncued.uncued_DLPFC_alpha = p;
p = signrank(uncued_DLPFC_beta(:, 1), uncued_DLPFC_beta(:, 2));
stats_uncued.uncued_DLPFC_beta = p;
p = signrank(uncued_DLPFC_gamma(:, 1), uncued_DLPFC_gamma(:, 2));
stats_uncued.uncued_DLPFC_gamma = p;

% SMA
p = signrank(uncued_SMA_theta(:, 1), uncued_SMA_theta(:, 2));
stats_uncued.uncued_SMA_theta = p;
p = signrank(uncued_SMA_alpha(:, 1), uncued_SMA_alpha(:, 2));
stats_uncued.uncued_SMA_alpha = p;
p = signrank(uncued_SMA_beta(:, 1), uncued_SMA_beta(:, 2));
stats_uncued.uncued_SMA_beta = p;
p = signrank(uncued_SMA_gamma(:, 1), uncued_SMA_gamma(:, 2));
stats_uncued.uncued_SMA_gamma = p;

% M1
p = signrank(uncued_M1_theta(:, 1), uncued_M1_theta(:, 2));
stats_uncued.uncued_M1_theta = p;
p = signrank(uncued_M1_alpha(:, 1), uncued_M1_alpha(:, 2));
stats_uncued.uncued_M1_alpha = p;
p = signrank(uncued_M1_beta(:, 1), uncued_M1_beta(:, 2));
stats_uncued.uncued_M1_beta = p;
p = signrank(uncued_M1_gamma(:, 1), uncued_M1_gamma(:, 2));
stats_uncued.uncued_M1_gamma = p;

% PPC
p = signrank(uncued_PPC_theta(:, 1), uncued_PPC_theta(:, 2));
stats_uncued.uncued_PPC_theta = p;
p = signrank(uncued_PPC_alpha(:, 1), uncued_PPC_alpha(:, 2));
stats_uncued.uncued_PPC_alpha = p;
p = signrank(uncued_PPC_beta(:, 1), uncued_PPC_beta(:, 2));
stats_uncued.uncued_PPC_beta = p;
p = signrank(uncued_PPC_gamma(:, 1), uncued_PPC_gamma(:, 2));
stats_uncued.uncued_PPC_gamma = p;

stats.stats_uncued = stats_uncued;

%% Wilcoxon signed-rank test - cued
% DLPFC
p = signrank(cued_DLPFC_theta(:, 1), cued_DLPFC_theta(:, 2));
stats_cued.cued_DLPFC_theta = p;
p = signrank(cued_DLPFC_alpha(:, 1), cued_DLPFC_alpha(:, 2));
stats_cued.cued_DLPFC_alpha = p;
p = signrank(cued_DLPFC_beta(:, 1), cued_DLPFC_beta(:, 2));
stats_cued.cued_DLPFC_beta = p;
p = signrank(cued_DLPFC_gamma(:, 1), cued_DLPFC_gamma(:, 2));
stats_cued.cued_DLPFC_gamma = p;

% SMA
p = signrank(cued_SMA_theta(:, 1), cued_SMA_theta(:, 2));
stats_cued.cued_SMA_theta = p;
p = signrank(cued_SMA_alpha(:, 1), cued_SMA_alpha(:, 2));
stats_cued.cued_SMA_alpha = p;
p = signrank(cued_SMA_beta(:, 1), cued_SMA_beta(:, 2));
stats_cued.cued_SMA_beta = p;
p = signrank(cued_SMA_gamma(:, 1), cued_SMA_gamma(:, 2));
stats_cued.cued_SMA_gamma = p;

% M1
p = signrank(cued_M1_theta(:, 1), cued_M1_theta(:, 2));
stats_cued.cued_M1_theta = p;
p = signrank(cued_M1_alpha(:, 1), cued_M1_alpha(:, 2));
stats_cued.cued_M1_alpha = p;
p = signrank(cued_M1_beta(:, 1), cued_M1_beta(:, 2));
stats_cued.cued_M1_beta = p;
p = signrank(cued_M1_gamma(:, 1), cued_M1_gamma(:, 2));
stats_cued.cued_M1_gamma = p;

% PPC
p = signrank(cued_PPC_theta(:, 1), cued_PPC_theta(:, 2));
stats_cued.cued_PPC_theta = p;
p = signrank(cued_PPC_alpha(:, 1), cued_PPC_alpha(:, 2));
stats_cued.cued_PPC_alpha = p;
p = signrank(cued_PPC_beta(:, 1), cued_PPC_beta(:, 2));
stats_cued.cued_PPC_beta = p;
p = signrank(cued_PPC_gamma(:, 1), cued_PPC_gamma(:, 2));
stats_cued.cued_PPC_gamma = p;

stats.stats_cued = stats_cued;

disp('Wilcoxon Signed Rank Test:');
disp(stats.stats_dual);
disp(stats.stats_uncued);
disp(stats.stats_cued);

save(strcat(statistics_path, '\stats_psd.mat'), 'stats');

