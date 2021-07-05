%% Statistical Analysis: EEG (ERD/ERS - Topoplots)
clear; clc; close all;

%% Initialize data
addpath('C:\Users\catar\OneDrive - Universidade do Porto\Twente\Combined-EEG-fNIRS-system\Analysis');
addpath('C:\Users\catar\OneDrive - Universidade do Porto\Twente\Combined-EEG-fNIRS-system\Analysis\Dual Task');
laptop = 'laptopCatarina';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);
results_path = 'C:\Users\catar\OneDrive - Universidade do Porto\Twente\Data Analysis\eeg\erders';
statistics_path = 'C:\Users\catar\OneDrive - Universidade do Porto\Twente\Data Analysis\statistics\eeg\erders (topoplots)';

load(fullfile(results_path, 'erders_allsubs.mat'), 'allsubs');

subrec = ["28" "04"; "02" "02"; "76" "01"; "64" "01"];

for subject=1:size(subrec, 1)
    sub = subrec(subject, 1);
    rec = subrec(subject, 2);
    
    % Subject 1 had no channels eliminated
    if subject==1
        load([mainpath_in, '\pre-processed\sub-', char(sub), '\eeg\sub-',...
            char(sub), '_rec-', char(rec), '_eeg_divided.mat']);
    end
    
    % Get the locations of the channels of interest
    locs = {EEG_divided.EEG_task.chanlocs.labels};
    % DLPFC
    F7_loc = find(contains(locs, 'F7'));
    F8_loc = find(contains(locs, 'F8'));
    % SMA
    FC1_loc = find(contains(locs, 'FC1'));
    FC2_loc = find(contains(locs, 'FC2'));
    % M1
    C3_loc = find(contains(locs, 'C3'));
    CP1_loc = find(contains(locs, 'CP1'));
    % PPC
    P3_loc = find(contains(locs, 'P3'));
    P4_loc = find(contains(locs, 'P4'));

    sub_struct = allsubs.(genvarname(strcat('sub', char(sub))));
    
    %% Get average value over regions for each frequency band - dual
    % Theta
    dualuncued_theta = sub_struct.dualuncued_ERD_ERS_theta;
    dualcued_theta = sub_struct.dualcued_ERD_ERS_theta;
    
    dual_theta_DLPFC(subject,1) = mean([dualuncued_theta(F7_loc)...
        dualuncued_theta(F8_loc)], 'omitnan');
    dual_theta_SMA(subject,1) = mean([dualuncued_theta(FC1_loc)...
        dualuncued_theta(FC2_loc)], 'omitnan');
    dual_theta_M1(subject,1) = mean([dualuncued_theta(C3_loc)...
        dualuncued_theta(CP1_loc)], 'omitnan');
    dual_theta_PPC(subject,1) = mean([dualuncued_theta(P3_loc)...
        dualuncued_theta(P4_loc)], 'omitnan');
    
    dual_theta_DLPFC(subject,2) = mean([dualcued_theta(F7_loc)...
        dualcued_theta(F8_loc)], 'omitnan');
    dual_theta_SMA(subject,2) = mean([dualcued_theta(FC1_loc)...
        dualcued_theta(FC2_loc)], 'omitnan');
    dual_theta_M1(subject,2) = mean([dualcued_theta(C3_loc)...
        dualcued_theta(CP1_loc)], 'omitnan');
    dual_theta_PPC(subject,2) = mean([dualcued_theta(P3_loc)...
        dualcued_theta(P4_loc)], 'omitnan');
   
    % Alpha
    dualuncued_alpha = sub_struct.dualuncued_ERD_ERS_alpha;
    dualcued_alpha = sub_struct.dualcued_ERD_ERS_alpha;
    
    dual_alpha_DLPFC(subject,1) = mean([dualuncued_alpha(F7_loc)...
        dualuncued_alpha(F8_loc)], 'omitnan');
    dual_alpha_SMA(subject,1) = mean([dualuncued_alpha(FC1_loc)...
        dualuncued_alpha(FC2_loc)], 'omitnan');
    dual_alpha_M1(subject,1) = mean([dualuncued_alpha(C3_loc)...
        dualuncued_alpha(CP1_loc)], 'omitnan');
    dual_alpha_PPC(subject,1) = mean([dualuncued_alpha(P3_loc)...
        dualuncued_alpha(P4_loc)], 'omitnan');
    
    dual_alpha_DLPFC(subject,2) = mean([dualcued_alpha(F7_loc)...
        dualcued_alpha(F8_loc)], 'omitnan');
    dual_alpha_SMA(subject,2) = mean([dualcued_alpha(FC1_loc)...
        dualcued_alpha(FC2_loc)], 'omitnan');
    dual_alpha_M1(subject,2) = mean([dualcued_alpha(C3_loc)...
        dualcued_alpha(CP1_loc)], 'omitnan');
    dual_alpha_PPC(subject,2) = mean([dualcued_alpha(P3_loc)...
        dualcued_alpha(P4_loc)], 'omitnan');
    
    % Beta
    dualuncued_beta = sub_struct.dualuncued_ERD_ERS_beta;
    dualcued_beta = sub_struct.dualcued_ERD_ERS_beta;
    
    dual_beta_DLPFC(subject,1) = mean([dualuncued_beta(F7_loc)...
        dualuncued_beta(F8_loc)], 'omitnan');
    dual_beta_SMA(subject,1) = mean([dualuncued_beta(FC1_loc)...
        dualuncued_beta(FC2_loc)], 'omitnan');
    dual_beta_M1(subject,1) = mean([dualuncued_beta(C3_loc)...
        dualuncued_beta(CP1_loc)], 'omitnan');
    dual_beta_PPC(subject,1) = mean([dualuncued_beta(P3_loc)...
        dualuncued_beta(P4_loc)], 'omitnan');
    
    dual_beta_DLPFC(subject,2) = mean([dualcued_beta(F7_loc)...
        dualcued_beta(F8_loc)], 'omitnan');
    dual_beta_SMA(subject,2) = mean([dualcued_beta(FC1_loc)...
        dualcued_beta(FC2_loc)], 'omitnan');
    dual_beta_M1(subject,2) = mean([dualcued_beta(C3_loc)...
        dualcued_beta(CP1_loc)], 'omitnan');
    dual_beta_PPC(subject,2) = mean([dualcued_beta(P3_loc)...
        dualcued_beta(P4_loc)], 'omitnan');
    
    % Gamma
    dualuncued_gamma = sub_struct.dualuncued_ERD_ERS_gamma;
    dualcued_gamma = sub_struct.dualcued_ERD_ERS_gamma;
    
    dual_gamma_DLPFC(subject,1) = mean([dualuncued_gamma(F7_loc)...
        dualuncued_gamma(F8_loc)], 'omitnan');
    dual_gamma_SMA(subject,1) = mean([dualuncued_gamma(FC1_loc)...
        dualuncued_gamma(FC2_loc)], 'omitnan');
    dual_gamma_M1(subject,1) = mean([dualuncued_gamma(C3_loc)...
        dualuncued_gamma(CP1_loc)], 'omitnan');
    dual_gamma_PPC(subject,1) = mean([dualuncued_gamma(P3_loc)...
        dualuncued_gamma(P4_loc)], 'omitnan');
    
    dual_gamma_DLPFC(subject,2) = mean([dualcued_gamma(F7_loc)...
        dualcued_gamma(F8_loc)], 'omitnan');
    dual_gamma_SMA(subject,2) = mean([dualcued_gamma(FC1_loc)...
        dualcued_gamma(FC2_loc)], 'omitnan');
    dual_gamma_M1(subject,2) = mean([dualcued_gamma(C3_loc)...
        dualcued_gamma(CP1_loc)], 'omitnan');
    dual_gamma_PPC(subject,2) = mean([dualcued_gamma(P3_loc)...
        dualcued_gamma(P4_loc)], 'omitnan');
    
    %% Get average value over regions for each frequency band - uncued
    % Theta
    singleuncued_theta = sub_struct.singleuncued_ERD_ERS_theta;
    
    uncued_theta_DLPFC(subject,1) = mean([singleuncued_theta(F7_loc)...
        singleuncued_theta(F8_loc)], 'omitnan');
    uncued_theta_SMA(subject,1) = mean([singleuncued_theta(FC1_loc)...
        singleuncued_theta(FC2_loc)], 'omitnan');
    uncued_theta_M1(subject,1) = mean([singleuncued_theta(C3_loc)...
        singleuncued_theta(CP1_loc)], 'omitnan');
    uncued_theta_PPC(subject,1) = mean([singleuncued_theta(P3_loc)...
        singleuncued_theta(P4_loc)], 'omitnan');
   
    uncued_theta_DLPFC(subject,2) = mean([dualuncued_theta(F7_loc)...
        dualuncued_theta(F8_loc)], 'omitnan');
    uncued_theta_SMA(subject,2) = mean([dualuncued_theta(FC1_loc)...
        dualuncued_theta(FC2_loc)], 'omitnan');
    uncued_theta_M1(subject,2) = mean([dualuncued_theta(C3_loc)...
        dualuncued_theta(CP1_loc)], 'omitnan');
    uncued_theta_PPC(subject,2) = mean([dualuncued_theta(P3_loc)...
        dualuncued_theta(P4_loc)], 'omitnan');
    
    % Alpha
    singleuncued_alpha = sub_struct.singleuncued_ERD_ERS_alpha;
    
    uncued_alpha_DLPFC(subject,1) = mean([singleuncued_alpha(F7_loc)...
        singleuncued_alpha(F8_loc)], 'omitnan');
    uncued_alpha_SMA(subject,1) = mean([singleuncued_alpha(FC1_loc)...
        singleuncued_alpha(FC2_loc)], 'omitnan');
    uncued_alpha_M1(subject,1) = mean([singleuncued_alpha(C3_loc)...
        singleuncued_alpha(CP1_loc)], 'omitnan');
    uncued_alpha_PPC(subject,1) = mean([singleuncued_alpha(P3_loc)...
        singleuncued_alpha(P4_loc)], 'omitnan');
    
    uncued_alpha_DLPFC(subject,2) = mean([dualuncued_alpha(F7_loc)...
        dualuncued_alpha(F8_loc)], 'omitnan');
    uncued_alpha_SMA(subject,2) = mean([dualuncued_alpha(FC1_loc)...
        dualuncued_alpha(FC2_loc)], 'omitnan');
    uncued_alpha_M1(subject,2) = mean([dualuncued_alpha(C3_loc)...
        dualuncued_alpha(CP1_loc)], 'omitnan');
    uncued_alpha_PPC(subject,2) = mean([dualuncued_alpha(P3_loc)...
        dualuncued_alpha(P4_loc)], 'omitnan');
    
    % Beta
    singleuncued_beta = sub_struct.singleuncued_ERD_ERS_beta;
    
    uncued_beta_DLPFC(subject,1) = mean([singleuncued_beta(F7_loc)...
        singleuncued_beta(F8_loc)], 'omitnan');
    uncued_beta_SMA(subject,1) = mean([singleuncued_beta(FC1_loc)...
        singleuncued_beta(FC2_loc)], 'omitnan');
    uncued_beta_M1(subject,1) = mean([singleuncued_beta(C3_loc)...
        singleuncued_beta(CP1_loc)], 'omitnan');
    uncued_beta_PPC(subject,1) = mean([singleuncued_beta(P3_loc)...
        singleuncued_beta(P4_loc)], 'omitnan');
    
    uncued_beta_DLPFC(subject,2) = mean([dualuncued_beta(F7_loc)...
        dualuncued_beta(F8_loc)], 'omitnan');
    uncued_beta_SMA(subject,2) = mean([dualuncued_beta(FC1_loc)...
        dualuncued_beta(FC2_loc)], 'omitnan');
    uncued_beta_M1(subject,2) = mean([dualuncued_beta(C3_loc)...
        dualuncued_beta(CP1_loc)], 'omitnan');
    uncued_beta_PPC(subject,2) = mean([dualuncued_beta(P3_loc)...
        dualuncued_beta(P4_loc)], 'omitnan');
    
    % Gamma
    singleuncued_gamma = sub_struct.singleuncued_ERD_ERS_gamma;
    
    uncued_gamma_DLPFC(subject,1) = mean([singleuncued_gamma(F7_loc)...
        singleuncued_gamma(F8_loc)], 'omitnan');
    uncued_gamma_SMA(subject,1) = mean([singleuncued_gamma(FC1_loc)...
        singleuncued_gamma(FC2_loc)], 'omitnan');
    uncued_gamma_M1(subject,1) = mean([singleuncued_gamma(C3_loc)...
        singleuncued_gamma(CP1_loc)], 'omitnan');
    uncued_gamma_PPC(subject,1) = mean([singleuncued_gamma(P3_loc)...
        singleuncued_gamma(P4_loc)], 'omitnan'); 
    
    uncued_gamma_DLPFC(subject,2) = mean([dualuncued_gamma(F7_loc)...
        dualuncued_gamma(F8_loc)], 'omitnan');
    uncued_gamma_SMA(subject,2) = mean([dualuncued_gamma(FC1_loc)...
        dualuncued_gamma(FC2_loc)], 'omitnan');
    uncued_gamma_M1(subject,2) = mean([dualuncued_gamma(C3_loc)...
        dualuncued_gamma(CP1_loc)], 'omitnan');
    uncued_gamma_PPC(subject,2) = mean([dualuncued_gamma(P3_loc)...
        dualuncued_gamma(P4_loc)], 'omitnan');
    
    %% Get average value over regions for each frequency band - cued
    % Theta
    singlecued_theta = sub_struct.singlecued_ERD_ERS_theta;
    
    cued_theta_DLPFC(subject,1) = mean([singlecued_theta(F7_loc)...
        singlecued_theta(F8_loc)], 'omitnan');
    cued_theta_SMA(subject,1) = mean([singlecued_theta(FC1_loc)...
        singlecued_theta(FC2_loc)], 'omitnan');
    cued_theta_M1(subject,1) = mean([singlecued_theta(C3_loc)...
        singlecued_theta(CP1_loc)], 'omitnan');
    cued_theta_PPC(subject,1) = mean([singlecued_theta(P3_loc)...
        singlecued_theta(P4_loc)], 'omitnan');
    
    cued_theta_DLPFC(subject,2) = mean([dualcued_theta(F7_loc)...
        dualcued_theta(F8_loc)], 'omitnan');
    cued_theta_SMA(subject,2) = mean([dualcued_theta(FC1_loc)...
        dualcued_theta(FC2_loc)], 'omitnan');
    cued_theta_M1(subject,2) = mean([dualcued_theta(C3_loc)...
        dualcued_theta(CP1_loc)], 'omitnan');
    cued_theta_PPC(subject,2) = mean([dualcued_theta(P3_loc)...
        dualcued_theta(P4_loc)], 'omitnan');
   
    % Alpha
    singlecued_alpha = sub_struct.singlecued_ERD_ERS_alpha;
    
    cued_alpha_DLPFC(subject,1) = mean([singlecued_alpha(F7_loc)...
        singlecued_alpha(F8_loc)], 'omitnan');
    cued_alpha_SMA(subject,1) = mean([singlecued_alpha(FC1_loc)...
        singlecued_alpha(FC2_loc)], 'omitnan');
    cued_alpha_M1(subject,1) = mean([singlecued_alpha(C3_loc)...
        singlecued_alpha(CP1_loc)], 'omitnan');
    cued_alpha_PPC(subject,1) = mean([singlecued_alpha(P3_loc)...
        singlecued_alpha(P4_loc)], 'omitnan');
    
    cued_alpha_DLPFC(subject,2) = mean([dualcued_alpha(F7_loc)...
        dualcued_alpha(F8_loc)], 'omitnan');
    cued_alpha_SMA(subject,2) = mean([dualcued_alpha(FC1_loc)...
        dualcued_alpha(FC2_loc)], 'omitnan');
    cued_alpha_M1(subject,2) = mean([dualcued_alpha(C3_loc)...
        dualcued_alpha(CP1_loc)], 'omitnan');
    cued_alpha_PPC(subject,2) = mean([dualcued_alpha(P3_loc)...
        dualcued_alpha(P4_loc)], 'omitnan');
    
    % Beta
    singlecued_beta = sub_struct.singlecued_ERD_ERS_beta;
    
    cued_beta_DLPFC(subject,1) = mean([singlecued_beta(F7_loc)...
        singlecued_beta(F8_loc)], 'omitnan');
    cued_beta_SMA(subject,1) = mean([singlecued_beta(FC1_loc)...
        singlecued_beta(FC2_loc)], 'omitnan');
    cued_beta_M1(subject,1) = mean([singlecued_beta(C3_loc)...
        singlecued_beta(CP1_loc)], 'omitnan');
    cued_beta_PPC(subject,1) = mean([singlecued_beta(P3_loc)...
        singlecued_beta(P4_loc)], 'omitnan');
    
    cued_beta_DLPFC(subject,2) = mean([dualcued_beta(F7_loc)...
        dualcued_beta(F8_loc)], 'omitnan');
    cued_beta_SMA(subject,2) = mean([dualcued_beta(FC1_loc)...
        dualcued_beta(FC2_loc)], 'omitnan');
    cued_beta_M1(subject,2) = mean([dualcued_beta(C3_loc)...
        dualcued_beta(CP1_loc)], 'omitnan');
    cued_beta_PPC(subject,2) = mean([dualcued_beta(P3_loc)...
        dualcued_beta(P4_loc)], 'omitnan');
    
    % Gamma
    singlecued_gamma = sub_struct.singlecued_ERD_ERS_gamma;
    
    cued_gamma_DLPFC(subject,1) = mean([singlecued_gamma(F7_loc)...
        singlecued_gamma(F8_loc)], 'omitnan');
    cued_gamma_SMA(subject,1) = mean([singlecued_gamma(FC1_loc)...
        singlecued_gamma(FC2_loc)], 'omitnan');
    cued_gamma_M1(subject,1) = mean([singlecued_gamma(C3_loc)...
        singlecued_gamma(CP1_loc)], 'omitnan');
    cued_gamma_PPC(subject,1) = mean([singlecued_gamma(P3_loc)...
        singlecued_gamma(P4_loc)], 'omitnan');
    
    cued_gamma_DLPFC(subject,2) = mean([dualcued_gamma(F7_loc)...
        dualcued_gamma(F8_loc)], 'omitnan');
    cued_gamma_SMA(subject,2) = mean([dualcued_gamma(FC1_loc)...
        dualcued_gamma(FC2_loc)], 'omitnan');
    cued_gamma_M1(subject,2) = mean([dualcued_gamma(C3_loc)...
        dualcued_gamma(CP1_loc)], 'omitnan');
    cued_gamma_PPC(subject,2) = mean([dualcued_gamma(P3_loc)...
        dualcued_gamma(P4_loc)], 'omitnan');
    
end

%% Wilcoxon signed-rank test - dual 
% Theta band
p = signrank(dual_theta_DLPFC(:, 1), dual_theta_DLPFC(:, 2));
stats_dual.dual_theta_DLPFC = p;
p = signrank(dual_theta_SMA(:, 1), dual_theta_SMA(:, 2));
stats_dual.dual_theta_SMA = p;
p = signrank(dual_theta_M1(:, 1), dual_theta_M1(:, 2));
stats_dual.dual_theta_M1 = p;
p = signrank(dual_theta_PPC(:, 1), dual_theta_PPC(:, 2));
stats_dual.dual_theta_PPC = p;

% Alpha band
p = signrank(dual_alpha_DLPFC(:, 1), dual_alpha_DLPFC(:, 2));
stats_dual.dual_alpha_DLPFC = p;
p = signrank(dual_alpha_SMA(:, 1), dual_alpha_SMA(:, 2));
stats_dual.dual_alpha_SMA = p;
p = signrank(dual_alpha_M1(:, 1), dual_alpha_M1(:, 2));
stats_dual.dual_alpha_M1 = p;
p = signrank(dual_alpha_PPC(:, 1), dual_alpha_PPC(:, 2));
stats_dual.dual_alpha_PPC = p;

% Beta band
p = signrank(dual_beta_DLPFC(:, 1), dual_beta_DLPFC(:, 2));
stats_dual.dual_beta_DLPFC = p;
p = signrank(dual_beta_SMA(:, 1), dual_beta_SMA(:, 2));
stats_dual.dual_beta_SMA = p;
p = signrank(dual_beta_M1(:, 1), dual_beta_M1(:, 2));
stats_dual.dual_beta_M1 = p;
p = signrank(dual_beta_PPC(:, 1), dual_beta_PPC(:, 2));
stats_dual.dual_beta_PPC = p;

% Gamma band
p = signrank(dual_gamma_DLPFC(:, 1), dual_gamma_DLPFC(:, 2));
stats_dual.dual_gamma_DLPFC = p;
p = signrank(dual_gamma_SMA(:, 1), dual_gamma_SMA(:, 2));
stats_dual.dual_gamma_SMA = p;
p = signrank(dual_gamma_M1(:, 1), dual_gamma_M1(:, 2));
stats_dual.dual_gamma_M1 = p;
p = signrank(dual_gamma_PPC(:, 1), dual_gamma_PPC(:, 2));
stats_dual.dual_gamma_PPC = p;

stats.stats_dual = stats_dual;

%% Wilcoxon signed-rank test - uncued
% Theta band
p = signrank(uncued_theta_DLPFC(:, 1), uncued_theta_DLPFC(:, 2));
stats_uncued.uncued_theta_DLPFC = p;
p = signrank(uncued_theta_SMA(:, 1), uncued_theta_SMA(:, 2));
stats_uncued.uncued_theta_SMA = p;
p = signrank(uncued_theta_M1(:, 1), uncued_theta_M1(:, 2));
stats_uncued.uncued_theta_M1 = p;
p = signrank(uncued_theta_PPC(:, 1), uncued_theta_PPC(:, 2));
stats_uncued.uncued_theta_PPC = p;

% Alpha band
p = signrank(uncued_alpha_DLPFC(:, 1), uncued_alpha_DLPFC(:, 2));
stats_uncued.uncued_alpha_DLPFC = p;
p = signrank(uncued_alpha_SMA(:, 1), uncued_alpha_SMA(:, 2));
stats_uncued.uncued_alpha_SMA = p;
p = signrank(uncued_alpha_M1(:, 1), uncued_alpha_M1(:, 2));
stats_uncued.uncued_alpha_M1 = p;
p = signrank(uncued_alpha_PPC(:, 1), uncued_alpha_PPC(:, 2));
stats_uncued.uncued_alpha_PPC = p;

% Beta band
p = signrank(uncued_beta_DLPFC(:, 1), uncued_beta_DLPFC(:, 2));
stats_uncued.uncued_beta_DLPFC = p;
p = signrank(uncued_beta_SMA(:, 1), uncued_beta_SMA(:, 2));
stats_uncued.uncued_beta_SMA = p;
p = signrank(uncued_beta_M1(:, 1), uncued_beta_M1(:, 2));
stats_uncued.uncued_beta_M1 = p;
p = signrank(uncued_beta_PPC(:, 1), uncued_beta_PPC(:, 2));
stats_uncued.uncued_beta_PPC = p;

% Gamma band
p = signrank(uncued_gamma_DLPFC(:, 1), uncued_gamma_DLPFC(:, 2));
stats_uncued.uncued_gamma_DLPFC = p;
p = signrank(uncued_gamma_SMA(:, 1), uncued_gamma_SMA(:, 2));
stats_uncued.uncued_gamma_SMA = p;
p = signrank(uncued_gamma_M1(:, 1), uncued_gamma_M1(:, 2));
stats_uncued.uncued_gamma_M1 = p;
p = signrank(uncued_gamma_PPC(:, 1), uncued_gamma_PPC(:, 2));
stats_uncued.uncued_gamma_PPC = p;

stats.stats_uncued = stats_uncued;

%% Wilcoxon signed-rank test - cued
% Theta band
p = signrank(cued_theta_DLPFC(:, 1), cued_theta_DLPFC(:, 2));
stats_cued.cued_theta_DLPFC = p;
p = signrank(cued_theta_SMA(:, 1), cued_theta_SMA(:, 2));
stats_cued.cued_theta_SMA = p;
p = signrank(cued_theta_M1(:, 1), cued_theta_M1(:, 2));
stats_cued.cued_theta_M1 = p;
p = signrank(cued_theta_PPC(:, 1), cued_theta_PPC(:, 2));
stats_cued.cued_theta_PPC = p;

% Alpha band
p = signrank(cued_alpha_DLPFC(:, 1), cued_alpha_DLPFC(:, 2));
stats_cued.cued_alpha_DLPFC = p;
p = signrank(cued_alpha_SMA(:, 1), cued_alpha_SMA(:, 2));
stats_cued.cued_alpha_SMA = p;
p = signrank(cued_alpha_M1(:, 1), cued_alpha_M1(:, 2));
stats_cued.cued_alpha_M1 = p;
p = signrank(cued_alpha_PPC(:, 1), cued_alpha_PPC(:, 2));
stats_cued.cued_alpha_PPC = p;

% Beta band
p = signrank(cued_beta_DLPFC(:, 1), cued_beta_DLPFC(:, 2));
stats_cued.cued_beta_DLPFC = p;
p = signrank(cued_beta_SMA(:, 1), cued_beta_SMA(:, 2));
stats_cued.cued_beta_SMA = p;
p = signrank(cued_beta_M1(:, 1), cued_beta_M1(:, 2));
stats_cued.cued_beta_M1 = p;
p = signrank(cued_beta_PPC(:, 1), cued_beta_PPC(:, 2));
stats_cued.cued_beta_PPC = p;

% Gamma band
p = signrank(cued_gamma_DLPFC(:, 1), cued_gamma_DLPFC(:, 2));
stats_cued.cued_gamma_DLPFC = p;
p = signrank(cued_gamma_SMA(:, 1), cued_gamma_SMA(:, 2));
stats_cued.cued_gamma_SMA = p;
p = signrank(cued_gamma_M1(:, 1), cued_gamma_M1(:, 2));
stats_cued.cued_gamma_M1 = p;
p = signrank(cued_gamma_PPC(:, 1), cued_gamma_PPC(:, 2));
stats_cued.cued_gamma_PPC = p;

stats.stats_cued = stats_cued;

disp(stats.stats_dual);
disp(stats.stats_uncued);
disp(stats.stats_cued);

save(strcat(statistics_path, '\stats_erders.mat'), 'stats');
