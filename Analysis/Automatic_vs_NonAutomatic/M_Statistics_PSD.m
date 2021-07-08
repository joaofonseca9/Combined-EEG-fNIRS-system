%% Statistic for ERD/ERS.

clear; clc; close all;
addpath('C:\Users\maria\OneDrive\Documentos\GitHub\Combined-EEG-fNIRS-system\Analysis');

laptop = 'laptopMariana';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);
results_path = 'C:\Users\maria\OneDrive\Ambiente de Trabalho\Automaticity Results\Separate Channels';
statistics_path = 'C:\Users\maria\OneDrive\Ambiente de Trabalho\Statistics\Separate Channels';

load(fullfile(results_path, 'psd_allsubs.mat'), 'allsubs');

subrec = ["28" "04"; "02" "02"; "76" "01"];

for subject=1:size(subrec, 1)
    sub = subrec(subject, 1);
    rec = subrec(subject, 2);
    
    sub_struct = allsubs.(genvarname(strcat('sub', char(sub))));
    
    %% Get average value over frequency bands for each region - automatic
    % sequence.
    
    % DLPFC.
    autouncued_DLPFC = sub_struct.autouncued_power_DLPFC;
    autocued_DLPFC = sub_struct.autocued_power_DLPFC;
    
    % Theta band - 4:0.25:8.
    auto_DLPFC_theta(subject, 1) = mean(autouncued_DLPFC(1:length(4:0.25:8)));
    auto_DLPFC_theta(subject, 2) = mean(autocued_DLPFC(1:length(4:0.25:8)));
    % Alpha band - 8:0.25:13.
    auto_DLPFC_alpha(subject, 1) = mean(autouncued_DLPFC(length(4:0.25:8):length(4:0.25:13)));
    auto_DLPFC_alpha(subject, 2) = mean(autocued_DLPFC(length(4:0.25:8):length(4:0.25:13)));
    % Beta band - 13:0.25:32.
    auto_DLPFC_beta(subject, 1) = mean(autouncued_DLPFC(length(4:0.25:13):length(4:0.25:32)));
    auto_DLPFC_beta(subject, 2) = mean(autocued_DLPFC(length(4:0.25:13):length(4:0.25:32)));
    % Gamma band - 32:0.25:48.
    auto_DLPFC_gamma(subject, 1) = mean(autouncued_DLPFC(length(4:0.25:32):length(4:0.25:48)));
    auto_DLPFC_gamma(subject, 2) = mean(autocued_DLPFC(length(4:0.25:32):length(4:0.25:48)));
    
    % SMA.
    autouncued_SMA = sub_struct.autouncued_power_SMA;
    autocued_SMA = sub_struct.autocued_power_SMA;
    
    % Theta band - 4:0.25:8.
    auto_SMA_theta(subject, 1) = mean(autouncued_SMA(1:length(4:0.25:8)));
    auto_SMA_theta(subject, 2) = mean(autocued_SMA(1:length(4:0.25:8)));
    % Alpha band - 8:0.25:13.
    auto_SMA_alpha(subject, 1) = mean(autouncued_SMA(length(4:0.25:8):length(4:0.25:13)));
    auto_SMA_alpha(subject, 2) = mean(autocued_SMA(length(4:0.25:8):length(4:0.25:13)));
    % Beta band - 13:0.25:32.
    auto_SMA_beta(subject, 1) = mean(autouncued_SMA(length(4:0.25:13):length(4:0.25:32)));
    auto_SMA_beta(subject, 2) = mean(autocued_SMA(length(4:0.25:13):length(4:0.25:32)));
    % Gamma band - 32:0.25:48.
    auto_SMA_gamma(subject, 1) = mean(autouncued_SMA(length(4:0.25:32):length(4:0.25:48)));
    auto_SMA_gamma(subject, 2) = mean(autocued_SMA(length(4:0.25:32):length(4:0.25:48)));
    
    % M1.
    autouncued_M1 = sub_struct.autouncued_power_M1;
    autocued_M1 = sub_struct.autocued_power_M1;
    
    % Theta band - 4:0.25:8.
    auto_M1_theta(subject, 1) = mean(autouncued_M1(1:length(4:0.25:8)));
    auto_M1_theta(subject, 2) = mean(autocued_M1(1:length(4:0.25:8)));
    % Alpha band - 8:0.25:13.
    auto_M1_alpha(subject, 1) = mean(autouncued_M1(length(4:0.25:8):length(4:0.25:13)));
    auto_M1_alpha(subject, 2) = mean(autocued_M1(length(4:0.25:8):length(4:0.25:13)));
    % Beta band - 13:0.25:32.
    auto_M1_beta(subject, 1) = mean(autouncued_M1(length(4:0.25:13):length(4:0.25:32)));
    auto_M1_beta(subject, 2) = mean(autocued_M1(length(4:0.25:13):length(4:0.25:32)));
    % Gamma band - 32:0.25:48.
    auto_M1_gamma(subject, 1) = mean(autouncued_M1(length(4:0.25:32):length(4:0.25:48)));
    auto_M1_gamma(subject, 2) = mean(autocued_M1(length(4:0.25:32):length(4:0.25:48)));
    
    % PPC.
    autouncued_PPC = sub_struct.autouncued_power_PPC;
    autocued_PPC = sub_struct.autocued_power_PPC;
    
    % Theta band - 4:0.25:8.
    auto_PPC_theta(subject, 1) = mean(autouncued_PPC(1:length(4:0.25:8)));
    auto_PPC_theta(subject, 2) = mean(autocued_PPC(1:length(4:0.25:8)));
    % Alpha band - 8:0.25:13.
    auto_PPC_alpha(subject, 1) = mean(autouncued_PPC(length(4:0.25:8):length(4:0.25:13)));
    auto_PPC_alpha(subject, 2) = mean(autocued_PPC(length(4:0.25:8):length(4:0.25:13)));
    % Beta band - 13:0.25:32.
    auto_PPC_beta(subject, 1) = mean(autouncued_PPC(length(4:0.25:13):length(4:0.25:32)));
    auto_PPC_beta(subject, 2) = mean(autocued_PPC(length(4:0.25:13):length(4:0.25:32)));
    % Gamma band - 32:0.25:48.
    auto_PPC_gamma(subject, 1) = mean(autouncued_PPC(length(4:0.25:32):length(4:0.25:48)));
    auto_PPC_gamma(subject, 2) = mean(autocued_PPC(length(4:0.25:32):length(4:0.25:48)));
    
    %% Get average value over frequency bands for each region - non-automatic
    % sequence.
    
    % DLPFC.
    nonautouncued_DLPFC = sub_struct.nonautouncued_power_DLPFC;
    nonautocued_DLPFC = sub_struct.nonautocued_power_DLPFC;
    
    % Theta band - 4:0.25:8.
    nonauto_DLPFC_theta(subject, 1) = mean(nonautouncued_DLPFC(1:length(4:0.25:8)));
    nonauto_DLPFC_theta(subject, 2) = mean(nonautocued_DLPFC(1:length(4:0.25:8)));
    % Alpha band - 8:0.25:13.
    nonauto_DLPFC_alpha(subject, 1) = mean(nonautouncued_DLPFC(length(4:0.25:8):length(4:0.25:13)));
    nonauto_DLPFC_alpha(subject, 2) = mean(nonautocued_DLPFC(length(4:0.25:8):length(4:0.25:13)));
    % Beta band - 13:0.25:32.
    nonauto_DLPFC_beta(subject, 1) = mean(nonautouncued_DLPFC(length(4:0.25:13):length(4:0.25:32)));
    nonauto_DLPFC_beta(subject, 2) = mean(nonautocued_DLPFC(length(4:0.25:13):length(4:0.25:32)));
    % Gamma band - 32:0.25:48.
    nonauto_DLPFC_gamma(subject, 1) = mean(nonautouncued_DLPFC(length(4:0.25:32):length(4:0.25:48)));
    nonauto_DLPFC_gamma(subject, 2) = mean(nonautocued_DLPFC(length(4:0.25:32):length(4:0.25:48)));
    
    % SMA.
    nonautouncued_SMA = sub_struct.nonautouncued_power_SMA;
    nonautocued_SMA = sub_struct.nonautocued_power_SMA;
    
    % Theta band - 4:0.25:8.
    nonauto_SMA_theta(subject, 1) = mean(nonautouncued_SMA(1:length(4:0.25:8)));
    nonauto_SMA_theta(subject, 2) = mean(nonautocued_SMA(1:length(4:0.25:8)));
    % Alpha band - 8:0.25:13.
    nonauto_SMA_alpha(subject, 1) = mean(nonautouncued_SMA(length(4:0.25:8):length(4:0.25:13)));
    nonauto_SMA_alpha(subject, 2) = mean(nonautocued_SMA(length(4:0.25:8):length(4:0.25:13)));
    % Beta band - 13:0.25:32.
    nonauto_SMA_beta(subject, 1) = mean(nonautouncued_SMA(length(4:0.25:13):length(4:0.25:32)));
    nonauto_SMA_beta(subject, 2) = mean(nonautocued_SMA(length(4:0.25:13):length(4:0.25:32)));
    % Gamma band - 32:0.25:48.
    nonauto_SMA_gamma(subject, 1) = mean(nonautouncued_SMA(length(4:0.25:32):length(4:0.25:48)));
    nonauto_SMA_gamma(subject, 2) = mean(nonautocued_SMA(length(4:0.25:32):length(4:0.25:48)));
    
    % M1.
    nonautouncued_M1 = sub_struct.nonautouncued_power_M1;
    nonautocued_M1 = sub_struct.nonautocued_power_M1;
    
    % Theta band - 4:0.25:8.
    nonauto_M1_theta(subject, 1) = mean(nonautouncued_M1(1:length(4:0.25:8)));
    nonauto_M1_theta(subject, 2) = mean(nonautocued_M1(1:length(4:0.25:8)));
    % Alpha band - 8:0.25:13.
    nonauto_M1_alpha(subject, 1) = mean(nonautouncued_M1(length(4:0.25:8):length(4:0.25:13)));
    nonauto_M1_alpha(subject, 2) = mean(nonautocued_M1(length(4:0.25:8):length(4:0.25:13)));
    % Beta band - 13:0.25:32.
    nonauto_M1_beta(subject, 1) = mean(nonautouncued_M1(length(4:0.25:13):length(4:0.25:32)));
    nonauto_M1_beta(subject, 2) = mean(nonautocued_M1(length(4:0.25:13):length(4:0.25:32)));
    % Gamma band - 32:0.25:48.
    nonauto_M1_gamma(subject, 1) = mean(nonautouncued_M1(length(4:0.25:32):length(4:0.25:48)));
    nonauto_M1_gamma(subject, 2) = mean(nonautocued_M1(length(4:0.25:32):length(4:0.25:48)));
    
    % PPC.
    nonautouncued_PPC = sub_struct.nonautouncued_power_PPC;
    nonautocued_PPC = sub_struct.nonautocued_power_PPC;
    
    % Theta band - 4:0.25:8.
    nonauto_PPC_theta(subject, 1) = mean(nonautouncued_PPC(1:length(4:0.25:8)));
    nonauto_PPC_theta(subject, 2) = mean(nonautocued_PPC(1:length(4:0.25:8)));
    % Alpha band - 8:0.25:13.
    nonauto_PPC_alpha(subject, 1) = mean(nonautouncued_PPC(length(4:0.25:8):length(4:0.25:13)));
    nonauto_PPC_alpha(subject, 2) = mean(nonautocued_PPC(length(4:0.25:8):length(4:0.25:13)));
    % Beta band - 13:0.25:32.
    nonauto_PPC_beta(subject, 1) = mean(nonautouncued_PPC(length(4:0.25:13):length(4:0.25:32)));
    nonauto_PPC_beta(subject, 2) = mean(nonautocued_PPC(length(4:0.25:13):length(4:0.25:32)));
    % Gamma band - 32:0.25:48.
    nonauto_PPC_gamma(subject, 1) = mean(nonautouncued_PPC(length(4:0.25:32):length(4:0.25:48)));
    nonauto_PPC_gamma(subject, 2) = mean(nonautocued_PPC(length(4:0.25:32):length(4:0.25:48)));
    
end

%% Automatic sequence.
% Wilcoxon signed-rank test if not normally distributed.
% T-test if normally distributed.

% DLPFC.
p = signrank(auto_DLPFC_theta(:, 1), auto_DLPFC_theta(:, 2));
stats_auto.auto_DLPFC_theta = p;
p = signrank(auto_DLPFC_alpha(:, 1), auto_DLPFC_alpha(:, 2));
stats_auto.auto_DLPFC_alpha = p;
p = signrank(auto_DLPFC_beta(:, 1), auto_DLPFC_beta(:, 2));
stats_auto.auto_DLPFC_beta = p;
p = signrank(auto_DLPFC_gamma(:, 1), auto_DLPFC_gamma(:, 2));
stats_auto.auto_DLPFC_gamma = p;

% SMA.
p = signrank(auto_SMA_theta(:, 1), auto_SMA_theta(:, 2));
stats_auto.auto_SMA_theta = p;
p = signrank(auto_SMA_alpha(:, 1), auto_SMA_alpha(:, 2));
stats_auto.auto_SMA_alpha = p;
p = signrank(auto_SMA_beta(:, 1), auto_SMA_beta(:, 2));
stats_auto.auto_SMA_beta = p;
p = signrank(auto_SMA_gamma(:, 1), auto_SMA_gamma(:, 2));
stats_auto.auto_SMA_gamma = p;

% M1.
p = signrank(auto_M1_theta(:, 1), auto_M1_theta(:, 2));
stats_auto.auto_M1_theta = p;
p = signrank(auto_M1_alpha(:, 1), auto_M1_alpha(:, 2));
stats_auto.auto_M1_alpha = p;
p = signrank(auto_M1_beta(:, 1), auto_M1_beta(:, 2));
stats_auto.auto_M1_beta = p;
p = signrank(auto_M1_gamma(:, 1), auto_M1_gamma(:, 2));
stats_auto.auto_M1_gamma = p;

% PPC.
p = signrank(auto_PPC_theta(:, 1), auto_PPC_theta(:, 2));
stats_auto.auto_PPC_theta = p;
p = signrank(auto_PPC_alpha(:, 1), auto_PPC_alpha(:, 2));
stats_auto.auto_PPC_alpha = p;
p = signrank(auto_PPC_beta(:, 1), auto_PPC_beta(:, 2));
stats_auto.auto_PPC_beta = p;
p = signrank(auto_PPC_gamma(:, 1), auto_PPC_gamma(:, 2));
stats_auto.auto_PPC_gamma = p;

stats.stats_auto = stats_auto;

%% Non-automatic sequence.
% Wilcoxon signed-rank test if not normally distributed.
% T-test if normally distributed.

% DLPFC.
p = signrank(nonauto_DLPFC_theta(:, 1), nonauto_DLPFC_theta(:, 2));
stats_nonauto.nonauto_DLPFC_theta = p;
p = signrank(nonauto_DLPFC_alpha(:, 1), nonauto_DLPFC_alpha(:, 2));
stats_nonauto.nonauto_DLPFC_alpha = p;
p = signrank(nonauto_DLPFC_beta(:, 1), nonauto_DLPFC_beta(:, 2));
stats_nonauto.nonauto_DLPFC_beta = p;
p = signrank(nonauto_DLPFC_gamma(:, 1), nonauto_DLPFC_gamma(:, 2));
stats_nonauto.nonauto_DLPFC_gamma = p;

% SMA.
p = signrank(nonauto_SMA_theta(:, 1), nonauto_SMA_theta(:, 2));
stats_nonauto.nonauto_SMA_theta = p;
p = signrank(nonauto_SMA_alpha(:, 1), nonauto_SMA_alpha(:, 2));
stats_nonauto.nonauto_SMA_alpha = p;
p = signrank(nonauto_SMA_beta(:, 1), nonauto_SMA_beta(:, 2));
stats_nonauto.nonauto_SMA_beta = p;
p = signrank(nonauto_SMA_gamma(:, 1), nonauto_SMA_gamma(:, 2));
stats_nonauto.nonauto_SMA_gamma = p;

% M1.
p = signrank(nonauto_M1_theta(:, 1), nonauto_M1_theta(:, 2));
stats_nonauto.nonauto_M1_theta = p;
p = signrank(nonauto_M1_alpha(:, 1), nonauto_M1_alpha(:, 2));
stats_nonauto.nonauto_M1_alpha = p;
p = signrank(nonauto_M1_beta(:, 1), nonauto_M1_beta(:, 2));
stats_nonauto.nonauto_M1_beta = p;
p = signrank(nonauto_M1_gamma(:, 1), nonauto_M1_gamma(:, 2));
stats_nonauto.nonauto_M1_gamma = p;

% PPC.
p = signrank(nonauto_PPC_theta(:, 1), nonauto_PPC_theta(:, 2));
stats_nonauto.nonauto_PPC_theta = p;
p = signrank(nonauto_PPC_alpha(:, 1), nonauto_PPC_alpha(:, 2));
stats_nonauto.nonauto_PPC_alpha = p;
p = signrank(nonauto_PPC_beta(:, 1), nonauto_PPC_beta(:, 2));
stats_nonauto.nonauto_PPC_beta = p;
p = signrank(nonauto_PPC_gamma(:, 1), nonauto_PPC_gamma(:, 2));
stats_nonauto.nonauto_PPC_gamma = p;

stats.stats_nonauto = stats_nonauto;

save(strcat(statistics_path, '\stats_psd.mat'), 'stats');