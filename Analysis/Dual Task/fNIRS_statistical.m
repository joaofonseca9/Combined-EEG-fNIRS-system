%% Statistical Analysis: NIRS
clear; clc; close all;

%% Initialize data
addpath('C:\Users\catar\OneDrive - Universidade do Porto\Twente\Combined-EEG-fNIRS-system\Analysis');
addpath('C:\Users\catar\OneDrive - Universidade do Porto\Twente\Combined-EEG-fNIRS-system\Analysis\Dual Task');
laptop = 'laptopCatarina';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);

results_path = 'C:\Users\catar\OneDrive - Universidade do Porto\Twente\Data Analysis\nirs';
statistics_path = 'C:\Users\catar\OneDrive - Universidade do Porto\Twente\Data Analysis\statistics\nirs';

ft_defaults;
[~, ftpath] = ft_version;

subrec = ["28" "02";"64" "01";"02" "02";"76" "01"];
conditions = [3 4 7 8];
taskname = {'Dual Cued', 'Single Cued', 'Dual Uncued', 'Single Uncued'};
    
for subject = 1:size(subrec, 1)
    sub = subrec(subject, 1);
    rec = subrec(subject, 2);
    
    % Load the subject's file
    load(fullfile(results_path, ['sub-', char(sub)],...
            'nirs_TLblc.mat'), 'nirs_TLblc');
        
    % Separate Hb02 from Hb    
    [nirs_TLO2Hb, nirs_TLHHb] = separateHbO2FromHb(conditions, nirs_TLblc, subject);
    
    % Extract the different regions of interest: DLPFC, SMA, M1 and PPC
    % Only for HbO2 signal
    [nirs_HbO2_DLPFC, nirs_HbO2_SMA, nirs_HbO2_M1, nirs_HbO2_PPC] =...
        extractROIs(nirs_TLO2Hb);
    
    % Separate into the different conditions
    % DLPFC
    [nirs_dualcued, nirs_singlecued, nirs_dualuncued, nirs_singleuncued]...
        = extractConditions(conditions, nirs_HbO2_DLPFC);
    nirs_dualcued_DLPFC{subject} = nirs_dualcued;
    nirs_singlecued_DLPFC{subject} = nirs_singlecued;
    nirs_dualuncued_DLPFC{subject} = nirs_dualuncued;
    nirs_singleuncued_DLPFC{subject} = nirs_singleuncued;
    
    % SMA
    [nirs_dualcued, nirs_singlecued, nirs_dualuncued, nirs_singleuncued]...
        = extractConditions(conditions, nirs_HbO2_SMA);
    nirs_dualcued_SMA{subject} = nirs_dualcued;
    nirs_singlecued_SMA{subject} = nirs_singlecued;
    nirs_dualuncued_SMA{subject} = nirs_dualuncued;
    nirs_singleuncued_SMA{subject} = nirs_singleuncued;
    
    % M1
    [nirs_dualcued, nirs_singlecued, nirs_dualuncued, nirs_singleuncued]...
        = extractConditions(conditions, nirs_HbO2_M1);
    nirs_dualcued_M1{subject} = nirs_dualcued;
    nirs_singlecued_M1{subject} = nirs_singlecued;
    nirs_dualuncued_M1{subject} = nirs_dualuncued;
    nirs_singleuncued_M1{subject} = nirs_singleuncued;
    
    % PPC
    [nirs_dualcued, nirs_singlecued, nirs_dualuncued, nirs_singleuncued]...
        = extractConditions(conditions, nirs_HbO2_PPC);
    nirs_dualcued_PPC{subject} = nirs_dualcued;
    nirs_singlecued_PPC{subject} = nirs_singlecued;
    nirs_dualuncued_PPC{subject} = nirs_dualuncued;
    nirs_singleuncued_PPC{subject} = nirs_singleuncued;
    
end

%% T-TEST
% Define the parameters for the statistical comparison
cfg = [];
cfg.latency = 'all';
cfg.avgoverchan = 'yes';
cfg.avgovertime = 'yes';
cfg.parameter = 'avg';
cfg.method = 'analytic';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.alpha = 0.05;
cfg.correctm = 'no';

Nsub = 4;
cfg.design(1,1:2*Nsub) = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub) = [1:Nsub 1:Nsub];
% The 1st row in cfg.design contains the independent variable
cfg.ivar = 1; 
% The 2nd row in cfg.design contains the subject number
cfg.uvar = 2; 

% Dual Uncued vs Cued
% DLPFC
stat_dual_DLPFC = ft_timelockstatistics(cfg, nirs_dualuncued_DLPFC{:},...
    nirs_dualcued_DLPFC{:});
stats.dual_DLPFC = stat_dual_DLPFC.prob;
% SMA
stat_dual_SMA = ft_timelockstatistics(cfg, nirs_dualuncued_SMA{:},...
    nirs_dualcued_SMA{:});
stats.dual_SMA = stat_dual_SMA.prob;
% M1
stat_dual_M1 = ft_timelockstatistics(cfg, nirs_dualuncued_M1{:},...
    nirs_dualcued_M1{:});
stats.dual_M1 = stat_dual_M1.prob;
% PPC
stat_dual_PPC = ft_timelockstatistics(cfg, nirs_dualuncued_PPC{:},...
    nirs_dualcued_PPC{:});
stats.dual_PPC = stat_dual_PPC.prob;

% Single Uncued vs Dual Uncued
% DLPFC
stat_uncued_DLPFC = ft_timelockstatistics(cfg, nirs_singleuncued_DLPFC{:},...
    nirs_dualuncued_DLPFC{:});
stats.uncued_DLPFC = stat_uncued_DLPFC.prob;
% SMA
stat_uncued_SMA = ft_timelockstatistics(cfg, nirs_singleuncued_SMA{:},...
    nirs_dualuncued_SMA{:});
stats.uncued_SMA = stat_uncued_SMA.prob;
% M1
stat_uncued_M1 = ft_timelockstatistics(cfg, nirs_singleuncued_M1{:},...
    nirs_dualuncued_M1{:});
stats.uncued_M1 = stat_uncued_M1.prob;
% PPC
stat_uncued_PPC = ft_timelockstatistics(cfg, nirs_singleuncued_PPC{:},...
    nirs_dualuncued_PPC{:});
stats.uncued_PPC = stat_uncued_PPC.prob;

% Single Cued vs Dual Cued
% DLPFC
stat_cued_DLPFC = ft_timelockstatistics(cfg, nirs_singlecued_DLPFC{:},...
    nirs_dualcued_DLPFC{:});
stats.cued_DLPFC = stat_cued_DLPFC.prob;
% SMA
stat_cued_SMA = ft_timelockstatistics(cfg, nirs_singlecued_SMA{:},...
    nirs_dualcued_SMA{:});
stats.cued_SMA = stat_cued_SMA.prob;
% M1
stat_cued_M1 = ft_timelockstatistics(cfg, nirs_singlecued_M1{:},...
    nirs_dualcued_M1{:});
stats.cued_M1 = stat_cued_M1.prob;
% PPC
stat_cued_PPC = ft_timelockstatistics(cfg, nirs_singlecued_PPC{:},...
    nirs_dualcued_PPC{:});
stats.cued_PPC = stat_cued_PPC.prob;

disp(stats);

save(strcat(statistics_path, '\stats_nirs.mat'), 'stats');
