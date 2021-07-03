clear; clc; close all;
addpath('C:\Users\maria\OneDrive\Documentos\GitHub\Combined-EEG-fNIRS-system\Analysis');

laptop = 'laptopMariana';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);

results_path = 'C:\Users\maria\OneDrive\Ambiente de Trabalho\Automaticity Results\Hemodynamic Response';
statistics_path = 'C:\Users\maria\OneDrive\Ambiente de Trabalho\Statistics\Hemodynamic Response';

subrec = ["28" "02"; "02" "02"; "76" "01"];
conditions = [2 4 6 8];
taskname = {'Auto Cued', 'Non-Auto Cued', 'Auto Uncued', 'Non-Auto Uncued'};

for subject = 1:size(subrec, 1)
    sub = subrec(subject, 1);
    rec = subrec(subject, 2);
    
    % Load the subject's file.
    load(fullfile(results_path, ['Sub-', char(sub)],...
            'Timelock Analysis\nirs_TLblc.mat'), 'nirs_TLblc');
        
    % Separate Hb02 from Hb.    
    [nirs_TLO2Hb, nirs_TLHHb] = separateHbO2FromHb(conditions, nirs_TLblc, subject);
    
    % Extract the different regions of interest: DLPFC, SMA, M1 and PPC.
    % Only for HbO2 signal.
    [nirs_HbO2_DLPFC, nirs_HbO2_SMA, nirs_HbO2_M1, nirs_HbO2_PPC] =...
        extractROIs(nirs_TLO2Hb);
    
    % Separate into the different conditions.
    % DLPFC.
    [nirs_autocued, nirs_nonautocued, nirs_autouncued, nirs_nonautouncued]...
        = extractConditions(conditions, nirs_HbO2_DLPFC);
    nirs_autocued_DLPFC{subject} = nirs_autocued;
    nirs_nonautocued_DLPFC{subject} = nirs_nonautocued;
    nirs_autouncued_DLPFC{subject} = nirs_autouncued;
    nirs_nonautouncued_DLPFC{subject} = nirs_nonautouncued;
    % SMA.
    [nirs_autocued, nirs_nonautocued, nirs_autouncued, nirs_nonautouncued]...
        = extractConditions(conditions, nirs_HbO2_SMA);
    nirs_autocued_SMA{subject} = nirs_autocued;
    nirs_nonautocued_SMA{subject} = nirs_nonautocued;
    nirs_autouncued_SMA{subject} = nirs_autouncued;
    nirs_nonautouncued_SMA{subject} = nirs_nonautouncued;
    % M1.
    [nirs_autocued, nirs_nonautocued, nirs_autouncued, nirs_nonautouncued]...
        = extractConditions(conditions, nirs_HbO2_M1);
    nirs_autocued_M1{subject} = nirs_autocued;
    nirs_nonautocued_M1{subject} = nirs_nonautocued;
    nirs_autouncued_M1{subject} = nirs_autouncued;
    nirs_nonautouncued_M1{subject} = nirs_nonautouncued;
    % PPC.
    [nirs_autocued, nirs_nonautocued, nirs_autouncued, nirs_nonautouncued]...
        = extractConditions(conditions, nirs_HbO2_PPC);
    nirs_autocued_PPC{subject} = nirs_autocued;
    nirs_nonautocued_PPC{subject} = nirs_nonautocued;
    nirs_autouncued_PPC{subject} = nirs_autouncued;
    nirs_nonautouncued_PPC{subject} = nirs_nonautouncued;
    
end

%% T-TEST.

% Define the parameters for the statistical comparison.
cfg = [];
cfg.latency = 'all';
cfg.avgoverchan = 'yes';
cfg.avgovertime = 'yes';
cfg.parameter = 'avg';
cfg.method = 'analytic';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.alpha = 0.05;
cfg.correctm = 'no';

Nsub = 3;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
% The 1st row in cfg.design contains the independent variable.
cfg.ivar = 1; 
% The 2nd row in cfg.design contains the subject number.
cfg.uvar = 2; 

% Auto Uncued vs Cued.
% DLPFC.
stat_auto_DLPFC = ft_timelockstatistics(cfg, nirs_autouncued_DLPFC{:},...
    nirs_autocued_DLPFC{:});
stats.auto_DLPFC = stat_auto_DLPFC.prob;
% SMA.
stat_auto_SMA = ft_timelockstatistics(cfg, nirs_autouncued_SMA{:},...
    nirs_autocued_SMA{:});
stats.auto_SMA = stat_auto_SMA.prob;
% M1.
stat_auto_M1 = ft_timelockstatistics(cfg, nirs_autouncued_M1{:},...
    nirs_autocued_M1{:});
stats.auto_M1 = stat_auto_M1.prob;
% PPC.
stat_auto_PPC = ft_timelockstatistics(cfg, nirs_autouncued_PPC{:},...
    nirs_autocued_PPC{:});
stats.auto_PPC = stat_auto_PPC.prob;

% Non-Auto Uncued vs Cued.
% DLPFC.
stat_nonauto_DLPFC = ft_timelockstatistics(cfg, nirs_nonautouncued_DLPFC{:},...
    nirs_nonautocued_DLPFC{:});
stats.nonauto_DLPFC = stat_nonauto_DLPFC.prob;
% SMA.
stat_nonauto_SMA = ft_timelockstatistics(cfg, nirs_nonautouncued_SMA{:},...
    nirs_nonautocued_SMA{:});
stats.nonauto_SMA = stat_nonauto_SMA.prob;
% M1.
stat_nonauto_M1 = ft_timelockstatistics(cfg, nirs_nonautouncued_M1{:},...
    nirs_nonautocued_M1{:});
stats.nonauto_M1 = stat_nonauto_M1.prob;
% PPC.
stat_nonauto_PPC = ft_timelockstatistics(cfg, nirs_nonautouncued_PPC{:},...
    nirs_nonautocued_PPC{:});
stats.nonauto_PPC = stat_nonauto_PPC.prob;

save(strcat(statistics_path, '\stats_nirs.mat'), 'stats');

%% Functions.
function [nirs_TLO2Hb, nirs_TLHHb] = separateHbO2FromHb(conditions, nirs_TLblc, subject)
% Separate O2Hb and HHb channels.

for con = 1:length(conditions)
    cfg = [];
    cfg.channel = '* [O2Hb]';
    nirs_TLO2Hb{con} = ft_selectdata(cfg, nirs_TLblc{con}{subject});
    
    % Rename labels such that they have the same name as HHb channels.
    for i = 1:length(nirs_TLO2Hb{con}.label)
        tmp = strsplit(nirs_TLO2Hb{con}.label{i});
        nirs_TLO2Hb{con}.label{i}=tmp{1};
    end
    
    % The same for HHb channels.
    cfg = [];
    cfg.channel = '* [HHb]';
    nirs_TLHHb{con} = ft_preprocessing(cfg, nirs_TLblc{con}{subject});
    for i=1:length(nirs_TLHHb{con}.label)
        tmp = strsplit(nirs_TLHHb{con}.label{i});
        nirs_TLHHb{con}.label{i}=tmp{1};
    end

end
end

function [nirs_HbO2_DLPFC, nirs_HbO2_SMA, nirs_HbO2_M1, nirs_HbO2_PPC] =...
    extractROIs(nirs_TLO2Hb)

% DLPFC: Rx5-Tx7, Rx5-Tx8, Rx7-Tx7, Rx7-Tx8, Rx9-Tx13, Rx9-Tx12, Rx11-Tx12,
% Rx11-Tx13.
cfg = [];
cfg.channel = {'Rx5-Tx7', 'Rx5-Tx8', 'Rx7-Tx7', 'Rx7-Tx8', 'Rx9-Tx13',...
    'Rx9-Tx12', 'Rx11-Tx12', 'Rx11-Tx13'};
nirs_HbO2_DLPFC{1} = ft_selectdata(cfg, nirs_TLO2Hb{1});
nirs_HbO2_DLPFC{2} = ft_selectdata(cfg, nirs_TLO2Hb{2});
nirs_HbO2_DLPFC{3} = ft_selectdata(cfg, nirs_TLO2Hb{3});
nirs_HbO2_DLPFC{4} = ft_selectdata(cfg, nirs_TLO2Hb{4});

% PMC/SMA: Rx4-Tx5, Rx3-Tx5, Rx4-Tx4.
cfg = [];
cfg.channel = {'Rx4-Tx5', 'Rx3-Tx5', 'Rx4-Tx4'};
nirs_HbO2_SMA{1} = ft_selectdata(cfg, nirs_TLO2Hb{1});
nirs_HbO2_SMA{2} = ft_selectdata(cfg, nirs_TLO2Hb{2});
nirs_HbO2_SMA{3} = ft_selectdata(cfg, nirs_TLO2Hb{3});
nirs_HbO2_SMA{4} = ft_selectdata(cfg, nirs_TLO2Hb{4});

% M1: Rx3-Tx2, Rx1-Tx2, Rx3-Tx3, Rx1-Tx3, Rx2-Tx4, Rx2-Tx3.
cfg = [];
cfg.channel = {'Rx3-Tx2', 'Rx1-Tx2', 'Rx3-Tx3', 'Rx1-Tx3', 'Rx2-Tx4',...
    'Rx2-Tx3'};
nirs_HbO2_M1{1} = ft_selectdata(cfg, nirs_TLO2Hb{1});
nirs_HbO2_M1{2} = ft_selectdata(cfg, nirs_TLO2Hb{2});
nirs_HbO2_M1{3} = ft_selectdata(cfg, nirs_TLO2Hb{3});
nirs_HbO2_M1{4} = ft_selectdata(cfg, nirs_TLO2Hb{4});

% PPC: Rx8-Tx10, Rx6-Tx9, Rx8-Tx9, Rx12-Tx15, Rx10-Tx14, Rx12-Tx14.
cfg = [];
cfg.channel = {'Rx8-Tx10', 'Rx6-Tx9', 'Rx8-Tx9', 'Rx12-Tx15',...
    'Rx10-Tx14', 'Rx12-Tx14'};
nirs_HbO2_PPC{1} = ft_selectdata(cfg, nirs_TLO2Hb{1});
nirs_HbO2_PPC{2} = ft_selectdata(cfg, nirs_TLO2Hb{2});
nirs_HbO2_PPC{3} = ft_selectdata(cfg, nirs_TLO2Hb{3});
nirs_HbO2_PPC{4} = ft_selectdata(cfg, nirs_TLO2Hb{4});

end

function [nirs_autocued, nirs_nonautocued, nirs_autouncued,...
    nirs_nonautouncued] = extractConditions(conditions, nirs_HbO2)
% Extract the different conditions to analyse.

for con = 1:length(conditions)
    
    if con==1
        nirs_autocued = nirs_HbO2{con};
    elseif con==2
        nirs_nonautocued = nirs_HbO2{con};
    elseif con==3
        nirs_autouncued = nirs_HbO2{con};
    elseif con==4
        nirs_nonautouncued = nirs_HbO2{con};
    end
end

end