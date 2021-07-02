%% Analysis of the fNIRS signal (Hemodynamic Response and Statistical Analysis)
clear; clc; close all;

%% Initialize data
addpath('C:\Users\catar\OneDrive - Universidade do Porto\Twente\Combined-EEG-fNIRS-system\Analysis');
addpath('C:\Users\catar\OneDrive - Universidade do Porto\Twente\Combined-EEG-fNIRS-system\Analysis\Dual Task');
laptop = 'laptopCatarina';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);
results_path = 'C:\Users\catar\OneDrive - Universidade do Porto\Twente\Data Analysis\nirs';

ft_defaults;
[~, ftpath] = ft_version;

subrec = ["28" "02";"64" "01";"02" "02";"76" "01"];
conditions = [3 4 7 8];

%% Load data + processing per subject
% Go through all subjects
for subject = 1:size(subrec, 1)
    sub = subrec(subject, 1);
    rec = subrec(subject, 2);
  
    % Load fNIRS preprocessed signal
    load([mainpath_in, '\pre-processed\sub-', char(sub), '\nirs\sub-',...
        char(sub), '_rec-', char(rec), '_nirs_preprocessed.mat'], 'nirs_preprocessed');
   
    % Load layout 
    load(fullfile(mainpath_out,['sub-',char(sub)],'3d','layout.mat'), 'layout');
    
    % Keep the trials of interest (Dual Cued, Single Cued, Dual Uncued, Single
    % Uncued
    nirs = keepTrialsInterest(nirs_preprocessed);
    
    %% Baseline correction (preprocessing) + topoplot (spatial representation) per subject
    % Get the baseline and topoplot of all conditions for the subject
    taskname = {'Dual Cued', 'Single Cued', 'Dual Uncued', 'Single Uncued'};
    h = multiplot_condition(nirs, layout, conditions, taskname,...
        'baseline', [-10 0], 'trials', false, 'topoplot', 'yes', 'ylim',...
        [-0.2 0.2]);
    taskbaseline = {['sub-',char(sub),'_rec-',char(rec),'_dualcued_timelock'], ['sub-',char(sub),'_rec-',char(rec),'_singlecued_timelock'], ['sub-',char(sub),'_rec-',char(rec),'_dualuncued_timelock'], ['sub-',char(sub),'_rec-',char(rec),'_singleuncued_timelock']};
    tasktopoplotO2Hb = {['sub-',char(sub),'_rec-',char(rec),'_dualcued_topoplotO2Hb'], ['sub-',char(sub),'_rec-',char(rec),'_singlecued_topoplotO2Hb'], ['sub-',char(sub),'_rec-',char(rec),'_dualuncued_topoplotO2Hb'], ['sub-',char(sub),'_rec-',char(rec),'_singleuncued_topoplotO2Hb']};
    tasktopoplotHHb = {['sub-',char(sub),'_rec-',char(rec),'_dualcued_topoplotHHb'], ['sub-',char(sub),'_rec-',char(rec),'_singlecued_topoplotHHb'], ['sub-',char(sub),'_rec-',char(rec),'_dualuncued_topoplotHHb'], ['sub-',char(sub),'_rec-',char(rec),'_singleuncued_topoplotHHb']};
    
    % Save figures
    cd(fullfile(results_path, ['sub-',char(sub)], 'baseline + topoplot'));
    saveas(h{1},taskbaseline{1},'png'); saveas(h{2},tasktopoplotO2Hb{1},'png'); saveas(h{3},tasktopoplotHHb{1},'png'); 
    saveas(h{4},taskbaseline{2},'png'); saveas(h{5},tasktopoplotO2Hb{2},'png'); saveas(h{6},tasktopoplotHHb{2},'png'); 
    saveas(h{7},taskbaseline{3},'png'); saveas(h{8},tasktopoplotO2Hb{3},'png'); saveas(h{9},tasktopoplotHHb{3},'png');
    saveas(h{10},taskbaseline{4},'png'); saveas(h{11},tasktopoplotO2Hb{4},'png'); saveas(h{12},tasktopoplotHHb{4},'png'); 

    %% Statistical testing per subject 
    cd(fullfile(results_path,['sub-',char(sub)],'statistics'));
    for i=1:4 % loop over the 4 conditions
      [stat_O2Hb, stat_HHb] = statistics_withinsubjects(nirs, 'nirs', layout, i, taskname{i}, char(sub), char(rec));
    end
       
    %% Timelock analysis per subject
    for con = 1:length(conditions) % 4 conditions
        cfg = [];
        cfg.trials = find(nirs.trialinfo(:,1) == conditions(con)); % average the data for given task
        nirs_TL{con} = ft_timelockanalysis(cfg, nirs);
    end
    cd(fullfile(results_path, ['sub-',char(sub)]));
    save('nirs_TL.mat', 'nirs_TL');
    
    %% Baseline correction per subject
    for con = 1:length(conditions)
        cfg = [];
        cfg.baseline = [-10 0]; % define the amount of seconds you want to use for the baseline
        nirs_TLblc{con}{subject} = ft_timelockbaseline(cfg, nirs_TL{con});
    end
    
    disp(['These are the results for subject ', char(sub), '.']);
    disp('Press any key to move onto the next subject.');
    %pause;
    close all;
end

disp('This was the end of individual subjects.');

%%
for con = 1:length(conditions)
for subject = 1:size(subrec, 1)
    sub = subrec(subject, 1);
    rec = subrec(subject, 2);
    if strcmp(sub,"02")
        nirs_TLblc{con}{subject}.label = nirs_TLblc{1}{2}.label;
        nirs_TLblc{con}{subject}.cfg = nirs_TLblc{1}{2}.cfg;
        nirs_TLblc{con}{subject}.dof(1:46,1:length(nirs_TLblc{con}{subject}.time)) = 10;
        nirs_TLblc{con}{subject}.avg = [nirs_TLblc{con}{subject}.avg([1:39-1],:);zeros(1,length(nirs_TLblc{con}{subject}.time));nirs_TLblc{con}{subject}.avg(39:end,:)];
        nirs_TLblc{con}{subject}.avg = [nirs_TLblc{con}{subject}.avg([1:40-1],:);zeros(1,length(nirs_TLblc{con}{subject}.time));nirs_TLblc{con}{subject}.avg(40:end,:)];
    end
    if strcmp(sub,"76")
        nirs_TLblc{con}{subject}.label = nirs_TLblc{1}{1}.label;
        nirs_TLblc{con}{subject}.cfg = nirs_TLblc{1}{1}.cfg;
        nirs_TLblc{con}{subject}.dof(1:46,1:length(nirs_TLblc{con}{subject}.time)) = 10;
        nirs_TLblc{con}{subject}.avg(45,:) = 0;
        nirs_TLblc{con}{subject}.avg(46,:) = 0;
    end
    cd(fullfile(results_path, ['sub-',char(sub)]));
    save('nirs_TLblc.mat','nirs_TLblc');
end
end

%% Average the hemodynamic responses over all subjects
% Store baseline and timelockanalysis data of all subjects into one cell array
clear nirs_all 
for subject = 1:size(subrec, 1)
    sub = subrec(subject, 1);
    rec = subrec(subject, 2);
    for con = 1:length(conditions)
        cd()
        load(fullfile(results_path, ['sub-',char(sub)], 'nirs_TLblc'));
        nirs_all{con}{subject} = nirs_TLblc{con}{subject};
    end 
end
 
% Average over all subjects (for each condition seperately)
for con = 1:length(conditions)
    cfg = [];
    subsavg{con} = ft_timelockgrandaverage(cfg, nirs_all{con}{:});
end
cd(results_path);
save('subsavg.mat', 'subsavg');

%% Plot averaged data
% Separate O2Hb and HHb channels
for con = 1:length(conditions)
    cfg = [];
    cfg.channel = '* [O2Hb]';
    nirs_TLO2Hb{con} = ft_selectdata(cfg, subsavg{con});
    
    % rename labels such that they have the same name as HHb channels
    for i = 1:length(nirs_TLO2Hb{con}.label)
        tmp = strsplit(nirs_TLO2Hb{con}.label{i});
        nirs_TLO2Hb{con}.label{i}=tmp{1};
    end
    cd(results_path);
    save('nirs_TLO2Hb.mat','nirs_TLO2Hb');
    
    % the same for HHb channels
    cfg = [];
    cfg.channel = '* [HHb]';
    nirs_TLHHb{con} = ft_preprocessing(cfg, subsavg{con});
    for i=1:length(nirs_TLHHb{con}.label)
        tmp = strsplit(nirs_TLHHb{con}.label{i});
        nirs_TLHHb{con}.label{i}=tmp{1};
    end
    cd(results_path);
    save('nirs_TLHHb.mat','nirs_TLHHb');
end

%% Topoplots for each condition
cfg          = [];
cfg.layout   = layout;
cfg.marker   = 'labels';
cfg.ylim     = [-0.2 0.2];
cfg.xlim     = [5 10];
cfg.zlim     = cfg.ylim/2;
% Choose the time window over which you want to average
for con=1:4
    figure;
    title(taskname{con})
    ft_topoplotER(cfg, nirs_TLO2Hb{con});
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf,fullfile(results_path,['topoplot_',taskname{con},'.png']))
end

%% Plot both (O2Hb and HHB) on the layout
cfg = [];
cfg.showlabels = 'yes';
cfg.layout = layout;
cfg.interactive = 'yes'; % this allows to select a subplot and interact with it
cfg.linecolor = 'rbrbmcmc'; % O2Hb is showed in red (cued) and magenta (uncued), HHb in blue (cued) and cyan (uncued)
cfg.linestyle = {'--', '--', '-', '-', ':', ':', '-.', '-.'}; 
cfg.comment = 'dual cued is dashed line, single cued is solid line, dual uncued is dotted line and single uncued is a dashed dot line';
figure;
ft_multiplotER(cfg, nirs_TLO2Hb{1}, nirs_TLHHb{1}, nirs_TLO2Hb{2}, nirs_TLHHb{2}, nirs_TLO2Hb{3}, nirs_TLHHb{3}, nirs_TLO2Hb{4}, nirs_TLHHb{4});
cd(results_path);
saveas(gcf, 'avg_timelock.png');

% Plot for each task separately
for con = 1:length(conditions)
    cfg = [];
    cfg.showlabels = 'yes';
    cfg.layout = layout;
    cfg.showoutline = 'yes';
    cfg.interactive = 'yes'; % this allows to select a subplot and interact with it
    cfg.linecolor = 'rb'; % O2Hb is showed in red, HHb in blue
    figure; 
    title(taskname{con}); 
    ft_multiplotER(cfg, nirs_TLO2Hb{con}, nirs_TLHHb{con})
    saveas(gcf, [char(taskname{con}) '_avg_timelock.png']);
end

%% Extract channels belonging to the specific regions
cd(fullfile(results_path, 'areas'));

% DLPFC: Rx5-Tx7, Rx5-Tx8, Rx7-Tx7, Rx7-Tx8, Rx9-Tx13, Rx9-Tx12, Rx11-Tx12
cfg = [];
cfg.channel = {'Rx5-Tx7', 'Rx5-Tx8', 'Rx7-Tx7', 'Rx7-Tx8', 'Rx9-Tx13',...
    'Rx9-Tx12', 'Rx11-Tx12'};
nirs_HbO2_DLPFC{1} = ft_selectdata(cfg, nirs_TLO2Hb{1});
nirs_HbO2_DLPFC{2} = ft_selectdata(cfg, nirs_TLO2Hb{2});
nirs_HbO2_DLPFC{3} = ft_selectdata(cfg, nirs_TLO2Hb{3});
nirs_HbO2_DLPFC{4} = ft_selectdata(cfg, nirs_TLO2Hb{4});
nirs_Hb_DLPFC{1} = ft_selectdata(cfg, nirs_TLHHb{1});
nirs_Hb_DLPFC{2} = ft_selectdata(cfg, nirs_TLHHb{2});
nirs_Hb_DLPFC{3} = ft_selectdata(cfg, nirs_TLHHb{3});
nirs_Hb_DLPFC{4} = ft_selectdata(cfg, nirs_TLHHb{4});

save('nirs_HbO2_DLPFC.mat', 'nirs_HbO2_DLPFC');
save('nirs_Hb_DLPFC.mat', 'nirs_Hb_DLPFC');

% SMA: Rx4-Tx5, Rx3-Tx5, Rx4-Tx4
cfg = [];
cfg.channel = {'Rx4-Tx5', 'Rx3-Tx5', 'Rx4-Tx4'};
nirs_HbO2_SMA{1} = ft_selectdata(cfg, nirs_TLO2Hb{1});
nirs_HbO2_SMA{2} = ft_selectdata(cfg, nirs_TLO2Hb{2});
nirs_HbO2_SMA{3} = ft_selectdata(cfg, nirs_TLO2Hb{3});
nirs_HbO2_SMA{4} = ft_selectdata(cfg, nirs_TLO2Hb{4});
nirs_Hb_SMA{1} = ft_selectdata(cfg, nirs_TLHHb{1});
nirs_Hb_SMA{2} = ft_selectdata(cfg, nirs_TLHHb{2});
nirs_Hb_SMA{3} = ft_selectdata(cfg, nirs_TLHHb{3});
nirs_Hb_SMA{4} = ft_selectdata(cfg, nirs_TLHHb{4});

save('nirs_HbO2_SMA.mat', 'nirs_HbO2_SMA');
save('nirs_Hb_SMA.mat', 'nirs_Hb_SMA');

% M1: Rx3-Tx2, Rx1-Tx2, Rx3-Tx3, Rx1-Tx3, Rx2-Tx4, Rx2-Tx3
cfg = [];
cfg.channel = {'Rx3-Tx2', 'Rx1-Tx2', 'Rx3-Tx3', 'Rx1-Tx3', 'Rx2-Tx4',...
    'Rx2-Tx3'};
nirs_HbO2_M1{1} = ft_selectdata(cfg, nirs_TLO2Hb{1});
nirs_HbO2_M1{2} = ft_selectdata(cfg, nirs_TLO2Hb{2});
nirs_HbO2_M1{3} = ft_selectdata(cfg, nirs_TLO2Hb{3});
nirs_HbO2_M1{4} = ft_selectdata(cfg, nirs_TLO2Hb{4});
nirs_Hb_M1{1} = ft_selectdata(cfg, nirs_TLHHb{1});
nirs_Hb_M1{2} = ft_selectdata(cfg, nirs_TLHHb{2});
nirs_Hb_M1{3} = ft_selectdata(cfg, nirs_TLHHb{3});
nirs_Hb_M1{4} = ft_selectdata(cfg, nirs_TLHHb{4});

save('nirs_HbO2_M1.mat', 'nirs_HbO2_M1');
save('nirs_Hb_M1.mat', 'nirs_Hb_M1');

% PPC: Rx8-Tx10, Rx6-Tx9, Rx8-Tx9, Rx12-Tx15, Rx10-Tx14, Rx12-Tx14
cfg = [];
cfg.channel = {'Rx3-Tx2', 'Rx1-Tx2', 'Rx3-Tx3', 'Rx1-Tx3', 'Rx2-Tx4',...
    'Rx2-Tx3'};
nirs_HbO2_PPC{1} = ft_selectdata(cfg, nirs_TLO2Hb{1});
nirs_HbO2_PPC{2} = ft_selectdata(cfg, nirs_TLO2Hb{2});
nirs_HbO2_PPC{3} = ft_selectdata(cfg, nirs_TLO2Hb{3});
nirs_HbO2_PPC{4} = ft_selectdata(cfg, nirs_TLO2Hb{4});
nirs_Hb_PPC{1} = ft_selectdata(cfg, nirs_TLHHb{1});
nirs_Hb_PPC{2} = ft_selectdata(cfg, nirs_TLHHb{2});
nirs_Hb_PPC{3} = ft_selectdata(cfg, nirs_TLHHb{3});
nirs_Hb_PPC{4} = ft_selectdata(cfg, nirs_TLHHb{4});

save('nirs_HbO2_PPC.mat', 'nirs_HbO2_PPC');
save('nirs_Hb_PPC.mat', 'nirs_Hb_PPC');

%% Average the responses over all channels of specific regions and plot them
% DLPFC
for con = 1:length(conditions)
    regionsavg_HbO2_DLPFC{con} = mean(nirs_HbO2_DLPFC{con}.avg, 1);
    regionsavg_Hb_DLPFC{con} = mean(nirs_Hb_DLPFC{con}.avg, 1);
    
    figure; title(char(taskname{con}));
    plot(nirs_HbO2_DLPFC{con}.time, regionsavg_HbO2_DLPFC{con}, 'r');
    hold on;
    plot(nirs_Hb_DLPFC{con}.time, regionsavg_Hb_DLPFC{con}, 'b');
    hold on;
    xline(0);
    hold off;
    legend('Hb02', 'Hb');
    xlim([-10 20]);
    ylim([-0.4 0.4])
    
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, [char(taskname{con}) '_DLPFC.png']);
    
end

% SMA
for con = 1:length(conditions)
    regionsavg_HbO2_SMA{con} = mean(nirs_HbO2_SMA{con}.avg, 1);
    regionsavg_Hb_SMA{con} = mean(nirs_Hb_SMA{con}.avg, 1);
    
    figure; title(char(taskname{con}));
    plot(nirs_HbO2_SMA{con}.time, regionsavg_HbO2_SMA{con}, 'r');
    hold on;
    plot(nirs_Hb_SMA{con}.time, regionsavg_Hb_SMA{con}, 'b');
    hold on;
    xline(0);
    hold off;
    legend('Hb02', 'Hb');
    xlim([-5 20]);
    ylim([-0.2 0.25])
    
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, [char(taskname{con}) '_SMA.png']);
    
end

% M1
for con = 1:length(conditions)
    regionsavg_HbO2_M1{con} = mean(nirs_HbO2_M1{con}.avg, 1);
    regionsavg_Hb_M1{con} = mean(nirs_Hb_M1{con}.avg, 1);
    
    figure; title(char(taskname{con}));
    plot(nirs_HbO2_M1{con}.time, regionsavg_HbO2_M1{con}, 'r');
    hold on;
    plot(nirs_Hb_M1{con}.time, regionsavg_Hb_M1{con}, 'b');
    hold on;
    xline(0);
    hold off;
    legend('Hb02', 'Hb');
    xlim([-5 20]);
    ylim([-0.15 0.2])
    
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, [char(taskname{con}) '_M1.png']);
    
end

% PPC
for con = 1:length(conditions)
    regionsavg_HbO2_PPC{con} = mean(nirs_HbO2_PPC{con}.avg, 1);
    regionsavg_Hb_PPC{con} = mean(nirs_Hb_PPC{con}.avg, 1);
    
    figure; title(char(taskname{con}));
    plot(nirs_HbO2_PPC{con}.time, regionsavg_HbO2_PPC{con}, 'r');
    hold on;
    plot(nirs_Hb_PPC{con}.time, regionsavg_Hb_PPC{con}, 'b');
    hold on;
    xline(0);
    hold off;
    legend('Hb02', 'Hb');
    xlim([-5 20]);
    ylim([-0.1 0.15])
    
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, [char(taskname{con}) '_PPC.png']);
    
end

disp('These are the results for the average of all subjects.');

