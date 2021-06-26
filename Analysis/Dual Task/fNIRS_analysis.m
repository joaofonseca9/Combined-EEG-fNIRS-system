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

subrec = ["28" "02"; "64" "01"];
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
    taskbaseline = {['sub-',char(sub),'_rec-',char(rec),'_dualcued_baseline'], ['sub-',char(sub),'_rec-',char(rec),'_singlecued_baseline'], ['sub-',char(sub),'_rec-',char(rec),'_dualuncued_baseline'], ['sub-',char(sub),'_rec-',char(rec),'_singleuncued_baseline']};
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
        nirs_TLblc{con} = ft_timelockbaseline(cfg, nirs_TL{con});
    end
    cd(fullfile(results_path, ['sub-',char(sub)]));
    save('nirs_TLblc.mat','nirs_TLblc');
    
    disp(['These are the results for subject ', char(sub), '.']);
    disp('Press any key to move onto the next subject.');
    pause;
    close all;
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
        nirs_all{con}{subject} = nirs_TLblc{con};
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

% Plot both (O2Hb and HHB) on the layout
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
saveas(gcf, '_avg_timelock.png');

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

%% Average statistical testing
cd(fullfile(results_path));
for i = 1:4 % loop over the 4 conditions
  [stat_O2Hb, stat_HHb] = statistics_withinsubjects(subsavg{i}, 'subsavg', layout, i, taskname{i});
end

disp('This was the end of individual subjects.');
disp('These are the results for the average of all subjects.');

