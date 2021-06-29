clear; clc; close all;
addpath('C:\Users\maria\OneDrive\Documentos\GitHub\Combined-EEG-fNIRS-system\Analysis');

laptop = 'laptopMariana';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);
results_path = 'C:\Users\maria\OneDrive\Ambiente de Trabalho\Automaticity Results\Hemodynamic Response';
analysis_path = 'C:\Users\maria\OneDrive\Documentos\GitHub\Combined-EEG-fNIRS-system\Analysis\Automatic_vs_NonAutomatic';
addpath(analysis_path);

subrec = ["28" "02"];
conditions = [2 4 6 8];

% Loop through every subject.
for subject = 1:size(subrec, 1)
    sub = subrec(subject, 1);
    rec = subrec(subject, 2);
  
    % Load the subject's fNIRS signals.
    load([mainpath_in, '\pre-processed\sub-', char(sub), '\nirs\sub-',...
        char(sub), '_rec-', char(rec), '_nirs_preprocessed.mat']);
    
    % Load the layout of the optode's template.
    load(fullfile(mainpath_out,['sub-',char(sub)],'3d','layout.mat'), 'layout');
    
    % Keep only the trials of interest (Auto Cued, Non-Auto Cued, Auto
    % Uncued, Non-Auto Uncued)
    nirs = keepTrialsInterest(nirs_preprocessed);
    
    % Get the hemodynamic response of all conditions for the specific
    % subject.
    %% Baseline correction + plots
    % Get the baseline and topoplot of all conditions for the subject
    taskname = {'Auto Cued', 'Non-Auto Cued', 'Auto Uncued', 'Non-Auto Uncued'};
    [h,~,~] = multiplot_condition_modified(nirs, layout, [2 4 6 8], taskname,...
        'baseline', [-10 0], 'trials', false, 'topoplot', 'yes', 'ylim',...
        [-0.5 1]);
    taskbaseline = {['sub-',char(sub),'_rec-',char(rec),'_dualcued_baseline'], ['sub-',char(sub),'_rec-',char(rec),'_singlecued_baseline'], ['sub-',char(sub),'_rec-',char(rec),'_dualuncued_baseline'], ['sub-',char(sub),'_rec-',char(rec),'_singleuncued_baseline']};
    tasktopoplotO2Hb = {['sub-',char(sub),'_rec-',char(rec),'_dualcued_topoplotO2Hb'], ['sub-',char(sub),'_rec-',char(rec),'_singlecued_topoplotO2Hb'], ['sub-',char(sub),'_rec-',char(rec),'_dualuncued_topoplotO2Hb'], ['sub-',char(sub),'_rec-',char(rec),'_singleuncued_topoplotO2Hb']};
    tasktopoplotHHb = {['sub-',char(sub),'_rec-',char(rec),'_dualcued_topoplotHHb'], ['sub-',char(sub),'_rec-',char(rec),'_singlecued_topoplotHHb'], ['sub-',char(sub),'_rec-',char(rec),'_dualuncued_topoplotHHb'], ['sub-',char(sub),'_rec-',char(rec),'_singleuncued_topoplotHHb']};
    
    %% Statistical testing per subject 
    cd(fullfile(results_path, ['Sub-',char(sub)], 'Statistics'));
    for i=1:4
        %% Get the options
        condition = conditions(i);
        baseline = [-10 0];
        taskstatisticsO2Hb = {['sub-',sub,'_rec-',rec,'_dualcued_statisticsO2Hb'], ['sub-',sub,'_rec-',rec,'_singlecued_statisticsO2Hb'], ['sub-',sub,'_rec-',rec,'_dualuncued_statisticsO2Hb'], ['sub-',sub,'_rec-',rec,'_singleuncued_statisticsO2Hb']};
        taskstatisticsHHb = {['sub-',sub,'_rec-',rec,'_dualcued_statisticsHHb'], ['sub-',sub,'_rec-',rec,'_singlecued_statisticsHHb'], ['sub-',sub,'_rec-',rec,'_dualuncued_statisticsHHb'], ['sub-',sub,'_rec-',rec,'_singleuncued_statisticsHHb']};
        
        %% apply baseline correction
        cfg                 = [];
        cfg.baseline        = baseline; % define the amount of seconds you want to use for the baseline
        cfg.parameter = 'trial';
        data_blc = ft_timelockbaseline(cfg, nirs);
        
        % Collect data for the O2Hb and HHb channel during each trial of the
        % complex condition
        cfg=[];
        cfg.channel='*[O2Hb]';
        cfg.trials=[data_blc.trialinfo(:,1)]'==condition;
        data_complex_O2Hb=ft_selectdata(cfg, data_blc);
        
        cfg=[];
        cfg.trials=[data_blc.trialinfo(:,1)]'==condition;
        cfg.channel='*[HHb]';
        data_complex_HHb=ft_selectdata(cfg, data_blc);
        
        %% t-tests with fieldtrip
        
        n_trials=length(data_complex_O2Hb.trialinfo);
        % for oxy channels
        cfg = [];
        cfg.latency     = [5 10];
        cfg.avgovertime = 'yes';
        cfg.parameter   = 'trial';
        cfg.method      = 'stats';
        cfg.statistic   = 'ttest';
        cfg.alpha       = 0.05;
        % cfg.alpha       = 0.05/32; % bonferroni correction for mulitple comparison (short channels are not included)
        cfg.tail        = 1;
        cfg.correctm = 'yes';
        cfg.design(1, 1:n_trials) = ones(1,n_trials);
        cfg.ivar =1;
        
        stat_O2Hb = ft_timelockstatistics(cfg,data_complex_O2Hb);
        % plot
        data_complex_O2Hb_lay=data_complex_O2Hb;
        for i=1:length(data_complex_O2Hb.label)
            split_label=strsplit(data_complex_O2Hb.label{i});
            data_complex_O2Hb_lay.label{i}=split_label{1};
        end
        cfg = [];
        cfg.style     = 'blank';
        cfg.layout    = layout;
        cfg.highlight = 'on';
        cfg.highlightchannel = find(stat_O2Hb.mask);
        cfg.comment   = 'no';
        h=figure;
        ft_topoplotER(cfg, data_complex_O2Hb_lay); title(sprintf('significant channels (p=0.05) for %s vs baseline [O2Hb]', condition_name))
        saveas(h,taskstatisticsO2Hb{condition},'png');
        
        
        
        
%      [stat_O2Hb, stat_HHb, h] = statistics_withinsubjects(nirs, 'nirs', layout, i, taskname{i}, char(sub), char(rec));
    end
    
  
    disp(['These are the results for subject ', char(sub), '.']);
    disp('Press any key to move onto the next subject.');
    pause;
    close all;

end

%% Functions

function nirs_output = keepTrialsInterest(nirs_input)
% Go through the different trials and keep only the ones of interest:
% 2 - Auto Cued
% 4 - Non-Auto Cued
% 6 - Auto Uncued
% 8 - Non-Auto Uncued

nirs_output = nirs_input;
numRemovedTrials = 0;

% Go through the different trials.
for i=1:length(nirs_input.trialinfo)
    % If it is not a trial of interest.
    if nirs_input.trialinfo(i) ~= 2 && nirs_input.trialinfo(i) ~= 4 ...
            && nirs_input.trialinfo(i) ~= 6 ...
            && nirs_input.trialinfo(i) ~= 8
        % Eliminate that trial and all the unecessary information.
        nirs_output.trial(i-numRemovedTrials) = [];
        nirs_output.time(i-numRemovedTrials) = [];
        nirs_output.trialinfo(i-numRemovedTrials) = [];
        nirs_output.sampleinfo(i-numRemovedTrials, :) = [];
        % Increase number of removed trials.
        numRemovedTrials = numRemovedTrials+1;
    end 
end

end

function plotScales(hlim, vlim, hpos, vpos, width, height)

% the placement of all elements is identical
placement = {'hpos', hpos, 'vpos', vpos, 'width', width, 'height', height, 'hlim', hlim, 'vlim', vlim};

ft_plot_box([hlim vlim], placement{:}, 'edgecolor', 'k');

if hlim(1)<=0 && hlim(2)>=0
  ft_plot_line([0 0], vlim, placement{:}, 'color', 'k');
  ft_plot_text(0, vlim(1), '0  ', placement{:}, 'rotation', 90, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 8);
end

if vlim(1)<=0 && vlim(2)>=0
  ft_plot_line(hlim, [0 0], placement{:}, 'color', 'k');
  ft_plot_text(hlim(1), 0, '0  ', placement{:}, 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'middle', 'FontSize', 8);
end

ft_plot_text(hlim(1), vlim(1), [num2str(hlim(1), 3) ' '], placement{:}, 'rotation', 90, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top',    'FontSize', 8);
ft_plot_text(hlim(2), vlim(1), [num2str(hlim(2), 3) ' '], placement{:}, 'rotation', 90, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 8);
ft_plot_text(hlim(1), vlim(1), [num2str(vlim(1), 3) ' '], placement{:}, 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'bottom', 'FontSize', 8);
ft_plot_text(hlim(1), vlim(2), [num2str(vlim(2), 3) ' '], placement{:}, 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'top',    'FontSize', 8);
end