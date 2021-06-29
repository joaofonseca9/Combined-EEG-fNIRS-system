%% Analysis of the fNIRS signals - hemodynamic response.

clear; clc; close all;
addpath('C:\Users\maria\OneDrive\Documentos\GitHub\Combined-EEG-fNIRS-system\Analysis');

laptop = 'laptopMariana';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);
results_path = 'C:\Users\maria\OneDrive\Ambiente de Trabalho\Automaticity Results\Hemodynamic Response';
analysis_path = 'C:\Users\maria\OneDrive\Documentos\GitHub\Combined-EEG-fNIRS-system\Analysis\Automatic_vs_NonAutomatic';
addpath(analysis_path);

subrec = ["28" "02"; "02" "02"];
conditions = [2 4 6 8];

% Loop through every subject.
for subject = 1:size(subrec, 1)
    sub = subrec(subject, 1);
    rec = subrec(subject, 2);
    
    % Load the subject's fNIRS signals.
    load([mainpath_in, '\pre-processed\sub-', char(sub), '\nirs\sub-',...
        char(sub), '_rec-', char(rec), '_nirs_preprocessed.mat']);
    
    % Load the layout of the optode's template.
    load(fullfile(mainpath_out,['sub-',char(sub)],'3d','layout.mat'),...
        'layout');
    
    % Keep only the trials of interest (Auto Cued, Non-Auto Cued, Auto
    % Uncued, Non-Auto Uncued).
    nirs = keepTrialsInterest(nirs_preprocessed);
    
    %% Baseline correction + plots
    % Get the baseline and topoplot of all conditions for the subject
    taskname = {'Auto Cued', 'Non-Auto Cued', 'Auto Uncued',...
        'Non-Auto Uncued'};
    h = multiplot_condition(nirs, layout, [2 4 6 8], taskname,...
        'baseline', [-10 0], 'trials', false, 'topoplot', 'yes',...
        'ylim', [-0.5 1]);
    taskbaseline = {['sub-', char(sub), '_rec-',char(rec),...
        '_autocued_baseline'], ['sub-', char(sub), '_rec-', char(rec),...
        '_nonautocued_baseline'], ['sub-', char(sub), '_rec-', char(rec),...
        '_autouncued_baseline'], ['sub-', char(sub), '_rec-', char(rec),...
        '_nonautouncued_baseline']};
    tasktopoplotO2Hb = {['sub-', char(sub), '_rec-', char(rec),...
        '_autocued_topoplotO2Hb'], ['sub-', char(sub), '_rec-', char(rec),...
        '_nonautocued_topoplotO2Hb'], ['sub-', char(sub), '_rec-',...
        char(rec), '_autouncued_topoplotO2Hb'], ['sub-', char(sub),...
        '_rec-', char(rec), '_nonautouncued_topoplotO2Hb']};
    tasktopoplotHHb = {['sub-', char(sub), '_rec-', char(rec),...
        '_autocued_topoplotHHb'], ['sub-', char(sub), '_rec-', char(rec),...
        '_nonautocued_topoplotHHb'], ['sub-', char(sub), '_rec-',...
        char(rec), '_autouncued_topoplotHHb'], ['sub-', char(sub),...
        '_rec-', char(rec), '_nonautouncued_topoplotHHb']};
    
    % Save figures.
    cd(fullfile(results_path, ['Sub-', char(sub)], 'Plots'));
    saveas(h{1}, taskbaseline{1}, 'png'); 
    saveas(h{2}, tasktopoplotO2Hb{1}, 'png'); 
    saveas(h{3}, tasktopoplotHHb{1}, 'png');
    saveas(h{4}, taskbaseline{2}, 'png'); 
    saveas(h{5}, tasktopoplotO2Hb{2}, 'png'); 
    saveas(h{6}, tasktopoplotHHb{2}, 'png');
    saveas(h{7}, taskbaseline{3}, 'png'); 
    saveas(h{8}, tasktopoplotO2Hb{3}, 'png'); 
    saveas(h{9}, tasktopoplotHHb{3}, 'png');
    saveas(h{10}, taskbaseline{4}, 'png'); 
    saveas(h{11}, tasktopoplotO2Hb{4}, 'png'); 
    saveas(h{12}, tasktopoplotHHb{4}, 'png');
    
    %% Statistical testing per subject
    cd(fullfile(results_path, ['Sub-',char(sub)], 'Statistics'));
    for i=1:4
        [stat_O2Hb, stat_HHb, h] = statistics_withinsubjects(nirs, 'nirs', layout, i, taskname{i}, char(sub), char(rec));
    end
    
    %% Timelock analysis per subject
    cd(fullfile(results_path, ['Sub-',char(sub)], 'Timelock Analysis'));
    for con = 1:4 % 4 conditions
        cfg = [];
        cfg.trials = find(nirs.trialinfo(:,1) == conditions(con)); % average the data for given task
        nirs_TL{con} = ft_timelockanalysis(cfg, nirs);
    end
    save('nirs_TL.mat', 'nirs_TL');
  
    %% Timelock analysis for baseline.
    for con = 1:length(conditions)
        cfg = [];
        cfg.baseline = [-10 0]; % define the amount of seconds you want to use for the baseline
        nirs_TLblc{con} = ft_timelockbaseline(cfg, nirs_TL{con});
    end
    save('nirs_TLblc.mat','nirs_TLblc');
    
    disp(['These are the results for subject ', char(sub), '.']);
    disp('Press any key to move onto the next subject.');
    pause;
    close all;
    
end

%% Average the hemodynamic responses over all subjects.
% Store baseline and timelockanalysis data of all subjects into one cell
% array.
clear nirs_all
for subject = 1:size(subrec, 1)
    sub = subrec(subject, 1);
    rec = subrec(subject, 2);
    for con = 1:length(conditions)
        cd()
        load(fullfile(results_path, ['Sub-', char(sub)],...
            'Timelock Analysis\nirs_TLblc.mat'), 'nirs_TLblc');
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
    
    % Rename labels such that they have the same name as HHb channels.
    for i = 1:length(nirs_TLO2Hb{con}.label)
        tmp = strsplit(nirs_TLO2Hb{con}.label{i});
        nirs_TLO2Hb{con}.label{i}=tmp{1};
    end
    cd(results_path);
    save('nirs_TLO2Hb.mat', 'nirs_TLO2Hb');
    
    % The same for HHb channels.
    cfg = [];
    cfg.channel = '* [HHb]';
    nirs_TLHHb{con} = ft_preprocessing(cfg, subsavg{con});
    for i=1:length(nirs_TLHHb{con}.label)
        tmp = strsplit(nirs_TLHHb{con}.label{i});
        nirs_TLHHb{con}.label{i}=tmp{1};
    end
    cd(results_path);
    save('nirs_TLHHb.mat', 'nirs_TLHHb');
end

% Plot both (O2Hb and HHB) on the layout.
cfg = [];
cfg.showlabels = 'yes';
cfg.layout = layout;
% This allows to select a subplot and interact with it.
cfg.interactive = 'yes';
% O2Hb is showed in red (cued) and magenta (uncued).
% HHb in blue (cued) and cyan (uncued).
cfg.linecolor = 'rbrbmcmc';
cfg.linestyle = {'--', '--', '-', '-', ':', ':', '-.', '-.'};
cfg.comment = 'auto cued is dashed line; non-auto cued is solid line, auto uncued is dotted line and non-auto uncued is a dashed dot line';
figure;
ft_multiplotER(cfg, nirs_TLO2Hb{1}, nirs_TLHHb{1}, nirs_TLO2Hb{2},...
    nirs_TLHHb{2}, nirs_TLO2Hb{3}, nirs_TLHHb{3}, nirs_TLO2Hb{4},...
    nirs_TLHHb{4});
cd(results_path);
saveas(gcf, 'avg_timelock.png');

% Plot for each task separately.
for con = 1:length(conditions)
    cfg = [];
    cfg.showlabels = 'yes';
    cfg.layout = layout;
    cfg.showoutline = 'yes';
    % This allows to select a subplot and interact with it.
    cfg.interactive = 'yes';
    % O2Hb is showed in red, HHb in blue.
    cfg.linecolor = 'rb'; 
    figure;
    title(taskname{con});
    ft_multiplotER(cfg, nirs_TLO2Hb{con}, nirs_TLHHb{con})
    saveas(gcf, [char(taskname{con}) '_avg_timelock.png']);
end

disp('This was the end of individual subjects.');
disp('These are the results for the average of all subjects.');

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