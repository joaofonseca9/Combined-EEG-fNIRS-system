%% Analysis of the EEG signals - separate channels.

clear; clc; close all;
addpath('C:\Users\maria\OneDrive\Documentos\GitHub\Combined-EEG-fNIRS-system\Analysis');

laptop = 'laptopMariana';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);
results_path = 'C:\Users\maria\OneDrive\Ambiente de Trabalho\Automaticity Results\Hemodynamic Response';

subrec = ["04" "01"];

% Loop through every subject.
for subject = 1:size(subrec, 1)
    sub = subrec(subject, 1);
    rec = subrec(subject, 2);
  
    % Load the subject's fNIRS signals.
    load([mainpath_in, '\pre-processed\sub-', char(sub), '\nirs\sub-',...
        char(sub), '_rec-', char(rec), '_nirs_preprocessed.mat']);
    load([mainpath_in, '\pre-processed\sub-', char(sub), '\nirs\sub-',...
        char(sub), '_rec-', char(rec), '_nirs_epoch.mat']);
    nirs_preprocessed.trialinfo = nirs_epoch.trialinfo;
    nirs_preprocessed.sampleinfo = nirs_epoch.sampleinfo;
    
    % Load the layout of the optode's template.
    load(fullfile(mainpath_out,['sub-',char(sub)],'3d','layout.mat'), 'layout');
    
    % Keep only the trials of interest (Auto Cued, Non-Auto Cued, Auto
    % Uncued, Non-Auto Uncued)
    nirs = keepTrialsInterest(nirs_preprocessed);
    
    % Get the hemodynamic response of all conditions for the specific
    % subject.
    taskname = {'Auto Cued', 'Non-Auto Cued', 'Auto Uncued', 'Non-Auto Uncued'};
    [h, xval_allconditions, yval_allconditions] =...
        multiplot_condition_modified(nirs, layout, [2 4 6 8], taskname,...
        'baseline', [-10 0], 'trials', false, 'topoplot', 'no', 'ylim',...
        [-0.5 1]);
    
    % Insert the subject's hemodynamic responses on a matrix of all
    % subjects.
    xval_allsubjects_allconditions(:, :, :, subject) = xval_allconditions; 
    yval_allsubjects_allconditions(:, :, :, :, subject) = yval_allconditions; 
    
    disp(['These are the results for subject ', char(sub), '.']);
    disp('Press any key to move onto the next subject.');
    pause;
    close all;

end

% Average the hemodynamic responses over all subjects.
xval_avgsubjects_allconditions = mean(xval_allsubjects_allconditions, 4); 
yval_avgsubjects_allconditions = mean(yval_allsubjects_allconditions, 5);

% For each condition.
for con=1:size(xval_avgsubjects_allconditions, 3)
    f=1;

    % Aelect the channels that are both in data and layout.
    data_labels=nirs.label;
    for i=1:length(data_labels)
        tmp = strsplit(data_labels{i});
        data_labels{i}=tmp{1};
    end
    % Take the subselection of channels that is contained.
    [selchan, sellay] = match_str(data_labels, layout.label);

    % Plot the data for the selected condition.
    chanWidth=layout.width(sellay)*2;
    chanHeight=layout.height(sellay)*2;
    chanX=layout.pos(sellay,1);
    chanY=layout.pos(sellay,2);
    chanLabel  = layout.label(sellay);
    xlim=[-20 20];
    ylim=[-0.5 1];
    xmin=xlim(1);
    xmax=xlim(2);
    ymin=ylim(1);
    ymax=ylim(2);
    h{f}=figure; f=f+1;
    for c=1:2:length(selchan)
        xval=xval_avgsubjects_allconditions(:, :, con);
        yval=yval_avgsubjects_allconditions(:, :, selchan(c), con);
        hold on;
        ft_plot_vector(xval, yval, 'width', chanWidth(c), 'height', chanHeight(c), 'hpos', chanX(c), 'vpos', chanY(c), 'hlim', [xmin xmax], 'vlim', [ymin ymax], 'color', 'rb','linewidth', 1, 'axis', 'yes');
    end
    title(char(taskname{con}))
    
    % Plot the layout, labels and outline.
    ft_plot_layout(layout, 'box', 'no', 'label', 'yes', 'outline', 'yes',...
        'point', 'no', 'mask', 'yes', 'fontsize', 8, 'labelyoffset',...
        1.4*median(layout.height/2), 'labelalignh', 'center', 'chanindx',...
        find(~ismember(layout.label, {'COMNT', 'SCALE'})), 'interpreter',...
        'none');

    % Plot scale.   
    l = find(strcmp(layout.label, 'SCALE'));
    x = layout.pos(l,1);
    y = layout.pos(l,2);
    plotScales([xmin xmax], [ymin ymax], x, y, chanWidth(1), chanHeight(1));
    
    % Save figure.
    set(gcf, 'Position', get(0, 'Screensize'));
    savetitle = strcat(char(taskname{con}), '_hemodynamicresponse');
    savetitle = savetitle(find(~isspace(savetitle)));
    savetitle = lower(savetitle);
    saveas(gcf, fullfile(results_path, savetitle),'png');
    
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