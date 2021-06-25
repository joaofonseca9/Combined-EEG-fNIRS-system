%% Analysis of the fNIRS signal (Hemodynamic Responses)
clear; clc; close all;

%% Initialize data
laptop = 'laptopCatarina';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);
results_path = 'C:\Users\catar\OneDrive - Universidade do Porto\Twente\Data Analysis\Hemodynamic Response';

subrec = ["28" "02"];

%% Load data + baseline correction (preprocessing), timelock analysis + spatial representation + statistical testing
% Go through all subjects
for subject = 1:size(subrec, 1)
    sub = subrec(subject, 1);
    rec = subrec(subject, 2);
  
    % Load fNIRS signals
    load([mainpath_in, '\pre-processed\sub-', char(sub), '\nirs\sub-',...
        char(sub), '_rec-', char(rec), '_nirs_preprocessed.mat']);
    load([mainpath_in, '\pre-processed\sub-', char(sub), '\nirs\sub-',...
        char(sub), '_rec-', char(rec), '_nirs_epoch.mat']);
    nirs_preprocessed.trialinfo = nirs_epoch.trialinfo;
    nirs_preprocessed.sampleinfo = nirs_epoch.sampleinfo;
    
    % Load layout 
    load(fullfile(mainpath_out,['sub-',char(sub)],'3d','layout.mat'), 'layout');
    
    % Keep the trials of interest (Dual Cued, Single Cued, Dual Uncued, Single
    % Uncued
    nirs = keepTrialsInterest(nirs_preprocessed);
    
    % Get the hemodynamic response of all conditions for the subject
    taskname = {'Dual Task Cued', 'Single Task Cued', 'Dual Task Uncued', 'Single Task Uncued'};
    [h, x_allconditions, y_allconditions] =...
        multiplot_condition(nirs, layout, [3 4 7 8], taskname,...
        'baseline', [-10 0], 'trials', false, 'topoplot', 'no', 'ylim',...
        [-0.5 1]);
    tasktimelock = {['sub-',char(sub),'_rec-',char(rec),'_dualcued_timelock'], ['sub-',char(sub),'_rec-',char(rec),'_singlecued_timelock'], ['sub-',char(sub),'_rec-',char(rec),'_dualuncued_timelock'], ['sub-',char(sub),'_rec-',char(rec),'_singleuncued_timelock']};
    taskspatialO2Hb = {['sub-',char(sub),'_rec-',char(rec),'_dualcued_spatialO2Hb'], ['sub-',char(sub),'_rec-',char(rec),'_singlecued_spatialO2Hb'], ['sub-',char(sub),'_rec-',char(rec),'_dualuncued_spatialO2Hb'], ['sub-',char(sub),'_rec-',char(rec),'_singleuncued_spatialO2Hb']};
    taskspatialHHb = {['sub-',char(sub),'_rec-',char(rec),'_dualcued_spatialHHb'], ['sub-',char(sub),'_rec-',char(rec),'_singlecued_spatialHHb'], ['sub-',char(sub),'_rec-',char(rec),'_dualuncued_spatialHHb'], ['sub-',char(sub),'_rec-',char(rec),'_singleuncued_spatialHHb']};
    
    % Save figures
    cd(results_path);
    if not(isfolder('timelock + spatial'))
       mkdir('timelock + spatial')
    end
    cd(fullfile(results_path, ['sub-',char(sub)], 'timelock + spatial'));
    saveas(h{1},tasktimelock{1},'png'); saveas(h{2},taskspatialO2Hb{1},'png'); saveas(h{3},taskspatialHHb{1},'png'); 
    saveas(h{4},tasktimelock{2},'png'); saveas(h{5},taskspatialO2Hb{2},'png'); saveas(h{6},taskspatialHHb{2},'png'); 
    saveas(h{7},tasktimelock{3},'png'); saveas(h{8},taskspatialO2Hb{3},'png'); saveas(h{9},taskspatialHHb{3},'png');
    saveas(h{10},tasktimelock{4},'png'); saveas(h{11},taskspatialO2Hb{4},'png'); saveas(h{12},taskspatialHHb{4},'png'); 

    % Statistical testing 
    cd(results_path);
    if not(isfolder('statistics'))
       mkdir('statistics')
    end
    cd(fullfile(results_path,['sub-',char(sub)],'statistics'));
    for i=1:4 % loop over the 4 conditions
      [stat_O2Hb, stat_HHb] = statistics_withinsubjects(nirs, 'nirs', layout, i, taskname{i}, char(sub), char(rec));
    end
       
    % Insert hemodynamic responses in a matrix of all subjects
    x_allsubjects_allconditions(:, :, :, subject) = x_allconditions; 
    y_allsubjects_allconditions(:, :, :, :, subject) = y_allconditions; 
    
    disp(['These are the results for subject ', char(sub), '.']);
    disp('Press any key to move onto the next subject.');
    pause;
    close all;
end

%% Average the hemodynamic responses over all subjects
x_avgsubjects_allconditions = mean(x_allsubjects_allconditions, 4); 
y_avgsubjects_allconditions = mean(y_allsubjects_allconditions, 5);

% For each condition
for condition=1:size(x_avgsubjects_allconditions, 3)
    f=1;

    % Select the channels that are both in data and layout
    data_labels=nirs.label;
    for i=1:length(data_labels)
        tmp = strsplit(data_labels{i});
        data_labels{i}=tmp{1};
    end
    % Take the subselection of channels that is contained
    [selchan, sellay] = match_str(data_labels, layout.label);

    % Plot the data for the selected condition
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
        xval=x_avgsubjects_allconditions(:, :, condition);
        yval=y_avgsubjects_allconditions(:, :, selchan(c), condition);
        hold on;
        ft_plot_vector(xval, yval, 'width', chanWidth(c), 'height', chanHeight(c), 'hpos', chanX(c), 'vpos', chanY(c), 'hlim', [xmin xmax], 'vlim', [ymin ymax], 'color', 'rb','linewidth', 1, 'axis', 'yes');
    end
    title(char(taskname{condition}))
    
    % Plot the layout, labels and outline
    ft_plot_layout(layout, 'box', 'no', 'label', 'yes', 'outline', 'yes',...
        'point', 'no', 'mask', 'yes', 'fontsize', 8, 'labelyoffset',...
        1.4*median(layout.height/2), 'labelalignh', 'center', 'chanindx',...
        find(~ismember(layout.label, {'COMNT', 'SCALE'})), 'interpreter',...
        'none');

    % Plot scale
    l = find(strcmp(layout.label, 'SCALE'));
    x = layout.pos(l,1);
    y = layout.pos(l,2);
    plotScales([xmin xmax], [ymin ymax], x, y, chanWidth(1), chanHeight(1));
    
    % Save figure
    set(gcf, 'Position', get(0, 'Screensize'));
    savetitle = strcat(char(taskname{condition}), '_avgsubs_hemodynamicresponse');
    savetitle = savetitle(find(~isspace(savetitle)));
    savetitle = lower(savetitle);
    saveas(gcf, fullfile(results_path, savetitle),'png');
    
end

disp('This was the end of individual subjects.');
disp('These are the results for the average of all subjects.');

