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
    
    taskname = {'Auto Cued', 'Non-Auto Cued', 'Auto Uncued', 'Non-Auto Uncued'};
    h = multiplot_condition(nirs, layout, [2 4 6 8], taskname,...
        'baseline', [-10 0], 'trials', false, 'topoplot', 'no', 'ylim',...
        [-0.5 1]);

%   tasktimelock = {['sub-',sub,'_rec-',rec,'_autodualcued_timelock'], ['sub-',sub,'_rec-',rec,'_autosinglecued_timelock'], ['sub-',sub,'_rec-',rec,'_nonautodualcued_timelock'], ['sub-',sub,'_rec-',rec,'_nonautosinglecued_timelock'],['sub-',sub,'_rec-',rec,'_autodualuncued_timelock'], ['sub-',sub,'_rec-',rec,'_autosingleuncued_timelock'], ['sub-',sub,'_rec-',rec,'_nonautodualuncued_timelock'], ['sub-',sub,'_rec-',rec,'_nonautosingleuncued_timelock']};
%   taskspatialO2Hb = {['sub-',sub,'_rec-',rec,'_autodualcued_spatialO2Hb'], ['sub-',sub,'_rec-',rec,'_autosinglecued_spatialO2Hb'], ['sub-',sub,'_rec-',rec,'_nonautodualcued_spatialO2Hb'], ['sub-',sub,'_rec-',rec,'_nonautosinglecued_spatialO2Hb'],['sub-',sub,'_rec-',rec,'_autodualuncued_spatialO2Hb'], ['sub-',sub,'_rec-',rec,'_autosingleuncued_spatialO2Hb'], ['sub-',sub,'_rec-',rec,'_nonautodualuncued_spatialO2Hb'], ['sub-',sub,'_rec-',rec,'_nonautosingleuncued_spatialO2Hb']};
%   taskspatialHHb = {['sub-',sub,'_rec-',rec,'_autodualcued_spatialHHb'], ['sub-',sub,'_rec-',rec,'_autosinglecued_spatialHHb'], ['sub-',sub,'_rec-',rec,'_nonautodualcued_spatialHHb'], ['sub-',sub,'_rec-',rec,'_nonautosinglecued_spatialHHb'],['sub-',sub,'_rec-',rec,'_autodualuncued_spatialHHb'], ['sub-',sub,'_rec-',rec,'_autosingleuncued_spatialHHb'], ['sub-',sub,'_rec-',rec,'_nonautodualuncued_spatialHHb'], ['sub-',sub,'_rec-',rec,'_nonautosingleuncued_spatialHHb']};
    
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