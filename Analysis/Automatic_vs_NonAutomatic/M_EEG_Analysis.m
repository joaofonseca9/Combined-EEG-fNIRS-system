clear; clc; close all;
addpath('C:\Users\maria\OneDrive\Documentos\GitHub\Combined-EEG-fNIRS-system\Analysis');

laptop = 'laptopMariana';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);
eeglab;
ft_defaults;

sub='04';
rec='01';

file = getFileNames(mainpath_out, sub, rec);
%%
load(file.EEG_divided, 'EEG_divided');

EEG_AutoUncued = EEG_divided.EEG_AutoNoCue;
EEG_NonAutoUncued = EEG_divided.EEG_NonAutoNoCue;
%%
% Start_AutoCue  = event_samp((strcmp({EEG_AutoCue.event.type}, sprintf('s%d',marker_table.StartAutoCue(1))))==1);
% Stop_AutoCue   = event_samp((strcmp({EEG.event.type}, sprintf('s%d',marker_table.StopAutoCue(1))))==1);

% AFFz - SMA
% figure;
% topoplot(EEG_AutoUncued.data(6, :), EEG_AutoUncued.chanlocs);
% figure;
% topoplot(EEG_NonAutoUncued.data(6, :), EEG_NonAutoUncued.chanlocs);

% F7 - DLPFC
% figure;
% topoplot(EEG_AutoUncued.data(4, :), EEG_AutoUncued.chanlocs);
% figure;
% topoplot(EEG_NonAutoUncued.data(4, :), EEG_NonAutoUncued.chanlocs);

% pop_spectopo(EEG_AutoUncued, 1);

%% Topoplots for Auto Uncued
event_samp  = [EEG_AutoUncued.event.latency];

startAutoUncued = find(strcmp({EEG_AutoUncued.event.type}, 's1703')==1);
endAutoUncued = find(strcmp({EEG_AutoUncued.event.type}, 's1711')==1);
keypresses = event_samp((strcmp({EEG_AutoUncued.event.type}, 's1777')==1));

for trial=1:10
    title = char(strcat('Trial_', string(trial)));
    startAutoUncued_times = event_samp(startAutoUncued(trial));
    endAutoUncued_times = event_samp(endAutoUncued(trial));
    keypresses_times = keypresses((keypresses>startAutoUncued_times & keypresses<endAutoUncued_times));
    keypresses_times = keypresses_times - startAutoUncued_times;
    EEG_trial = pop_select(EEG_AutoUncued, 'point',...
        [startAutoUncued_times endAutoUncued_times]);
    pop_topoplot(EEG_trial, 1, keypresses_times, title); 
    
    % EEG_trial_key = pop_select(EEG_trial, 'point',...
        [keypresses_times(1)-2 keypresses_times(1)+2]);
    % figure; pop_plottopo(EEG_trial_key, [1:30] , 'preprocessed', 0, 'ydir',1);
end

%% Topoplots for Non-Auto Uncued
event_samp  = [EEG_NonAutoUncued.event.latency];

startNonAutoUncued = find(strcmp({EEG_NonAutoUncued.event.type}, 's1705')==1);
endNonAutoUncued = find(strcmp({EEG_NonAutoUncued.event.type}, 's1713')==1);
keypresses = event_samp((strcmp({EEG_NonAutoUncued.event.type}, 's1777')==1));

for trial=1:1
    title = char(strcat('Trial_', string(trial)));
    startNonAutoUncued_times = event_samp(startNonAutoUncued(trial));
    endNonAutoUncued_times = event_samp(endNonAutoUncued(trial));
    keypresses_times = keypresses((keypresses>startNonAutoUncued_times & keypresses<endNonAutoUncued_times));
    keypresses_times = keypresses_times - startNonAutoUncued_times;
    EEG_trial = pop_select(EEG_NonAutoUncued, 'point',...
        [startNonAutoUncued_times endNonAutoUncued_times]);
    pop_topoplot(EEG_trial, 1, keypresses_times, title);    
end