%% Comparison of the performance under the dual-tasking condition

clear; clc; close all;
addpath('C:\Users\maria\OneDrive\Documentos\GitHub\Combined-EEG-fNIRS-system\Analysis');

laptop = 'laptopMariana';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);

sub = '02';
rec = '02';

% Load the results
load([mainpath_in, '\source\sub-', sub, '\stim\results_sub-', sub,...
    '_rec-', rec, '.mat']);

%% Automatic sequence
% Check if answers from automatic sequence were correct
autodual = events_autodual.trial;
autodual_correct = zeros(1, length(autodual));

for trial = 1:length(autodual)
    stimuli = autodual(trial).stimuli;
    if count(cell2mat(stimuli.value),'G') == str2num(cell2mat(stimuli.response))
        autodual_correct(trial) = 1;
    end
end

% Check average number of mistakes per condition (cued and uncued)
autodual_mistakescued = 0;
autodual_mistakesuncued = 0;

for trial = 1:length(autodual)
    if autodual_correct(trial) == 0 && autodual(trial).cue == 1
        autodual_mistakescued = autodual_mistakescued + 1;
    elseif autodual_correct(trial) == 0 && autodual(trial).cue == 0
        autodual_mistakesuncued = autodual_mistakesuncued + 1;
    end
end
autodual_avgmistakescued = autodual_mistakescued/(length(autodual)/2);
autodual_avgmistakesuncued = autodual_mistakesuncued/(length(autodual)/2);

%% Non-automatic sequence
% Check if answers from non-automatic sequence were correct.
nonautodual = events_nonautodual.trial;
nonautodual_correct = zeros(1, length(nonautodual));

for trial = 1:length(nonautodual)
    stimuli = nonautodual(trial).stimuli;
    if count(cell2mat(stimuli.value),'G') == str2num(cell2mat(stimuli.response))
        nonautodual_correct(trial) = 1;
    end
end

% Check average number of mistakes per condition (cued and uncued)
nonautodual_mistakescued = 0;
nonautodual_mistakesuncued = 0;

for trial = 1:length(nonautodual)
    if nonautodual_correct(trial) == 0 && nonautodual(trial).cue == 1
        nonautodual_mistakescued = nonautodual_mistakescued + 1;
    elseif nonautodual_correct(trial) == 0 && nonautodual(trial).cue == 0
        nonautodual_mistakesuncued = nonautodual_mistakesuncued + 1;
    end
end
nonautodual_avgmistakescued = nonautodual_mistakescued/(length(nonautodual)/2);
nonautodual_avgmistakesuncued = nonautodual_mistakesuncued/(length(nonautodual)/2);
