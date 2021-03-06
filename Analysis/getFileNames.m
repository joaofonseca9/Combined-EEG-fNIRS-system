% Get the filenames to save/load preprocessed data.
%
% @file getFileNames.m
%
% From the main path, the subject id and the number of the recording, get
% all the necessary paths to save the different steps of preprocessing.

function [file] = getFileNames(mainpath, ID, rec)

file.filtered = fullfile([mainpath, '\sub-', ID, '\eeg\sub-', ID, '_rec-',...
    rec, '_eeg_filtered.mat']);
file.removedBadChannels = fullfile([mainpath, '\sub-', ID, '\eeg\sub-', ID,...
    '_rec-', rec, '_eeg_removedbadchannels.mat']);
file.preICA = fullfile([mainpath, '\sub-', ID, '\eeg\sub-', ID, '_rec-',...
    rec, '_eeg_preica.mat']);
file.pstICA = fullfile([mainpath, '\sub-', ID, '\eeg\sub-', ID, '_rec-',...
    rec, '_eeg_pstica.mat']);
file.preprocessed = fullfile([mainpath, '\sub-', ID, '\eeg\sub-', ID, '_rec-',...
    rec, '_eeg_preprocessed.mat']);
file.LETTER = fullfile([mainpath, '\sub-', ID, '\eeg\sub-',...
    ID, '_rec-', rec, '_eeg_letter.mat']);
file.EEG_divided = fullfile([mainpath, '\sub-', ID, '\eeg\sub-',...
    ID, '_rec-', rec, '_eeg_divided.mat']);

end