%% addFolders.m
% set paths to all necessary folders for running anlaysis.
% set mother_folder (main folder for analysis and data)
% 
% OUTPUT:
% - folder = struct with all stings to pathnames

if strcmp(laptop, 'laptopJoao')
    % add folder for analyses / data
    addpath(genpath('C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\EEG scripts (introductory)\dryEEG_analysisscripts_Janne'));
%     cd 'C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\EEG scripts (introductory)\dryEEG_analysisscripts_Janne';
    
    % add EEGlab
    addpath('C:\Users\joaop\Downloads\eeglab2021.0');
    % add fieldtrip
    addpath('C:\Users\joaop\Downloads\fieldtrip-20210212\fieldtrip-20210212');
    
elseif strcmp(laptop, 'laptopMariana')
    % add folder for analyses / data
%    addpath(genpath('C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\EEG scripts (introductory)\dryEEG_analysisscripts_Janne'));
%     cd 'C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\EEG scripts (introductory)\dryEEG_analysisscripts_Janne';
    
    % add EEGlab
    addpath('C:\Users\maria\OneDrive\Ambiente de Trabalho\eeglab2021.0');
    % add fieldtrip
    addpath('C:\Users\maria\OneDrive\Documentos\GitHub\fieldtrip');
 
elseif strcmp(laptop, 'laptopCatarina')
    % add folder for analyses / data
    addpath(genpath('C:\Users\catar\OneDrive - Universidade do Porto\Internship\After Experiment\Previous Analysis Scripts\dryEEG_analysisscripts_Janne'));
    % cd 'C:\Users\catar\OneDrive - Universidade do Porto\Internship\After Experiment\Previous Analysis Scripts\dryEEG_analysisscripts_Janne';
    
    % add EEGlab
    addpath('C:\Program Files\Matlab\R2020b\toolbox\eeglab_current\eeglab2021.0');
    % add fieldtrip
    addpath('C:\Program Files\Matlab\R2020b\toolbox\fieldtrip-20210308');
        
end

