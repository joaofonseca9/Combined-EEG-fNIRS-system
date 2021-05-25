%% addFolders.m
% set paths to all necessary folders for running anlaysis.
% set mother_folder (main folder for analysis and data)
% 
% OUTPUT:
% - folder = struct with all stings to pathnames

function [mainpath_in, mainpath_out]=addFolders(laptop)

    if strcmp(laptop, 'laptopJoao')
        % why?
        % add folder for analyses / data
        addpath(genpath('C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\EEG scripts (introductory)\dryEEG_analysisscripts_Janne'));
    %     cd 'C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\EEG scripts (introductory)\dryEEG_analysisscripts_Janne';

        % add EEGlab
        addpath('C:\Users\joaop\Downloads\eeglab2021.0');
        
        % add fieldtrip
        addpath('C:\Users\joaop\Downloads\fieldtrip-20210212\fieldtrip-20210212');
        
        % add the path of the analysis scripts
        % por colocar
        
        % set paths for the input and output data
        mainpath_in='C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\Experiment\Data\Pilots';
        mainpath_out='C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\Experiment\Data\Pilots\pre-processed\';
        
    elseif strcmp(laptop, 'laptopMariana')
        % add folder for analyses / data
    %    addpath(genpath('C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\EEG scripts (introductory)\dryEEG_analysisscripts_Janne'));
    %     cd 'C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\EEG scripts (introductory)\dryEEG_analysisscripts_Janne';

        % add EEGlab
        addpath('C:\Users\maria\OneDrive\Ambiente de Trabalho\eeglab2021.0');
        
        % add fieldtrip
        addpath('C:\Users\maria\OneDrive\Documentos\GitHub\fieldtrip');
        
        % add the path of the analysis scripts
        addpath('C:\Users\maria\OneDrive\Documentos\GitHub\Combined-EEG-fNIRS-system\Analysis');
        
        % set paths for the input and output data
        mainpath_in='C:\Users\maria\Universidade do Porto\João Pedro Barbosa Fonseca - Internship\Experiment\Data\Pilots';
        mainpath_out='C:\Users\maria\Universidade do Porto\João Pedro Barbosa Fonseca - Internship\Experiment\Data\Pilots\pre-processed\';

    elseif strcmp(laptop, 'laptopCatarina')
        % why?
        % add folder for analyses / data
        % addpath(genpath('C:\Users\catar\OneDrive - Universidade do Porto\Internship\After Experiment\Previous Analysis Scripts\dryEEG_analysisscripts_Janne'));
        % cd 'C:\Users\catar\OneDrive - Universidade do Porto\Internship\After Experiment\Previous Analysis Scripts\dryEEG_analysisscripts_Janne';

        % add EEGlab
        addpath('C:\Program Files\Matlab\R2020b\toolbox\eeglab_current\eeglab2021.0');
        
        % add fieldtrip
        addpath('C:\Program Files\Matlab\R2020b\toolbox\fieldtrip-20210308');
        
        % add the path of the analysis scripts
        addpath('C:\Users\catar\OneDrive - Universidade do Porto\Twente\Combined-EEG-fNIRS-system\Analysis');
        
        % set paths for the input and output data
        mainpath_in='C:\Users\catar\OneDrive - Universidade do Porto\Internship\Experiment\Data\Pilots';
        mainpath_out='C:\Users\catar\OneDrive - Universidade do Porto\Internship\Experiment\Data\Pilots\pre-processed';
        
    end
end

