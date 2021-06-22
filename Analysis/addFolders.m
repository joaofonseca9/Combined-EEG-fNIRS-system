%% addFolders.m
% set paths to all necessary folders for running anlaysis.
% set mother_folder (main folder for analysis and data)
% 
% OUTPUT:
% - folder = struct with all stings to pathnames

function [mainpath_in, mainpath_out, eeglab_path]=addFolders(laptop)

    if strcmp(laptop, 'laptopJoao')
        % add EEGlab
        eeglab_path = 'C:\\Users\\joaop\\Downloads\\eeglab2021.0';
        addpath(eeglab_path);
        
        % add fieldtrip
        addpath('C:\Users\joaop\Downloads\fieldtrip-20210212\fieldtrip-20210212');
        
        % set paths for the input and output data
        mainpath_in='C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\Experiment\Data\Exp';
        mainpath_out='C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\Experiment\Data\Exp\pre-processed\';
        
    elseif strcmp(laptop, 'laptopMariana')
        % add EEGlab
        eeglab_path = 'C:\\Users\\maria\\OneDrive\\Ambiente de Trabalho\\eeglab2021.0';
        addpath(eeglab_path);
        
        % add fieldtrip
        addpath('C:\Users\maria\OneDrive\Documentos\GitHub\fieldtrip');
        
        % set paths for the input and output data
        mainpath_in='C:\Users\maria\Universidade do Porto\João Pedro Barbosa Fonseca - Internship\Experiment\Data\Pilots';
        mainpath_out='C:\Users\maria\Universidade do Porto\João Pedro Barbosa Fonseca - Internship\Experiment\Data\Pilots\pre-processed\';

    elseif strcmp(laptop, 'laptopCatarina')
        % add EEGlab
        eeglab_path = 'C:\\Program Files\\Matlab\\R2020b\\toolbox\\eeglab_current\\eeglab2021.0';
        addpath(eeglab_path);
        
        % add fieldtrip
        addpath('C:\Program Files\Matlab\R2020b\toolbox\fieldtrip-20210308');
        
        % set paths for the input and output data
        mainpath_in='C:\Users\catar\OneDrive - Universidade do Porto\Internship\Experiment\Data\Exp';
        mainpath_out='C:\Users\catar\OneDrive - Universidade do Porto\Internship\Experiment\Data\Exp\pre-processed';
        
    end
end

