%% Statistical Analysis: EEG (ERP)
clear; clc; close all;

%% Initialize data
addpath('C:\Users\catar\OneDrive - Universidade do Porto\Twente\Combined-EEG-fNIRS-system\Analysis');
addpath('C:\Users\catar\OneDrive - Universidade do Porto\Twente\Combined-EEG-fNIRS-system\Analysis\Dual Task');
laptop = 'laptopCatarina';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);
results_path = 'C:\Users\catar\OneDrive - Universidade do Porto\Twente\Data Analysis\eeg\erp';
statistics_path = 'C:\Users\catar\OneDrive - Universidade do Porto\Twente\Data Analysis\statistics\eeg\erp';

subrec = ["28" "04"; "02" "02"; "76" "01"; "64" "01"];

% List of the 30 channels in the cap
list_channels = ["Fp1"; "Fpz"; "Fp2"; "F7"; "F3"; "AFFz"; "F4"; "F8";...
    "FC5"; "FC1"; "FC2"; "FC6"; "T7"; "C3"; "Cz"; "C4"; "T8"; "CP5";...
    "CP1"; "CP2"; "CP6"; "P7"; "P3"; "Pz"; "P4"; "P8"; "POz"; "O1";...
    "Oz"; "O2"];

for subject=1:size(subrec, 1)
    sub = subrec(subject, 1);
    rec = subrec(subject, 2);
    
    %% Load and collect subject data
    % Load data
    load([results_path,'\sub-',char(sub),'\ERP_DualCued.mat'], 'ERPdual'); dualcued = ERPdual; clear ERPdual;
    load([results_path,'\sub-',char(sub),'\ERP_DualUncued.mat'], 'ERPdual'); dualuncued = ERPdual; clear ERPdual;
    load([results_path,'\sub-',char(sub),'\ERP_SingleCued.mat'], 'ERPdual'); singlecued = ERPdual; clear ERPdual;
    load([results_path,'\sub-',char(sub),'\ERP_SingleUncued.mat'], 'ERPdual'); singleuncued = ERPdual; clear ERPdual;
    
    % Collect ERP and GFPt from all subjects
    for iChan = 1:length(list_channels)
        dualcuedloc = find(strcmp(list_channels(iChan), dualcued{1, subject}.chans));
        if ~isempty(dualcuedloc);        allERP.dualcued(iChan,:,subject) = dualcued{1, subject}.ERP(dualcuedloc,:);
        elseif isempty(dualcuedloc);     allERP.dualcued(iChan,:,subject) = NaN;
        end
        singlecuedloc = find(strcmp(list_channels(iChan), singlecued{1, subject}.chans));
        if ~isempty(singlecuedloc);        allERP.singlecued(iChan,:,subject) = singlecued{1, subject}.ERP(singlecuedloc,:);
        elseif isempty(singlecuedloc);     allERP.singlecued(iChan,:,subject) = NaN;
        end
        singleuncuedloc = find(strcmp(list_channels(iChan), singleuncued{1, subject}.chans));
        if ~isempty(singleuncuedloc);        allERP.singleuncued(iChan,:,subject) = singleuncued{1, subject}.ERP(singleuncuedloc,:);
        elseif isempty(singleuncuedloc);     allERP.singleuncued(iChan,:,subject) = NaN;
        end
        dualuncuedloc = find(strcmp(list_channels(iChan), dualuncued{1, subject}.chans));
        if ~isempty(dualuncuedloc);        allERP.dualuncued(iChan,:,subject) = dualuncued{1, subject}.ERP(dualuncuedloc,:);
        elseif isempty(dualuncuedloc);     allERP.dualuncued(iChan,:,subject) = NaN;
        end
        clear dualcuedloc singlecuedloc dualuncuedloc singleuncuedloc
    end
    allGFP.dualcued(subject,:) = dualcued{1, subject}.GFPt;
    allGFP.singlecued(subject,:) = singlecued{1, subject}.GFPt;
    allGFP.dualuncued(subject,:) = dualuncued{1, subject}.GFPt;
    allGFP.singleuncued(subject,:) = singleuncued{1, subject}.GFPt;
      
    %% Spearman's Rank Correlation (COR) - Fiedler (2015)
    % Dual
    for iL = 1:length(dualcued{1, subject}.GFPt)
        Dg(iL) = dualuncued{1, subject}.GFPt(iL)-mean(dualuncued{1, subject}.GFPt);
        Dd(iL) = dualcued{1, subject}.GFPt(iL)-mean(dualcued{1, subject}.GFPt);
    end
    cor_dual(subject,1) = sum(Dg.*Dd) / sqrt( sum(Dg.^2).*sum(Dd.^2) ); 
        
    % Uncued
    for iL = 1:length(dualuncued{1, subject}.GFPt)
        Dg(iL) = singleuncued{1, subject}.GFPt(iL)-mean(singleuncued{1, subject}.GFPt);
        Dd(iL) = dualuncued{1, subject}.GFPt(iL)-mean(dualuncued{1, subject}.GFPt);
    end
    cor_uncued(subject,1) = sum(Dg.*Dd) / sqrt( sum(Dg.^2).*sum(Dd.^2) ); 
    
    % Cued
    for iL = 1:length(dualcued{1, subject}.GFPt)
        Dg(iL) = singlecued{1, subject}.GFPt(iL)-mean(singlecued{1, subject}.GFPt);
        Dd(iL) = dualcued{1, subject}.GFPt(iL)-mean(dualcued{1, subject}.GFPt);
    end
    cor_cued(subject,1) = sum(Dg.*Dd) / sqrt( sum(Dg.^2).*sum(Dd.^2) ); 
    
    clear D Dg Dd dualcued dualuncued
    
end

%% Average Spearman's Rank Correlation
stats.cor_dual(1,1) = mean(cor_dual);
stats.cor_dual(1,2) = 2*std(cor_dual);
stats.cor_uncued(1,1) = mean(cor_uncued);
stats.cor_uncued(1,2) = 2*std(cor_uncued);
stats.cor_cued(1,1) = mean(cor_cued);
stats.cor_cued(1,2) = 2*std(cor_cued);

disp(stats.cor_dual);
disp(stats.cor_uncued);
disp(stats.cor_cued);

save(strcat(statistics_path, '\stats_erp.mat'), 'stats');
