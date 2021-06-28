%% Analysis of the EEG signal (ERP: topoplots)
clear; clc; close all;

%% Initialize data
addpath('C:\Users\catar\OneDrive - Universidade do Porto\Twente\Combined-EEG-fNIRS-system\Analysis');
addpath('C:\Users\catar\OneDrive - Universidade do Porto\Twente\Combined-EEG-fNIRS-system\Analysis\Dual Task');
laptop = 'laptopCatarina';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);
results_path = 'C:\Users\catar\OneDrive - Universidade do Porto\Twente\Data Analysis\eeg\erp';

eeglab;

subrec = ["28" "04"];

% List of the 30 channels in the cap
list_channels = ["Fp1"; "Fpz"; "Fp2"; "F7"; "F3"; "AFFz"; "F4"; "F8";...
    "FC5"; "FC1"; "FC2"; "FC6"; "T7"; "C3"; "Cz"; "C4"; "T8"; "CP5";...
    "CP1"; "CP2"; "CP6"; "P7"; "P3"; "Pz"; "P4"; "P8"; "POz"; "O1";...
    "Oz"; "O2"];

%% Load data + processing per subject
% Go through all subjects
for subject = 1:size(subrec, 1)
    sub = subrec(subject, 1);
    rec = subrec(subject, 2);
  
    % Set settings
    checkBadChan = false;
    runBadTrial = true;
    
    % Load EEG preprocessed signals
    load([mainpath_in, '\pre-processed\sub-', char(sub), '\eeg\sub-',...
        char(sub), '_rec-', char(rec), '_eeg_divided.mat'], 'EEG_divided');
    
    % Separate into the four different tasks
    EEG_DualUncued = EEG_divided.EEG_NonAutoDualNoCue;
    EEG_SingleUncued = EEG_divided.EEG_NonAutoNoCue;
    EEG_DualCued = EEG_divided.EEG_NonAutoDualCue;
    EEG_SingleCued = EEG_divided.EEG_NonAutoCue;
    
    % Create a new dataset in eeglab
    % Dual Uncued
    [ALLEEG,EEG_DU,~] = pop_newset(ALLEEG, EEG_DualUncued, 1,'setname','eegdata','gui','off');
    EEG_DU.data = double(EEG_DU.data);
    clear EEG_DualUncued; 
    % Single Uncued
    [ALLEEG,EEG_SU,~] = pop_newset(ALLEEG, EEG_SingleUncued, 1,'setname','eegdata','gui','off');
    EEG_SU.data = double(EEG_SU.data);
    clear EEG_SingleUncued; 
    % Dual Cued
    [ALLEEG,EEG_DC,~] = pop_newset(ALLEEG, EEG_DualCued, 1,'setname','eegdata','gui','off');
    EEG_DC.data = double(EEG_DC.data);
    clear EEG_DualCued; 
    % Single Cued
    [ALLEEG,EEG_SC,~] = pop_newset(ALLEEG, EEG_SingleCued, 1,'setname','eegdata','gui','off');
    EEG_SC.data = double(EEG_SC.data);
    clear EEG_SingleCued; 
    
    EEG.conditions = {EEG_DU EEG_SU EEG_DC EEG_SC};
      
    for i = 1:length(EEG.conditions) % for each condition
    %% Remove bad channels 
    % Before filtering the data, bad channels will be removed if amplitude
    % is > 50mV (noise) or < 1mV (flat line)
    [~,badchans] = removeBadChansAmp(EEG.conditions{1,i}, checkBadChan);
    % Save bad channels in struct and remove channels from data
    EEG.conditions{1,i}.badchan = badchans;
    EEG.conditions{1,i}.badchan = [EEG.conditions{1,i}.badchan; badchans'];
    if ~isempty(badchans)
        [EEG.conditions{1,i}] = pop_select(EEG.conditions{1,i}, 'nochannel', badchans);
    end
    clear badchans
    
    %% Extract trials
    % Create data-matrix with trials [channel x time x trial]
    % Create time vector [sec]
    [data, time] = extractEpochs(EEG.conditions{1,i});
    EEG.conditions{1,i}.datatrials = data;
    EEG.conditions{1,i}.timetrials = time;
%     % Check for correct number of trials (267)
%     if size(EEG.conditions{1,i}.datatrials,3)~=267
%         error('INCORRECT NUMBER OF TRIALS: %d',size(EEG.conditions{1,i}.datatrials,3));
%     else
%         EEG.conditions{1,i}.datatrials(:,:,[268:size(EEG.conditions{1,i}.datatrials,3)]) = [];  
%     end
    
    %% Remove bad channels within condition
    % Remove bad channels if:
    % amplitude is higher than mean+2xSD
    % power is higher than mean+2xSD for all 3/5 frequency bands
    % (delta/theta/alpha/beta/gamma)
    % Repeat operation (checkforBadChans) and recalculate threshold (mean+2xSD)
    % as long as bad channels are found (large outlier can influence threshold
    % and as a consequence other bad channels will not be removed)
    [badchans] = badChansAmpPow(EEG.conditions{1,i}.datatrials, EEG.conditions{1,i}, checkBadChan);
    EEG.conditions{1,i}.badchan = badchans;
    %EEG.conditions{1,i}.badchan = [EEG.conditions{1,i}.badchan; badchans'];
    if ~isempty(badchans)
        for j = 1:size(badchans,1)
            idx(j,1) = find(strcmp({EEG.conditions{1,i}.chanlocs.labels},badchans(j)));
        end
        EEG.conditions{1,i}.datatrials(idx,:,:) = [];
        [EEG.conditions{1,i}] = pop_select(EEG.conditions{1,i}, 'nochannel', badchans);
    end
    clear badchans idx EEG.conditions{1,i}.datatrials
    
    %% Frequencies of Interest (FOI): 1-30Hz
    % FOI for ERP is 1-30 Hz using bandpass filter 
    EEG.conditions{1,i}.data = double(EEG.conditions{1,i}.data);
    [data_filt] = bandpassFilter(EEG.conditions{1,i}.data, EEG.conditions{1,i}.srate, [1 30]);
    % Save data in EEG structure and visualize raw signal
    EEG.conditions{1,i}.data = data_filt;   clear data_filt;
    if checkBadChan == true
        pop_eegplot(EEG.conditions{1,i},1,0,0);
    end
    
    %% Extract filtered (FOI) trials
    % Create data-matrix with trials [channel x time x trial]
    % (data is filtered for frequencies of interest)
    [EEG.conditions{1,i}.datatrials, ~] = extractEpochs(EEG.conditions{1,i});
%     % Check for correct number of trials (267)
%     if size(EEG.conditions{1,i}.datatrials,3)~=267
%         error('INCORRECT NUMBER OF TRIALS: %d',size(EEG.conditions{1,i}.datatrials,3));
%     else
%         EEG.conditions{1,i}.datatrials(:,:,[268:size(EEG.conditions{1,i}.datatrials,3)]) = []; 
%     end
    
    %% Baseline correction
    % Extract average of channel's baseline from channel's signal per trial
    [EEG.conditions{1,i}.datatrials] = baselineCorrection(EEG.conditions{1,i}.datatrials, EEG.conditions{1,i}.timetrials);
    
    %% Remove bad channels from trials
    % Remove bad channels from trial if:
    % amplitude is higher than mean+2xSD of channel
    % power is higher than mean+2xSD of channel in 2/4 frequency bands
    % (delta/theta/alpha/beta)
    % Replace bad channels within trial with NaN
    if runBadTrial == true
        [EEG.conditions{1,i}.badtrial] = badChansTrial(EEG.conditions{1,i}.datatrials, EEG.conditions{1,i});
        datanan = EEG.conditions{1,i}.datatrials;
        for iTrial = 1:size(datanan,3)
            for iChan = 1:size(datanan,1)
                if EEG.conditions{1,i}.badtrial(iChan,iTrial) == 1
                    datanan(iChan,:,iTrial) = NaN;
                end
            end
        end
        EEG.conditions{1,i}.datatrials = datanan; clear datanan;
    end
    
    %% Average ERP over trials
    % Average over trials, use omit NaNs (bad channels within trials)
    ERP = mean(EEG.conditions{1,i}.datatrials,3,'omitnan');
    
    %% Define structure for ERP
    EEG.conditions{1,i}.ERP.data = EEG.conditions{1,i}.datatrials;
    EEG.conditions{1,i}.ERP.ERP = ERP;
    EEG.conditions{1,i}.ERP.time = EEG.conditions{1,i}.timetrials;
    EEG.conditions{1,i}.ERP.chans = {EEG.conditions{1,i}.chanlocs.labels};
    EEG.conditions{1,i}.ERPBAD.chans = {EEG.conditions{1,i}.chanlocs.labels};
    EEG.conditions{1,i}.ERPBAD.badchans = EEG.conditions{1,i}.badchan;
    EEG.conditions{1,i}.ERPBAD.badtrial = EEG.conditions{1,i}.badtrial; 
    
    end
    
    % Save the values onto a allSubjects variable
    dualuncued_erp_allSubjects(:, subject) = EEG.conditions{1,1};
    singleuncued_erp_allSubjects(:, subject) = EEG.conditions{1,2};
    dualcued_erp_allSubjects(:, subject) = EEG.conditions{1,3};
    singlecued_erp_allSubjects(:, subject) = EEG.conditions{1,4};
    
    disp(['ERP processing done for subject ', char(sub), '.']);
    disp('Press any key to move onto the next subject.');
    pause;
    close all; clc; 
  
end

%% Topoplot of averaged ERP for all subjects per condition
% Average ERPS of all subjects
dualuncued_avgerp_allSubjects = mean(dualuncued_erp_allSubjects.ERP.ERP,1, 'omitnan');
singleuncued_avgerp_allSubjects = mean(singleuncued_erp_allSubjects.ERP.ERP,1, 'omitnan');
dualcued_avgerp_allSubjects = mean(dualcued_erp_allSubjects.ERP.ERP,1, 'omitnan');
singlecued_avgerp_allSubjects = mean(singlecued_erp_allSubjects.ERP.ERP,1, 'omitnan');

figure; 
axh(1) = subplot(2,2,1);     topoplot(dualuncued_avgerp_allSubjects, EEG_DU.chanlocs,'electrodes','off');   colorbar; caxlim(1,:) = caxis; set(gca,'FontSize',9);
axh(2) = subplot(2,2,2);     topoplot(singleuncued_avgerp_allSubjects, EEG_SU.chanlocs,'electrodes','off');   colorbar; caxlim(2,:) = caxis; set(gca,'FontSize',9);
axh(3) = subplot(2,2,3);     topoplot(dualcued_avgerp_allSubjects, EEG_DC.chanlocs,'electrodes','off');   colorbar; caxlim(3,:) = caxis; set(gca,'FontSize',9);
axh(4) = subplot(2,2,4);     topoplot(singlecued_avgerp_allSubjects, EEG_SC.chanlocs,'electrodes','off');   colorbar; caxlim(4,:) = caxis; set(gca,'FontSize',9);
  
set(axh,'clim',[-max(caxlim(:,2)) max(caxlim(:,2))]);
set(gcf,'Position',[590 470 470 320]);

saveas(gcf,fullfile([results_path,'ALLERP_topo.jpg']));

figure; 
pop_timtopo(dualuncued_erp_allSubjects);
pop_timtopo(singleuncued_erp_allSubjects);
pop_timtopo(dualcued_erp_allSubjects);
pop_timtopo(singlecued_erp_allSubjects);