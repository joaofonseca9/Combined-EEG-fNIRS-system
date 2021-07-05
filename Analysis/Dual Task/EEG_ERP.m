%% Analysis of the EEG signal (ERP: topoplots)
clear; clc; close all;

%% Initialize data
addpath('C:\Users\catar\OneDrive - Universidade do Porto\Twente\Combined-EEG-fNIRS-system\Analysis');
addpath('C:\Users\catar\OneDrive - Universidade do Porto\Twente\Combined-EEG-fNIRS-system\Analysis\Dual Task');
laptop = 'laptopCatarina';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);
results_path = 'C:\Users\catar\OneDrive - Universidade do Porto\Twente\Data Analysis\eeg\erp';

eeglab;

subrec = ["28" "04";"64" "01"; "02" "02";"76" "01"];

% List of the 30 channels in the cap
list_channels = ["Fp1"; "Fpz"; "Fp2"; "F7"; "F3"; "AFFz"; "F4"; "F8";...
    "FC5"; "FC1"; "FC2"; "FC6"; "T7"; "C3"; "Cz"; "C4"; "T8"; "CP5";...
    "CP1"; "CP2"; "CP6"; "P7"; "P3"; "Pz"; "P4"; "P8"; "POz"; "O1";...
    "Oz"; "O2"];

taskname = {'Dual Uncued', 'Single Uncued', 'Dual Cued', 'Single Cued'};

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
    
    %% Re-reference to M2 and create a new dataset in eeglab
    % Dual Uncued
    locs_DU = {EEG_DualUncued.chanlocs.labels};
    M2_loc_DU = find(contains(locs_DU, 'M2'));
    [EEG_DU] = pop_reref(EEG_DualUncued, M2_loc_DU);
    [ALLEEG, EEG_DU, CURRENTSET] = pop_newset(ALLEEG, EEG_DU, 1,'setname','eegdata','gui','off');
    EEG_DU.data = double(EEG_DU.data);
    clear EEG_DualUncued; 
    
    % Single Uncued
    locs_SU = {EEG_SingleUncued.chanlocs.labels};
    M2_loc_SU = find(contains(locs_SU, 'M2'));
    [EEG_SU] = pop_reref(EEG_SingleUncued, M2_loc_SU);
    [ALLEEG, EEG_SU, CURRENTSET] = pop_newset(ALLEEG, EEG_SU, 1,'setname','eegdata','gui','off');
    EEG_SU.data = double(EEG_SU.data);
    clear EEG_SingleUncued; 
    
    % Dual Cued
    locs_DC = {EEG_DualCued.chanlocs.labels};
    M2_loc_DC = find(contains(locs_DC, 'M2'));
    [EEG_DC] = pop_reref(EEG_DualCued, M2_loc_DC);
    [ALLEEG, EEG_DC, CURRENTSET] = pop_newset(ALLEEG, EEG_DC, 1,'setname','eegdata','gui','off');
    EEG_DC.data = double(EEG_DC.data);
    clear EEG_DualCued; 
    
    % Single Cued
    locs_SC = {EEG_SingleCued.chanlocs.labels};
    M2_loc_SC = find(contains(locs_SC, 'M2'));
    [EEG_SC] = pop_reref(EEG_SingleCued, M2_loc_SC);
    [ALLEEG, EEG_SC, CURRENTSET] = pop_newset(ALLEEG, EEG_SC, 1,'setname','eegdata','gui','off');
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
    
    %% Global field power in time domain (GFPt)
    % calculate difference between two channels, per time point 
    % take the sum of the differences, per time point
    % gives the GFPt over time
    [GFPt]  = globalFieldPotential(ERP);
    
    %% Plotting
    figure
    subplot(1,2,1); plot(time,mean(ERP,3,'omitnan'),'b');
    subplot(1,2,2); hold on; h1 = fill([time,fliplr(time)], [GFPt,fliplr(GFPt)],'b','LineStyle','none');
    set(h1,'FaceAlpha',0.4); plot(time,GFPt,'r','LineWidth',1.5);
    
    subplot(1,2,1); set(gca,'FontSize',11); box on;
    ylabel('Potential (\muV)','FontSize',14); title('ERP','FontSize',14); ylim([-6 6])
    xticks([-0.1:0.1:0.4]); yticks([-5:5:5]); line([0 0],[-6 6],'Color','k');xlabel('Time (s)','FontSize',14);
    subplot(1,2,2); set(gca,'FontSize',11); box on;
    ylabel('GFPt (\muV)','FontSize',14); xlabel('Time (s)','FontSize',14); ylim([-2 26])
    xticks([-0.1:0.1:0.4]); yticks([0:10:20]); line([0 0],[-2 25],'Color','k');  title('GFPt');
    
    
    %% Save structure per subject per condition
    ERPdual{i}.data        = EEG.conditions{1,i}.datatrials;
    ERPdual{i}.ERP         = ERP;
    ERPdual{i}.GFPt        = GFPt;
    ERPdual{i}.time        = EEG.conditions{1,i}.timetrials;
    ERPdual{i}.chans       = {EEG.conditions{1,i}.chanlocs.labels};
    ERPdual{i}.condition = EEG.conditions{1,i};
       
    save(fullfile(results_path, ['sub-',char(sub)],['ERP_',taskname{1,i}]));
    
    end 
    
    disp(['ERP processing done for subject ', char(sub), '.']);
    disp('Press any key to move onto the next subject.');
    %pause;
    close all; clc; 
  
end
disp('This was the end of individual subjects.');

%% Collect ERP of all subjects 
for subject = 1:size(subrec, 1)
    sub = subrec(subject, 1);
    rec = subrec(subject, 2);
   
    for i = 1:length(EEG.conditions)
        % Collect ERP of all subjects [chan x time x subject]
        allERP{1,i}.ERPdual = zeros(length(ERPdual{1,i}.chans),length(ERPdual{1,i}.ERP));
    
        for chan = 1:length(list_channels)
            loc = find(strcmp(list_channels(chan), ERPdual{1,i}.chans));
            if ~isempty(loc);        allERP{1,i}.ERPdual(chan,:,subject) = ERPdual{1,i}.ERP(loc,:);
            elseif isempty(loc);     allERP{1,i}.ERPdual(chan,:,subject) = NaN;
            end
            clear loc
        end

        % Collect GFPt of all subjects [subject x time]
        allGFP{i}.ERPdual{i}(subject,:) = ERPdual{i}.GFPt;
        
        clear ERPdual{i}
    end
end

%% Plot mean per channel and condition over subjects + mean and std of GFPt
time = -0.1 : 1/1024 : 0.7;

for i = 1:length(EEG.conditions)
    ERPdualGFP1 = mean(allGFP{i}.ERPdual{i},1) + 2*std(allGFP{i}.ERPdual{i},[],1);
    ERPdualGFP2 = mean(allGFP{i}.ERPdual{i},1) - 2*std(allGFP{i}.ERPdual{i},[],1);
    
    figure; 
    subplot(1,2,1); plot(time,mean(allERP{i}.ERPdual,3,'omitnan'),'b');
    subplot(1,2,2); hold on; h1 = fill([time,fliplr(time)], [ERPdualGFP1,fliplr(ERPdualGFP2)],'b','LineStyle','none');
    set(h1,'FaceAlpha',0.4); plot(time,mean(allGFP{i}.ERPdual{i},1),'b','LineWidth',1.5); 
    
    subplot(1,2,1); set(gca,'FontSize',11); box on;
    ylabel('Potential (\muV)','FontSize',14); title('ERP for ',taskname{i},'FontSize',14); ylim([-3 3])
    xticks([-0.1:0.1:0.7]); line([0 0],[-3 3],'Color','k'); 
    subplot(1,2,2); set(gca,'FontSize',11); box on;
    ylabel('GFPt (\muV)','FontSize',14); xlabel('Time (s)','FontSize',14); ylim([-2 22])
    xticks([-0.1:0.1:0.7]); yticks([0:10:20]); line([0 0],[-2 22],'Color','k'); 
    set(gcf,'Position',[400 420 1325 420]);

    clear ERPdualGFP1 ERPdualGFP2 h1

    saveas(gcf,fullfile(results_path, sprintf('ALLERP_%s.png',taskname{1,i})));
end

%% Topoplot of main ERP components per condition
load('chanlocs.mat')

for i = 1:length(EEG.conditions)
% Find the main ERP components in the average GFPt 
% Average GFPt and ERP over subjects
avgGFP{i}.ERPdual{i} = mean(allGFP{i}.ERPdual{i},1);
avgERP{i}.ERPdual{i} = mean(allERP{i}.ERPdual,3, 'omitnan');

% Get main component N0
[ERPdualNO{i}.loc, ERPdualNO{i}.time, ERPdualNO{i}.amp] = get_mainVEPcomponents (time, avgGFP{i}.ERPdual{i}, [-0.02 0.03]);

% Get main component P3: 0.25-0.5s (P300 peak)
[ERPdualP3{i}.loc, ERPdualP3{i}.time, ERPdualP3{i}.amp] = get_mainVEPcomponents (time, avgGFP{i}.ERPdual{i}, [0.25 0.4]);

figure; 
subplot(1,2,1); title('N0'); hold on;
topoplot(avgERP{i}.ERPdual{i}(:,ERPdualNO{i}.loc), chanlocs,'electrodes','off');   colorbar; caxlim(1,:) = caxis; set(gca,'FontSize',9);
ax(1) = gca; axis on; hold on;
axis off; c1 = colorbar; caxlim(1,:) = caxis; 
set(ax,'clim',[-max(caxlim(:,2)) max(caxlim(:,2))]); 
clear caxlim

subplot(1,2,2); title('P300'); hold on;
topoplot(avgERP{i}.ERPdual{i}(:,ERPdualP3{i}.loc), chanlocs,'electrodes','off');   colorbar; caxlim(1,:) = caxis; set(gca,'FontSize',9);
ax(1) = gca; axis on; hold on;
axis off; c2 = colorbar; caxlim(1,:) = caxis; 
set(ax,'clim',[-max(caxlim(:,2)) max(caxlim(:,2))]); 

ylabel(c1,'Potential (\muV)','FontSize',9); 
ylabel(c2,'Potential (\muV)','FontSize',9);

set(gcf,'Position',[400 420 1325 420]);

saveas(gcf,fullfile(results_path, sprintf('ALLERPtopoplot_%s.png',taskname{1,i})));

end

disp('These are the results for the average of all subjects.');
