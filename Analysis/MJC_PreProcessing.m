clear all;
close all;

%% Initialize FieldTrip and EEGLAB
<<<<<<< Updated upstream
% laptop='laptopCatarina';
% laptop='laptopMariana';
laptop='laptopJoao';
=======
laptop='laptopCatarina';
% laptop='laptopMariana';
% laptop='laptopJoao';
>>>>>>> Stashed changes
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);
        
eeglab;
ft_defaults;
[~, ftpath] = ft_version;

sub='28';
rec_nirs='02';
rec_eeg='04';

file_nirs = getFileNames(mainpath_out, sub, rec_nirs);
file_eeg = getFileNames(mainpath_out, sub, rec_eeg);

%% Select data folder 

sub_path = fullfile(mainpath_in,'source',['sub-',sub]);
eeg_path = fullfile(sub_path,'eeg');
nirs_path = fullfile(sub_path,'nirs');
sub_vhdr = fullfile(['sub-',sub,'_rec-',rec_eeg,'_eeg.vhdr']);
nirspre_path = fullfile(mainpath_out, ['sub-',sub], 'nirs');

% Before changing directory to the subpath, add current directory to access
% the function files
addpath(pwd)
cd(sub_path);

oxyfile = fullfile(nirs_path,['sub-',sub,'_rec-',rec_nirs,'_nirs.oxy4']);

%% Upload files
correct = input('Load data (Y/N)? If not, its assumed the .set files have been generated \n', 's');
if strcmpi(correct, 'y')
    done = 1;
    data_loaded = 0;
elseif strcmpi(correct, 'n')
    data_loaded = 1;
end

%% Read data 
if data_loaded == 0
    % FIELDTRIP - load the eeg&nirs data
    cd(nirs_path);
    cfg = [];
    cfg.dataset = oxyfile;
    nirs_raw = ft_preprocessing(cfg);
    nirs_events = ft_read_event(cfg.dataset);
       
    cd(nirspre_path);
    save(['sub-',sub,'_rec-',rec_nirs,'_nirs.mat'], 'nirs_raw');
    save(['sub-',sub,'_rec-',rec_nirs,'_nirs_events.mat'], 'nirs_events');
        
    % EEGLAB load eeg only data
    [EEG,~]         = pop_loadbv(fullfile(sub_path,'eeg'), sub_vhdr);
    [ALLEEG,EEG,~]  = pop_newset(ALLEEG, EEG, 1,'setname','eeg_raw','gui','off','savenew',fullfile(eeg_path,['sub-',sub,'_rec-',rec_eeg,'_eeg']));

else % If data has been loaded and the datasets created, load the structs
    cd(nirs_path);
    load(fullfile(nirspre_path,['sub-',sub,'_rec-',rec_nirs,'_nirs.mat'])); % Avoids call to ft_preprocessing
    load(fullfile(nirspre_path,['sub-',sub,'_rec-',rec_nirs,'_nirs_events.mat'])); % Avoids call to ft_readevents
    [EEG]  = pop_loadset(['sub-',sub,'_rec-',rec_eeg,'_eeg.set'],fullfile(sub_path,'eeg'));
end

%% Read stimuli results
results = load(fullfile(sub_path, 'stim', ['results_sub-',sub,'_rec-',rec_eeg]));
[marker_table,eeg_length,nirs_length]=checkMarkers(EEG, nirs_raw, nirs_events);

%% EEG: Eliminate error moments
% Eliminate moments where the experience wasn't running due to an error.

[EEG] = MJC_rejectErrorMoments(EEG, sub);

%% EEG: Load MNI coordinates
% Load channel coordinates/positions of the standard MNI model of eeglab: 
% MNI dipfit channel positions

[EEG] = pop_chanedit(EEG, 'lookup', join([eeglab_path,...
        '\\plugins\\dipfit\\standard_BESA\\standard-10-5-cap385.elp']),...
        'lookup', join([eeglab_path,...
        '\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc']));

%% EEG: Filter - 50 Hz noise and harmonics

% Filter the signal to obtain the desired frequencies and to eliminate the
% 50 Hz noise
if ~isfile(file_eeg.filtered) 
    % Determine the power spectrum of the raw data
    eeg_raw = EEG.data;
    % [P_raw, f_raw] = periodogram(eeg_raw', [], [] , EEG.srate);
    
    % Filter the data
    eeg_filtered = filterNoise(double(eeg_raw), EEG, 4);
    EEG.data = eeg_filtered;
    [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1, 'setname', 'filtData',...
        'gui', 'off');
    save(file_eeg.filtered, 'EEG');
else
    load(file_eeg.filtered, 'EEG');
    eeg_filtered = EEG.data;
    [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1, 'setname', 'filtData',...
        'gui', 'off');
end

% Determine the power spectrum of the filtered data
% [P_filt, f_filt] = periodogram(eeg_filtered', [], [] , EEG.srate);

% Plot the power spectrums
% figure;
% subplot(1, 3, 1); plot(f_raw, P_raw); 
% xlim([0 200]); ylim([0 7e5]); title('EEG raw data');
% subplot(1, 3, 2); plot(f_filt, P_filt); 
% xlim([0 200]); ylim([0 7e5]); title('EEG filtered data - same scale');
% subplot(1, 3, 3); plot(f_filt, P_filt); 
% xlim([0 50]); title('EEG filtered data - different scale');

%% EEG: Remove bad channels 
% Visually inspect the signals and choose if a signals is too bad that it
% needs to be removed.
% First see the power spectrum and then check if the signal is actually bad
% on the plot.

if ~isfile(file_eeg.removedBadChannels) 
    figure; 
    pop_spectopo(EEG, 1, [0 EEG.pnts], 'EEG', 'percent', 50, 'freqrange',...
        [2 75], 'electrodes', 'off');
    pop_eegplot(EEG);
    RC = input('Remove channel [nr/no]: ','s');
    while ~strcmp(RC, 'no')
        [EEG] = pop_select(EEG, 'nochannel', eval(RC));
        figure;
        pop_spectopo(EEG, 1, [0 EEG.pnts], 'EEG', 'percent', 50,...
            'freqrange', [2 75], 'electrodes', 'off');
        pop_eegplot(EEG);
        RC = input('Remove channel [nr/no]: ','s');
    end
    [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1, 'setname',...
        'removedBadChannels', 'gui', 'off');
    save(file_eeg.removedBadChannels, 'EEG');
else
    load(file_eeg.removedBadChannels, 'EEG');
    [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1, 'setname',...
        'removedBadChannels', 'gui', 'off');
end

%% EEG: Removal of eye blinks - preICA
% Identify the different independent components in the signal

if ~isfile(file_eeg.preICA)  
    [EEG] = pop_runica(EEG,'icatype', 'runica', 'extended', 1,...
        'interrupt', 'on');
    [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1, 'setname', 'preICA',...
        'gui', 'off');
    save(file_eeg.preICA, 'EEG');
else
    load(file_eeg.preICA, 'EEG');
    [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1, 'setname', 'preICA',...
        'gui', 'off');
end

%% EEG: Removal of eye blinks - pstICA
% Visual analysis to remove the component corresponding to eye blinks

if ~isfile(file_eeg.pstICA)
    [EEG] = run_postICA(EEG);
    [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1, 'setname', 'pstICA',...
        'gui', 'off');
    save(file_eeg.pstICA, 'EEG');
else                          
    load(file_eeg.pstICA, 'EEG');
    [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1, 'setname', 'pstICA',...
        'gui', 'off');
end

%% EEG: Set reference
% Re-reference the system to linked mastoids

if ~isfile(file_eeg.preprocessed)
    locs = {EEG.chanlocs.labels};
    M1_loc = find(contains(locs, 'M1'));
    M2_loc = find(contains(locs, 'M2'));
    [EEG] = pop_reref(EEG, [M1_loc, M2_loc]);
    [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1, 'setname',...
        'preprocessed', 'gui', 'off');
    save(file_eeg.preprocessed, 'EEG');
else                          
    load(file_eeg.preprocessed, 'EEG');
    [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1, 'setname',...
        'preprocessed', 'gui', 'off');
end

%% EEG: Extract task data
[EEG_divided, file_eeg] = extractTaskData_EEG(EEG,marker_table, results, file_eeg, mainpath_out);
save(file_eeg.EEG_divided ,'EEG_divided');
[ALLEEG,EEG,~]  = pop_newset(ALLEEG, EEG_divided.EEG_task, 1,'setname','taskData','gui','off');

%% NIRS: Show layout of optode template
cfg = [];
cfg.layout = fullfile(mainpath_out,['sub-',sub],'3d','layout.mat');
ft_layoutplot(cfg);

%% NIRS: Select channels and read events
% Selects channels that are also present in the original template
load(cfg.layout);
names_channels = layout.label(1:(length(layout.label)-2)); % remove 'COMNT' and 'SCALE' from channel labels
TF = startsWith(nirs_raw.label,names_channels);
channels_opt = nirs_raw.label(TF);

% Add all the short channels that are connected to receivers that are
% part of the template
SC = nirs_raw.label(contains(nirs_raw.label, {'a ', 'b ', 'c ', 'd '})); % all short channels
Rx_channels = cellfun(@(x)(strsplit(x, 'T')), names_channels, 'UniformOutput', false); % get the receivers that are part of the template
Rx_channels = cellfun(@(x)(x{1}), Rx_channels, 'UniformOutput', false);
channels_short = SC(startsWith(SC, Rx_channels));

% Select data of correct channels
cd(nirspre_path);
cfg = [];
cfg.channel=[channels_opt; channels_short];
nirs_chan = ft_selectdata(cfg, nirs_raw);
save(['sub-',sub,'_rec-',rec_nirs,'_nirs_chan.mat'], 'nirs_chan');

% Find the event triggers and event types of the dataset (not necessary)
cd(nirs_path);
cfg.dataset = oxyfile;
cfg.trialdef = [];
cfg.trialdef.eventtype = '?';
cfg.trialdef.chanindx = -1;
ft_definetrial(cfg);

% Read events
cd(nirspre_path);
cfg = [];
cfg.dataset = oxyfile;
nirs_events = ft_read_event(cfg.dataset, 'chanindx', -1, 'type', 'event'); % filter out ADC channels (chanindx=-1) and other kind of events ('type'=event)
if strcmp(sub,'02')
    nirs_events(1302:1305)=[];
end
save(['sub-',sub,'_rec-',rec_nirs,'_nirs_events.mat'], 'nirs_events')

% % Select only 0-3000 seconds(to get rid of NaNs)
% cfg = [];
% cfg.latency = [0 3000];
% nirs_chan = ft_selectdata(cfg, nirs_chan);
% save('nirs_raw.mat', 'nirs_raw');

%% NIRS: Downsample the data (save memory, make processing faster and low pass filtering)
% The hemodynamic response takes about 5-10s (0.2-0.1Hz) to reach its peak.
% A 250 Hz measurement is much faster than needed so we need to downsample to 10Hz.  
% If resampling factor is larger than 10 -> resample multiple times
% New frequency must be higher than frequency of trigger

cfg = [];
cfg.resamplefs = 10;
nirs_down = ft_resampledata(cfg, nirs_chan);
save(['sub-',sub,'_rec-',rec_nirs,'_nirs_down.mat'], 'nirs_down');

% % Plot downsampled data
% cfg                = [];
% cfg.preproc.demean = 'yes';
% cfg.viewmode       = 'vertical';
% cfg.continuous     = 'no';
% cfg.ylim           = [ -0.003   0.003 ];
% cfg.channel        = 'Rx*'; % only show channels starting with Rx
% ft_databrowser(cfg, nirs_down);

% % Frequency analysis
% cfg = [];
% cfg.output = 'pow';
% cfg.method = 'mtmfft';
% cfg.taper = 'hanning';
% spectr = ft_freqanalysis(cfg, nirs_down);
% figure;
% hold on;
% plot(spectr.freq, (spectr.powspctrm));
% xlabel('Frequency (Hz)')
% ylabel('Power')
% 
% P=[];F=[];
% for i=1:length(nirs_down.trial)
%     [p, f] = pwelch(nirs_down.trial{i},[],[],[],nirs_down.fsample);
%     P=[P p];
%     F=[F f];
% end
% figure; plot(F(:,1), P); xlabel('Frequency (Hz)'); ylabel('Power');

%% NIRS: Remove bad channels - Inspect the raw data visually
% Show channels with low SCI
addpath(fullfile(ftpath, 'external', 'artinis')); % add artinis functions

cfg = [];
cfg.keepchannel = 'nan';
nirs_sci = ft_nirs_scalpcouplingindex(cfg, nirs_down);

% Show names of bad channels
idx = find(~ismember(nirs_down.label, nirs_sci.label));
bad_nirschannels = nirs_down.label(idx);
disp('The following channels are removed from the data set:');
disp(bad_nirschannels);
cd(nirspre_path);
save(['sub-',sub,'_rec-',rec_nirs,'_nirs_sci.mat'], 'nirs_sci');

databrowser_nirs(nirs_down, 'bad_chan', bad_nirschannels);

% Reject bad channels and trials visually 
cfg = [];
cfg.keepchannel = 'nan';
cfg.method = 'summary';
nirs_reject = ft_rejectvisual(cfg, nirs_sci);

% Make sure that if one channel is rejected, the corresponding channel with
% the other wavelength is rejected too
nirschan_names = cellfun(@(x)(strsplit(x)), nirs_reject.label,'UniformOutput', false);
nirschan_names = cellfun(@(x)(x{1}), nirschan_names, 'UniformOutput', false);
idx = 1; keepnirschan = [];
while idx<length(nirschan_names)-1
  if  strcmp(nirschan_names{idx}, nirschan_names{idx+1})
    keepnirschan = [keepnirschan idx idx+1];
    idx = idx+2;
  else
    idx = idx+1;
  end
end
cfg = [];
cfg.channel = keepnirschan;
nirs_reject = ft_selectdata(cfg, nirs_reject);
save(['sub-',sub,'_rec-',rec_nirs,'_nirs_reject.mat'], 'nirs_reject');

% Visualize
databrowser_nirs(nirs_reject)

%% NIRS: Short channel regression
% Substact the closely located short channel data from the long channels 
% and thus canceling out some noise.
cfg = [];
cfg.method = 'OLS'; %ordinary least square
cfg.verbose = true;
cfg.nearest = false;
<<<<<<< Updated upstream
nirs_reg = shortchannel_regression(cfg, nirs_down);
% save(['sub-',sub,'_rec-',rec_nirs,'_nirs_reg.mat'],'nirs_reg'); 
=======
nirs_reg = shortchannel_regression(cfg, nirs_reject);
save(['sub-',sub,'_rec-',rec_nirs,'_nirs_reg.mat'],'nirs_reg'); 
>>>>>>> Stashed changes

%% NIRS: Detrend and low-pass filtering
% Detrend + low-pass filter (low-pass filter data below the frequency of 
% the heart beat (1 Hz))
cfg = [];
cfg.detrend = 'yes';
cfg.lpfilter = 'yes';
cfg.lpfreq = 0.5; % look at frequency plot
nirs_lpf = ft_preprocessing(cfg, nirs_reg);
cd(nirspre_path);
% save(['sub-',sub,'_rec-',rec_nirs,'_nirs_lpf.mat'],'nirs_lpf'); 

% Visualize
databrowser_nirs(nirs_lpf);

% % Plot Low-pass filtered hemoglobin concentrations 
% idx = find(nirs_lpf.trialinfo==2, 1, 'first'); % check trials
% cfg          = [];
% cfg.channel  = 'Rx*';
% cfg.trials   = idx;
% cfg.baseline = 'yes';
% ft_singleplotER(cfg, nirs_lpf)

%% NIRS: Transform optical densities to concentration changes
cfg = [];
cfg.target = {'O2Hb', 'HHb'};
cfg.channel = 'nirs'; % e.g. one channel incl. wildcards, you can also use ?all? to select all nirs channels
nirs_preprocessed = ft_nirs_transform_ODs(cfg, nirs_lpf);
nirs_preprocessed.sampleinfo = nirs_lpf.sampleinfo;
nirs_preprocessed.hdr = nirs_lpf.hdr;
%nirs_preprocessed.trialinfo = nirs_lpf.trialinfo;
%nirs_trans.time = nirs_lpf.time;

cd(nirspre_path);
save(['sub-',sub,'_rec-',rec_nirs,'_nirs_preprocessed.mat'],'nirs_preprocessed'); 

% Visualize
databrowser_nirs(nirs_preprocessed)

% % Plot hemoglobin concentration over time averaged over all channels for the epoch around the first deviant
% idx = find(nirs_trans.trialinfo==2, 1, 'first'); % check trials
% cfg          = [];
% cfg.channel  = 'Rx*';
% cfg.trials   = idx;
% cfg.baseline = 'yes';
% ft_singleplotER(cfg, nirs_trans)

%% NIRS: Extract task data (epoch)
% Because of the resampling of the data we can't use the current events for 
% epoching so first need to convert them
clear pre post offset trl sel smp

cd(nirspre_path);
% Extract the data
<<<<<<< Updated upstream
[nirs_epoch] = extractTaskData_NIRS(nirs_raw, nirs_lpf, nirs_events, marker_table, sub, rec_nirs);
% save(['sub-',sub,'_rec-',rec_nirs,'_nirs_epoch.mat'], 'nirs_epoch');
=======
[nirs_epoch] = extractTaskData_NIRS(nirs_raw, nirs_preprocessed, nirs_events, marker_table, sub, rec_nirs);
save(['sub-',sub,'_rec-',rec_nirs,'_nirs_epoch.mat'], 'nirs_epoch');
>>>>>>> Stashed changes

% % Plot epoched optical density data around the first deviant stimulus
% idx = find(nirs_epoch.trialinfo==2, 1, 'first'); % check trials
% cfg          = [];
% cfg.channel  = 'Rx*';
% cfg.trials   = idx;
% cfg.baseline = 'yes';
% ft_singleplotER(cfg, nirs_epoch)

% Visualize
% databrowser_nirs(nirs_epoch);

<<<<<<< Updated upstream
%% NIRS: Remove bad channels - Inspect the raw data visually
% Show channels with low SCI
addpath(fullfile(ftpath, 'external', 'artinis')); % add artinis functions

cfg = [];
cfg.keepchannel = 'nan';
nirs_sci = ft_nirs_scalpcouplingindex(cfg, nirs_epoch);

% Show names of bad channels
idx = find(~ismember(nirs_epoch.label, nirs_sci.label));
bad_nirschannels = nirs_down.label(idx);
disp('The following channels are removed from the data set:');
disp(bad_nirschannels);
cd(nirspre_path);
save(['sub-',sub,'_rec-',rec_nirs,'_nirs_sci.mat'], 'nirs_sci');

databrowser_nirs(nirs_epoch, 'bad_chan', bad_nirschannels);

% Reject bad channels and trials visually 
cfg = [];
cfg.keepchannel = 'nan';
cfg.method = 'summary';
nirs_reject = ft_rejectvisual(cfg, nirs_sci);

% Make sure that if one channel is rejected, the corresponding channel with
% the other wavelength is rejected too
nirschan_names = cellfun(@(x)(strsplit(x)), nirs_reject.label,'UniformOutput', false);
nirschan_names = cellfun(@(x)(x{1}), nirschan_names, 'UniformOutput', false);
idx = 1; keepnirschan = [];
while idx<length(nirschan_names)-1
  if  strcmp(nirschan_names{idx}, nirschan_names{idx+1})
    keepnirschan = [keepnirschan idx idx+1];
    idx = idx+2;
  else
    idx = idx+1;
  end
end
cfg = [];
cfg.channel = keepnirschan;
nirs_reject = ft_selectdata(cfg, nirs_reject);
% save(['sub-',sub,'_rec-',rec_nirs,'_nirs_reject.mat'], 'nirs_reject');

% Visualize
databrowser_nirs(nirs_reject)

% % Frequency analysis
% cfg = [];
% cfg.output = 'pow';
% cfg.method = 'mtmfft';
% cfg.taper = 'hanning';
% spectr = ft_freqanalysis(cfg, nirs_epoch);
% figure;
% hold on;
% plot(spectr.freq, (spectr.powspctrm));
% xlabel('Frequency (Hz)')
% ylabel('Power')
% 
% P=[];F=[];
% for i=1:length(nirs_epoch.trial)
%     [p, f] = pwelch(nirs_epoch.trial{i},[],[],[],nirs_epoch.fsample);
%     P=[P p];
%     F=[F f];
% end
% figure; plot(F(:,1), P); xlabel('Frequency (Hz)'); ylabel('Power');

=======
>>>>>>> Stashed changes
% %% NIRS: Short channel regression: Concatenate trials and time data 
% % This step is necessary for the short channel regression to work
% % properly. If not using short channel regression (when
% % plotting only short channel data) then skip this step.
% nirs_all = nirs_reject;
% nirs_all.trial = {cat(2,nirs_reject.trial{:})};
% nirs_all.time = {cat(2,nirs_reject.time{:})};

% %% NIRS: Short channel regression
% % Substact the closely located short channel data from the long channels 
% % and thus canceling out some noise.
% cfg = [];
% cfg.method = 'OLS'; %ordinary least square
% cfg.verbose = true;
% cfg.nearest = false;
% nirs_reg = shortchannel_regression(cfg, nirs_all);
% save(['sub-',sub,'_rec-',rec_nirs,'_nirs_reg.mat'],'nirs_reg'); 

% %% NIRS: Short channel regression: Redefine the trial and time data
% % For proper preprocessing, the trial and time data have to be
% % redefined as previously estimated 
% cfg = [];
% cfg.length = length(nirs_reg.time{1})/length(nirs_reg.trialinfo)/nirs_reg.fsample;
% nirs_new = ft_redefinetrial(cfg, nirs_reg);
% nirs_new.trialinfo = nirs_reject.trialinfo;
% nirs_new.sampleinfo = nirs_reg.sampleinfo;
% nirs_new.time = nirs_reject.time;
% 
% if strcmp(sub,'02')
%     nirs_new.trial{1,80}=zeros(44,296);
% end
% 
% if strcmp(sub,'76')
%     nirs_new.trial{1,80}=zeros(42,323);
% end
% 
% % Get times in trials right
% for i = 1:length(nirs_new.trial)-1
%     if    length(nirs_new.trial{1,i})~=length(nirs_reject.trial{1,i})    
%         
%     for j = 1:abs(length(nirs_new.trial{1,i})-length(nirs_reject.trial{1,i}))
%         while length(nirs_new.trial{1,i})~=length(nirs_reject.trial{1,i})
%             if length(nirs_new.trial{1,i})-length(nirs_reject.trial{1,i}) > 0
%                 s = nirs_new.trial{1,i+1}(:,i);
%                 nirs_new.trial{1,i+1}=[nirs_new.trial{1,i+1} s];
%                 nirs_new.trial{1,i+1}(:,i) = nirs_new.trial{1,i}(:,abs(length(nirs_new.trial{1,i})-length(nirs_reject.trial{1,i}))+1-j);
%                 nirs_new.trial{1,i}(:,abs(length(nirs_new.trial{1,i})-length(nirs_reject.trial{1,i}))+1-j) = [];
%             else
%                 if length(nirs_new.trial{1,i+1}) >=length(nirs_new.trial)
%                            
%                 s = nirs_new.trial{1,i+1}(:,i);
%                 nirs_new.trial{1,i}=[nirs_new.trial{1,i} s];
%                 nirs_new.trial{1,i+1}(:,1) = [];
%                 else
%                 s = nirs_new.trial{1,i+2}(:,i);
%                 nirs_new.trial{1,i+1}=[nirs_new.trial{1,i+1} s];
%                 nirs_new.trial{1,i+2}(:,1) = [];
%                 
%                 s = nirs_new.trial{1,i+1}(:,i);
%                 nirs_new.trial{1,i}=[nirs_new.trial{1,i} s];
%                 nirs_new.trial{1,i+1}(:,1) = [];
%                 end
%             end
%         end
%     end
%     end
% end
% 
% if strcmp(sub,'02')
%     for i = 1:297
%     nirs_new.trial{1,80}(:,1)=[];
%     i=i+1;
%     end
%     d1=abs(length(nirs_reg.trial{1,1})-293);
%     d2=abs(length(nirs_reject.trial{1,80})-size(nirs_new.trial{1,80},2));
%     for i = size(nirs_new.trial{1,80},2):length(nirs_reject.trial{1,80})
%         for j = d1:d1+d2
%             nirs_new.trial{1,80}(:,i)=nirs_reg.trial{1,1}(:,j);
%         end
%     end
% elseif strcmp(sub,'76')
%      for i = 1:324
%     nirs_new.trial{1,80}(:,1)=[];
%     i=i+1;
%     end
%     d1=abs(length(nirs_reg.trial{1,1})-324);
%     d2=abs(length(nirs_reject.trial{1,80})-size(nirs_new.trial{1,80},2));
%     for i = size(nirs_new.trial{1,80},2):length(nirs_reject.trial{1,80})
%         for j = d1:d1+d2
%             nirs_new.trial{1,80}(:,i)=nirs_reg.trial{1,1}(:,j);
%         end
%     end
% else
%     d1=abs(length(nirs_reg.trial{1,1})-length(nirs_new.trial{1,80}));
%     d2=abs(length(nirs_reject.trial{1,80})-length(nirs_new.trial{1,80}));
%     for i = length(nirs_new.trial{1,80}):length(nirs_reject.trial{1,80})
%         for j = d1:d1+d2
%             nirs_new.trial{1,80}(:,i)=nirs_reg.trial{1,1}(:,j);
%         end
%     end
% end
% 
% % Visualize
% databrowser_nirs(nirs_new)

% %% NIRS: Detrend and low-pass filtering
% % Detrend + low-pass filter (low-pass filter data below the frequency of 
% % the heart beat (1 Hz))
% cfg = [];
% cfg.detrend = 'yes';
% cfg.lpfilter = 'yes';
% cfg.lpfreq = 0.5; % look at frequency plot
% nirs_lpf = ft_preprocessing(cfg, nirs_new);
% cd(nirspre_path);
% save(['sub-',sub,'_rec-',rec_nirs,'_nirs_lpf.mat'],'nirs_lpf'); 
% 
% % Visualize
% databrowser_nirs(nirs_lpf);
% 
% % % Plot Low-pass filtered hemoglobin concentrations 
% % idx = find(nirs_lpf.trialinfo==2, 1, 'first'); % check trials
% % cfg          = [];
% % cfg.channel  = 'Rx*';
% % cfg.trials   = idx;
% % cfg.baseline = 'yes';
% % ft_singleplotER(cfg, nirs_lpf)
<<<<<<< Updated upstream

%% NIRS: Transform optical densities to concentration changes
cfg = [];
cfg.target = {'O2Hb', 'HHb'};
cfg.channel = 'nirs'; % e.g. one channel incl. wildcards, you can also use ?all? to select all nirs channels
nirs_preprocessed = ft_nirs_transform_ODs(cfg, nirs_reject);
nirs_preprocessed.sampleinfo = nirs_reject.sampleinfo;
nirs_preprocessed.hdr = nirs_reject.hdr;
nirs_preprocessed.trialinfo = nirs_reject.trialinfo;
%nirs_preprocessed.time = nirs_reject.time;

cd(nirspre_path);
% save(['sub-',sub,'_rec-',rec_nirs,'_nirs_preprocessed.mat'],'nirs_preprocessed'); 

% Visualize
% databrowser_nirs(nirs_preprocessed)

% % Plot hemoglobin concentration over time averaged over all channels for the epoch around the first deviant
% idx = find(nirs_preprocessed.trialinfo==2, 1, 'first'); % check trials
% cfg          = [];
% cfg.channel  = 'Rx*';
% cfg.trials   = idx;
% cfg.baseline = 'yes';
% ft_singleplotER(cfg, nirs_preprocessed)

=======
>>>>>>> Stashed changes
