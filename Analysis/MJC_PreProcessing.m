clear;
close all;

%% Initialize FieldTrip and EEGLAB
laptop='laptopCatarina';
% laptop='laptopMariana';
% laptop='laptopJoao';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);

eeglab;
ft_defaults;

sub='04';
rec='01';

file = getFileNames(mainpath_out, sub, rec);

%% Select data folder 

sub_path = fullfile(mainpath_in,'source',['sub-',sub]);
eeg_path = fullfile(sub_path,'eeg');
nirs_path = fullfile(sub_path,'nirs');
sub_vhdr = fullfile(['sub-',sub,'_rec-',rec,'_eeg.vhdr']);

% Before changing directory to the subpath, add current directory to access
% the function files
addpath(pwd)
cd(sub_path);

oxyfile = fullfile(nirs_path,['sub-',sub,'_rec-',rec,'_nirs.oxy3']);

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
    if strcmp(sub,'03') && strcmp(rec,'02')
        nirs_events=ft_filter_event(nirs_events,'minsample',72787);
    end
    
    if strcmp(sub,'02') && strcmp(rec,'02')
        save(['sub-',sub,'_rec-',rec,'_nirseeg.mat'], 'nirs_raw');
        save(['sub-',sub,'_rec-',rec,'_nirseeg_events.mat'], 'nirs_events');
    else
        save(['sub-',sub,'_rec-',rec,'_nirs.mat'], 'nirs_raw');
        save(['sub-',sub,'_rec-',rec,'_nirs_events.mat'], 'nirs_events');
    end
    
    % EEGLAB load eeg only data
    [EEG,~]         = pop_loadbv(fullfile(sub_path,'eeg'), sub_vhdr);
    [ALLEEG,EEG,~]  = pop_newset(ALLEEG, EEG, 1,'setname','eeg_raw','gui','off','savenew',fullfile(eeg_path,['sub-',sub,'_rec-',rec,'_eeg']));

else % If data has been loaded and the datasets created, load the structs
    if strcmp(sub,'02') && strcmp(rec,'02')
        load(nirs_path,['sub-',sub,'_rec-',rec,'_nirseeg.mat']); % Avoids call to ft_preprocessing
        load(nirs_path,['sub-',sub,'_rec-',rec,'_nirseeg_events.mat']); % Avoids call to ft_readevents
    else
        load(fullfile(nirs_path,['sub-',sub,'_rec-',rec,'_nirs.mat'])); % Avoids call to ft_preprocessing
        load(fullfile(nirs_path,['sub-',sub,'_rec-',rec,'_nirs_events.mat'])); % Avoids call to ft_readevents
    end
    [EEG]  = pop_loadset(['sub-',sub,'_rec-',rec,'_eeg.set'],fullfile(sub_path,'eeg'));
end

%% Read stimuli results
if strcmp(sub,'02') && strcmp(rec,'02')
    nirs_raw = data_raw;
    nirs_events=eeg_fnirs_events;
end

results = load(fullfile(sub_path, 'stim', ['results_sub-',sub,'_rec-',rec]));
marker_table = checkMarkers(EEG, nirs_raw, nirs_events);

%% EEG: Load MNI coordinates
% Load channel coordinates/positions of the standard MNI model of eeglab: 
% MNI dipfit channel positions

[EEG] = pop_chanedit(EEG, 'lookup', join([eeglab_path,...
        '\\plugins\\dipfit\\standard_BESA\\standard-10-5-cap385.elp']),...
        'lookup', join([eeglab_path,...
        '\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc']));

%% EEG: Filter - 50 Hz noise and harmonics

% Determine the power spectrum of the raw data
eeg_raw = EEG.data;
% [P_raw, f_raw] = periodogram(eeg_raw', [], [] , EEG.srate);

% Filter the signal to obtain the desired frequencies and to eliminate the
% 50 Hz noise
if ~isfile(file.filtered) 
    eeg_filtered = filterNoise(double(eeg_raw), EEG, 4);
    EEG.data = eeg_filtered;
    [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1, 'setname', 'filtData',...
        'gui', 'off');
    save(file.filtered, 'EEG');
else
    load(file.filtered, 'EEG');
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

if ~isfile(file.removedBadChannels) 
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
    save(file.removedBadChannels, 'EEG');
else
    load(file.removedBadChannels, 'EEG');
    [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1, 'setname',...
        'removedBadChannels', 'gui', 'off');
end

%% EEG: Removal of eye blinks - preICA
% Identify the different independent components in the signal

if ~isfile(file.preICA)  
    [EEG] = pop_runica(EEG,'icatype', 'runica', 'extended', 1,...
        'interrupt', 'on');
    [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1, 'setname', 'preICA',...
        'gui', 'off');
    save(file.preICA, 'EEG');
else
    load(file.preICA, 'EEG');
    [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1, 'setname', 'preICA',...
        'gui', 'off');
end

%% EEG: Removal of eye blinks - pstICA
% Visual analysis to remove the component corresponding to eye blinks

if ~isfile(file.pstICA)
    [EEG] = run_postICA(EEG);
    [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1, 'setname', 'pstICA',...
        'gui', 'off');
    save(file.pstICA, 'EEG');
else                          
    load(file.pstICA, 'EEG');
    [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1, 'setname', 'pstICA',...
        'gui', 'off');
end

%% EEG: Set reference
% Re-reference the system to linked mastoids

locs = {EEG.chanlocs.labels};
M1_loc = find(contains(locs, 'M1'));
M2_loc = find(contains(locs, 'M2'));

if ~isfile(file.preprocessed)
    [EEG] = pop_reref(EEG, [M1_loc, M2_loc]);
    [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1, 'setname',...
        'preprocessed', 'gui', 'off');
    save(file.preprocessed, 'EEG');
else                          
    load(file.preprocessed, 'EEG');
    [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1, 'setname',...
        'preprocessed', 'gui', 'off');
end

%% EEG: Extract task data
[EEG_divided, file] = extractTaskData_EEG(EEG,marker_table, results, file, mainpath_out);
save(file.EEG_divided,'EEG_divided');
% [ALLEEG,EEG,~]  = pop_newset(ALLEEG, EEG_task, 1,'setname','taskData','gui','off');

%% NIRS: Show layout of optode template
cfg = [];
cfg.layout = fullfile(mainpath_out,['sub-',sub],'3d','layout.mat');
ft_layoutplot(cfg);

%% NIRS: Select channels
% Selects channels that are also present in the original template
load(cfg.layout);
channels_conc = layout.label(1:(length(layout.label)-2)); % remove 'COMNT' and 'SCALE' from channel labels
for i = 1:length(channels_conc)
    tmp = strsplit(channels_conc{i});
    names_channels{i}=char(tmp(1));
end
TF = startsWith(nirs_raw.label,names_channels);
channels_opt = nirs_raw.label(TF);

% Select data of correct channels
cd(nirs_path);
cfg = [];
cfg.inputfile = ['sub-',sub,'_rec-',rec,'_nirs.mat']; 
cfg.channel = channels_opt;
nirs_chan = ft_selectdata(cfg);
cd(nirs_path);
save('nirs_chan.mat', 'nirs_chan');

%% NIRS: Downsample the data (save memory, make processing faster and low pass filtering)
% The hemodynamic response takes about 5-10s (0.2-0.1Hz) to reach its peak.
% A 250 Hz measurement is much faster than needed so we need to downsample to 10Hz.  
% If resampling factor is larger than 10 -> resample multiple times
% New frequency must be higher than frequency of trigger

if ~exist(fullfile(['sub-',sub,'_rec-',rec,'_nirs.mat'],'nirs_raw'))
    cfg = [];
    cfg.inputfile = 'nirs_chan.mat';
    cfg.resamplefs = 10;
    nirs_down = ft_resampledata(cfg);
    cd(nirs_path);
    save('nirs_down.mat', 'nirs_down');
else load('nirs_down.mat')
end

% % Plot downsampled data
% cfg                = [];
% cfg.preproc.demean = 'yes';
% cfg.viewmode       = 'vertical';
% cfg.continuous     = 'no';
% cfg.ylim           = [ -0.003   0.003 ];
% cfg.channel        = 'Rx*'; % only show channels starting with Rx
% ft_databrowser(cfg, nirs_down);

%% NIRS: Extract task data (epoch)
clear pre post offset trl sel smp
load('nirs_down.mat');
load(['sub-',sub,'_rec-',rec,'_nirs.mat']);

if ~exist('nirs_epoch.mat') 
    % Find the event triggers and event types of the dataset
    cfg.dataset = oxyfile;
    cfg.trialdef = [];
    cfg.trialdef.eventtype = '?';
    cfg.trialdef.chanindx = -1;
    ft_definetrial(cfg);
    cfg = [];
    cfg.dataset = oxyfile;
    event = ft_read_event(cfg.dataset, 'chanindx', -1, 'type', 'event'); % filter out ADC channels (chanindx=-1) and other kind of events ('type'=event)
    
    % Extract the data
    [nirs_epoch] = extractTaskData_NIRS(nirs_raw, nirs_down, event, marker_table);
    cd(nirs_path);
    save('nirs_epoch.mat', 'nirs_epoch');
    
else
    load('nirs_epoch.mat')
end

% % Plot epoched optical density data around the first deviant stimulus
% idx = find(nirs_epoch.trialinfo==2, 1, 'first'); % check trials
% cfg          = [];
% cfg.channel  = 'Rx*';
% cfg.trials   = idx;
% cfg.baseline = 'yes';
% ft_singleplotER(cfg, nirs_epoch)

%% NIRS: Remove bad channels - Inspect the raw data visually 
% Show channels with low SCI
load('nirs_down.mat')
cfg = [];
nirs_sci = ft_nirs_scalpcouplingindex(cfg, nirs_down);

% Show names of bad channels
idx = find(~ismember(nirs_down.label, nirs_sci.label));
bad_nirschannels = nirs_down.label(idx);
disp('The following channels are removed from the data set:');
disp(bad_nirschannels);
cd(nirs_path);
save('nirs_sci.mat', 'nirs_sci');

% Inspect the data visually per channel
load('nirs_chan.mat'); % it's not possible to plot events on data that has been resampled
cfg = [];
cfg.preproc.demean = 'yes'; % substracts the mean value in the plot
cfg.viewmode = 'vertical';
cfg.event = ft_read_event(oxyfile, 'chanindx', -1, 'type', 'event');
cfg.ploteventlabels= 'colorvalue';
cfg.plotlabels= 'yes';
cfg.fontsize = 5;
cfg.continuous = 'yes'; 
cfg.blocksize = 300;
cfg.nirsscale = 10;
cfg.channel = 1:2; % [channels_opt(1), channels_opt(2)];
cfg.linecolor = 'brmm';
cfg.colorgroups = repmat([1 2],1, length(nirs_chan.label)/2);
ft_databrowser(cfg, nirs_chan);

% Inspect the data visually for each trial
load('nirs_epoch.mat'); 
cfg = [];
cfg.preproc.demean = 'yes'; % substracts the mean value in the plot so that the channels are well visualized above each other
cfg.viewmode = 'butterfly';
cfg.continuous = 'no';
cfg.ylim = [ -0.003   0.003 ];
cfg.channel = 'Rx*';
cfg.linecolor = 'br';
cfg.colorgroups = repmat([1 2],1, length(nirs_chan.label)/2);
ft_databrowser(cfg, nirs_epoch);

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

%% NIRS: Transform optical densities to concentration changes
load('nirs_epoch.mat');
cfg = [];
cfg.target = {'O2Hb', 'HHb'};
cfg.channel = 'nirs'; % e.g. one channel incl. wildcards, you can also use ?all? to select all nirs channels
nirs_conc = ft_nirs_transform_ODs(cfg, nirs_epoch);
cd(nirs_path);
save('nirs_conc.mat','nirs_conc'); 

% %Plot hemoglobin concentration over time averaged over all channels for the epoch around the first deviant
% idx = find(nirs_conc.trialinfo==2, 1, 'first'); % check trials
% cfg          = [];
% cfg.channel  = 'Rx*';
% cfg.trials   = idx;
% cfg.baseline = 'yes';
% ft_singleplotER(cfg, nirs_conc)

%% NIRS: Look at the data per condition
% multiplot_condition doesn't exist??

%% NIRS: Low-pass filtering
% Low-pass filter data below the frequency of the heart beat (1 Hz).
cfg = [];
cfg.inputfile = 'nirs_conc.mat'; 
cfg.lpfilter = 'yes';
cfg.lpfreq = 0.2; % look at frequency plot
nirs_lpf = ft_preprocessing(cfg);
cd(nirs_path);
save('nirs_lpf.mat','nirs_lpf'); 

% %Plot Low-pass filtered hemoglobin concentrations 
% idx = find(nirs_lpf.trialinfo==2, 1, 'first'); % check trials
% cfg          = [];
% cfg.channel  = 'Rx*';
% cfg.trials   = idx;
% cfg.baseline = 'yes';
% ft_singleplotER(cfg, nirs_lpf)

% to have a look at data multiplot_condition doesn't exist??

%% NIRS: Timelock analysis and baseline correction
% Timelock analysis
load('nirs_lpf.mat')
for task = 1:8 % 8 conditions 
    cfg = [];
    cfg.trials = find(nirs_lpf.trialinfo(:,1) == task); % average the data for given task
    nirs_TL{task} = ft_timelockanalysis(cfg, nirs_lpf);
end
cd(nirs_path);
save('nirs_TL.mat', 'nirs_TL');

% Apply baseline correction
for task = 1:8
    cfg = [];
    cfg.baseline = [-10 0]; % define the amount of seconds you want to use for the baseline
    nirs_TL_blc{task} = ft_timelockbaseline(cfg, nirs_TL{task});
end
cd(nirs_path);
save('nirs_TL_blc.mat','nirs_TL_blc');

%% NIRS: Separate O2Hb and HHb channels and plot tasks on the layout
close all;
load('nirs_TL_blc.mat')
for task = 1:8
    cfg = [];
    cfg.channel='* [O2Hb]';
    nirs_TL_O2Hb{task} = ft_preprocessing(cfg, nirs_TL_blc{task});
    % rename labels so that they have the same name as HHb channels
    for i = 1:length(nirs_TL_O2Hb{task}.label)
        tmp = strsplit(nirs_TL_O2Hb{task}.label{i});
        nirs_TL_O2Hb{task}.label{i} = tmp{1};
    end
    save('nirs_TL_O2Hb.mat','nirs_TL_O2Hb');
    
    % same for HHb channels
    cfg = [];
    cfg.channel = '* [HHb]';
    nirs_TL_HHb{task} = ft_preprocessing(cfg, nirs_TL_blc{task});
    for i = 1:length(nirs_TL_HHb{task}.label)
        tmp = strsplit(nirs_TL_HHb{task}.label{i});
        nirs_TL_HHb{task}.label{i} = tmp{1};
    end
    cd(nirs_path);
    save('nirs_TL_HHb.mat','nirs_TL_HHb');
end

load('nirs_TL_O2Hb.mat');
load('nirs_TL_HHb.mat');
for i = [5 1]
  cfg = [];
  cfg.showlabels = 'yes';
  cfg.layout = layout;
  cfg.ylim = [-1 1.5];
  cfg.interactive = 'yes'; % this allows to select a subplot and interact with it
  cfg.linecolor  = 'rbmcgykw'; % O2Hb is showed in red/magenta/green/black and HHb in blue/cyan/yellow/white
  cfg.comment = 'auto dual task is red and blue line\n auto single task is magenta and cyan line\n non auto dual task is green and yellow\n non auto single task is black and white';
  taskshort = {'multiplotuncued','','','','multiplotcued'};
  ft_multiplotER(cfg, nirs_TL_O2Hb{i}, nirs_TL_HHb{i}, nirs_TL_O2Hb{i+1}, nirs_TL_HHb{i+1}, nirs_TL_O2Hb{i+2}, nirs_TL_HHb{i+2},nirs_TL_O2Hb{i+3}, nirs_TL_HHb{i+3});
  saveas(gcf, ['sub-',sub,'_rec-',rec,'_',char(taskshort(i)), '_timelock.jpg']);
end

% change white line to orange, edit title and legend of the plot and save
% the image
savefigure = input('Enter any key once the plots are edited and saved: \n', 's');

% Plot for each task separately
taskname = {'Auto Dual Task Cued', 'Auto Single Task Cued', 'Non Auto Dual Task Cued', 'Non Auto Single Task Cued', 'Auto Dual Task Uncued', 'Auto Single Task Uncued', 'Non Auto Dual Task Uncued', 'Non Auto Single Task Uncued'};
taskshort = {'autodualcued', 'autosinglecued', 'nonautodualcued', 'nonautosinglecued','autodualuncued', 'autosingleuncued', 'nonautodualuncued', 'nonautosingleuncued'};
for task = 1:8
    cfg = [];
    cfg.showlabels = 'yes';
    cfg.layout = layout;
    cfg.interactive = 'yes'; % this allows to select a subplot and interact with it
    cfg.linecolor = 'rb'; % O2Hb is in red and HHb in blue
    ft_multiplotER(cfg, nirs_TL_O2Hb{task}, nirs_TL_HHb{task});
    saveas(gcf, ['sub-',sub,'_rec-',rec,'_',char(taskshort(task)), '_timelock.jpg']);
end

%% NIRS: Spatial representation
close all;
nirs_NAvA_autodualcued = nirs_TL_blc{1,1};
nirs_NAvA_autosinglecued = nirs_TL_blc{1,2};
nirs_NAvA_nonautodualcued = nirs_TL_blc{1,3};
nirs_NAvA_nonautosinglecued = nirs_TL_blc{1,4};
nirs_NAvA_autodualuncued = nirs_TL_blc{1,5};
nirs_NAvA_autosingleuncued = nirs_TL_blc{1,6};
nirs_NAvA_nonautodualuncued = nirs_TL_blc{1,7};
nirs_NAvA_nonautosingleuncued = nirs_TL_blc{1,8};

nirs_NAvA = [nirs_NAvA_autodualcued, nirs_NAvA_autosinglecued, nirs_NAvA_nonautodualcued, nirs_NAvA_nonautosinglecued,nirs_NAvA_autodualuncued, nirs_NAvA_autosingleuncued,nirs_NAvA_nonautodualuncued,nirs_NAvA_nonautosingleuncued];

% the plot doesn't look like the layout??
% Plot the spatial representation
load('layout.mat')
figure; ft_plot_layout(layout)
cd(nirs_path);
for task = 1:8
    cfg = [];
    cfg.layout = sort(layout.label); 
    cfg.marker = 'labels';
    cfg.xlim = [5 10]; % choose the time window over which you want to average
    cfg.zlim = [-0.5 0.5]; % choose the window over which you want to scale based on the plots
    cfg.channel  = '* [O2Hb]';
    nirs_NAvA(task).label = sort(nirs_NAvA(task).label);
    
    % Plot [O2Hb] channels
    figure; subplot(1,2,1);
    ft_topoplotER(cfg, nirs_NAvA(task));
    title('Topoplot average [O2Hb]'); 
    % Plot [HHb] channels
    cfg.channel  = '* [HHb]';
    subplot(1,2,2);
    ft_topoplotER(cfg, nirs_NAvA(task));
    title('Topoplot average [HHb]');
    saveas(gcf, ['sub-',sub,'_rec-',rec,'_',char(taskshort(task)), '_spatial.jpg']);
end
