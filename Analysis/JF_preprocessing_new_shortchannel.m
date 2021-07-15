clear all;
close all;

%% Initialize FieldTrip and EEGLAB
% laptop='laptopCatarina';
% laptop='laptopMariana';
laptop='laptopJoao';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);

eeglab;
ft_defaults;
[~, ftpath] = ft_version;
% select ID number and cap
subjects=[{'02','28','64','76'}];
for iSub = 1:size(subjects,2)
    sub = char(subjects(iSub));
    switch sub
        case '28'
            rec='02';
            rec_nirs=rec;
            rec_eeg=rec;
        case '64'
            rec='01';
            rec_nirs=rec;
            rec_eeg=rec;
        case '02'
            rec='02';
            rec_nirs=rec;
            rec_eeg=rec;
        case '76'
            rec='01';
            rec_nirs='01';
            rec_eeg=rec;
    end

file_nirs = getFileNames(mainpath_out, sub, rec_nirs);
file_eeg = getFileNames(mainpath_out, sub, rec_eeg);

%% Select data folder 
load('C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\Experiment\Data\Pilots\pre-processed\sub-03\stim\sub-03_rec-02_marker_table.mat')

sub_path = fullfile(mainpath_in,'source',['sub-',sub]);
eeg_path = fullfile(sub_path,'eeg');
nirs_path = fullfile(sub_path,'nirs');
sub_vhdr = fullfile(['sub-',sub,'_rec-',rec_eeg,'_eeg.vhdr']);
nirspre_path = fullfile(mainpath_out, ['sub-',sub], 'nirs');

% Before changing directory to the subpath, add current directory to access
% the function files
addpath(pwd)
cd(nirspre_path);

%Load nirs_down
load(['sub-',sub,'_rec-',rec_nirs,'_nirs_reject.mat'])
load(['sub-',sub,'_rec-',rec_nirs,'_nirs.mat'])
load(['sub-',sub,'_rec-',rec_nirs,'_nirs_events.mat'])
%% NIRS: Short channel regression
% Substact the closely located short channel data from the long channels 
% and thus canceling out some noise.
cfg = [];
cfg.method = 'OLS'; %ordinary least square
cfg.verbose = true;
cfg.nearest = false;
nirs_reg = shortchannel_regression(cfg, nirs_reject);
save(['sub-',sub,'_rec-',rec_nirs,'_nirs_reg.mat'],'nirs_reg'); 

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
[nirs_epoch] = extractTaskData_NIRS(nirs_raw, nirs_preprocessed, nirs_events, marker_table, sub, rec_nirs);
save(['sub-',sub,'_rec-',rec_nirs,'_nirs_epoch.mat'], 'nirs_epoch');

% % Plot epoched optical density data around the first deviant stimulus
% idx = find(nirs_epoch.trialinfo==2, 1, 'first'); % check trials
% cfg          = [];
% cfg.channel  = 'Rx*';
% cfg.trials   = idx;
% cfg.baseline = 'yes';
% ft_singleplotER(cfg, nirs_epoch)

% Visualize
% databrowser_nirs(nirs_epoch);



end