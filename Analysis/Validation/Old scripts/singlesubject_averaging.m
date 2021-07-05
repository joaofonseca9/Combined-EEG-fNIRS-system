%% 1. Define parameters
clear all
root_dir= 'C:\Users\mtabo\Documents\Data_Internship';

sub='sub-90';
filename=fullfile(root_dir, 'source_data', sub, 'nirs', sprintf('%s_task-automaticity_acq-online_rec-01_nirs.oxy4', sub));

load(fullfile(root_dir, 'layout.mat')); % load layout
addpath(fullfile(root_dir, 'pilot_Arne-master')); % add own functions
[~, ftpath]=ft_version;

% save data in processed directory
proc_folder=fullfile(root_dir, 'processedshortchanreg', sub);
if ~exist(proc_folder)
  mkdir(proc_folder)
end
cd(proc_folder)


%% 2. Read data & events
% a) Read data
cfg = [];
cfg.dataset = filename;
[data_raw] = ft_preprocessing(cfg);
save('data_raw.mat', 'data_raw');

% b) Show layout of optode template
cfg = [];
cfg.layout= layout;
ft_layoutplot(cfg);

% c) Select channels that were used in the experiment (too many channels
% are read in data_raw)
% ! This only selects channels that are also present in the original
% template! -> so the layout template?
names_channels = layout.label(1:(length(layout.label)-2)); % remove 'COMNT' and 'SCALE' from channellabels
TF=startsWith(data_raw.label,names_channels);
channels_opt=data_raw.label(TF);
 % also add all short channels that are connected to receivers that are
 % part of the templates
SC=data_raw.label(contains(data_raw.label, {'a ', 'b ', 'c ', 'd '})); % all short channels
Rx_channels=cellfun(@(x)(strsplit(x, 'T')), names_channels, 'UniformOutput', false); % get the receivers that are part of the template
Rx_channels=cellfun(@(x)(x{1}), Rx_channels, 'UniformOutput', false);
channels_short=SC(startsWith(SC, Rx_channels));

% d) Select data of correct channels
cfg=[];
cfg.channel=[channels_opt; channels_short];
%cfg.channel=[channels_opt]; %for plotting short channels only!!
data_chan=ft_selectdata(cfg, data_raw);
save('data_chan.mat', 'data_chan');

% e) read in events
% % Find the event triggers and event types of the dataset (no necessary
% % step)
% cfg.dataset=filename;
% cfg.trialdef = [];
% cfg.trialdef.eventtype = '?';
% cfg.trialdef.chanindx=-1;
% ft_definetrial(cfg);

% read in events
cfg = [];
cfg.dataset = filename;
event=ft_read_event(cfg.dataset, 'chanindx', -1, 'type', 'event'); % filter out ADC channels (chanindx=-1) and other kind of events ('type'=event)
save('event.mat', 'event')

% % b)Select only 0-3000 seconds(to get rid of NaNs)
%    cfg=[];
% cfg.latency=[0 3000];
% data_chan=ft_selectdata(cfg, data_chan);
% save('data_raw.mat', 'data_raw');


%% 3. Downsample the data (to save memory, make processing faster and low pass filtering)
% HDR = 5-10" --> freq 0.2-0.1 Hz
% ! if resampling factor is larger than 10 --> resample multiple times
% ! new frequency must be higher than frequency of trigger
cfg                   = [];
cfg.resamplefs        = 10;
data_down             = ft_resampledata(cfg, data_chan);
save('data_down.mat', 'data_down');

 %% 4. Epoch 
 % because of resampling of the data, we cannot use the current events for epoching, but first need to convert them
 clear pre post offset trl sel smp
 
 % Define events + epochs
 foot_auto = find(strcmp({event.value}, 'LSL foot_auto'));
 foot_nonauto =  find(strcmp({event.value}, 'LSL foot_nonauto'));
 finger_auto= find(strcmp({event.value}, 'LSL finger_auto'));
 finger_nonauto= find(strcmp({event.value}, 'LSL finger_nonauto'));
 % get the sample number in the original dataset
 % note that we transpose them to get columns
 smp.foot_auto = [event(foot_auto).sample]';
 smp.foot_nonauto = [event(foot_nonauto).sample]';
 smp.finger_auto = [event(finger_auto).sample]';
 smp.finger_nonauto = [event(finger_nonauto).sample]';
 
 factor = data_raw.fsample / data_down.fsample;
 
 % get the sample number after downsampling
 % ! maybe convert this structure-array into a matrix?
 smp.foot_auto = round((smp.foot_auto-1)/factor+1);
 smp.foot_nonauto= round((smp.foot_nonauto-1)/factor+1);
 smp.finger_auto = round((smp.finger_auto-1)/factor+1);
 smp.finger_nonauto= round((smp.finger_nonauto-1)/factor+1);
 
 pre    =  round(10*data_down.fsample); % seconds pre-stimulus (100=10 seconden)
 post   =  round(20*data_down.fsample); % seconds post-stimulus
 offset = -pre; % see ft_definetrial
 
 trl.foot_auto = [smp.foot_auto-pre smp.foot_auto+post];
 trl.foot_nonauto = [smp.foot_nonauto-pre smp.foot_nonauto+post];
 trl.finger_auto = [smp.finger_auto-pre smp.finger_auto+post];
 trl.finger_nonauto = [smp.finger_nonauto-pre smp.finger_nonauto+post];
 
 % add the offset
 trl.foot_auto(:,3) = offset;
 trl.foot_nonauto(:,3) = offset;
 trl.finger_auto(:,3) = offset;
 trl.finger_nonauto(:,3) = offset;
 
 % trialinfo
 trl.foot_auto(:,4) = 3;
 trl.foot_nonauto(:,4) = 4;
 trl.finger_auto(:,4) = 1;
 trl.finger_nonauto(:,4) = 2;
 
 % concatenate the four conditions and sort them
 trl = sortrows([trl.foot_auto; trl.foot_nonauto; trl.finger_auto; trl.finger_nonauto]);
 
 % remove trials that stretch beyond the end of the recording
 sel = trl(:,2)<size(data_down.trial{1},2);
 trl = trl(sel,:);
 
 cfg     = [];
 cfg.trl = trl;
 data_epoch = ft_redefinetrial(cfg, data_down);
 save('data_epoch.mat', 'data_epoch');

 % visualize
 databrowser_nirs(data_epoch);
 %ft_databrowser(cfg, data_epoch);


%% 5. remove bad channels
% replaces data of bad channels with nans
% a) Show channels with low SCI
addpath(fullfile(ftpath, 'external', 'artinis')); % add artinis functions
cfg                 = [];
cfg.keepchannel = 'nan';
data_sci            = ft_nirs_scalpcouplingindex(cfg,data_epoch);
% Show names of bad channels
idx=find(~ismember(data_epoch.label,data_sci.label));
bad_channels=data_down.label(idx);
disp('The following channels are removed from the data set:');
disp(bad_channels);
save('data_sci.mat', 'data_sci');

%databrowser_nirs(data_epoch, 'bad_chan', bad_channels);
databrowser_nirs(data_epoch, 'bad_chan', bad_channels);
%%
% b) reject visually bad channels and trials
cfg=[];
cfg.keepchannel = 'nan';
cfg.method='summary';

data_reject=ft_rejectvisual(cfg, data_sci);
% in this case, I rejected: 
% channels: 28, 38, 47,48
% trials: 3, 9, 11, 12, 17
% that is quite conservative

% make sure that if one channel is rejected, the corresponding channel with
% the other wavelength is rejected too
chan_names=cellfun(@(x)(strsplit(x)), data_reject.label,'UniformOutput', false);
chan_names=cellfun(@(x)(x{1}), chan_names, 'UniformOutput', false);
idx=1; keepchan=[];
while idx<length(chan_names)-1
  if  strcmp(chan_names{idx}, chan_names{idx+1})
    keepchan=[keepchan idx idx+1];
    idx=idx+2;
  else
    idx=idx+1;
  end
end
cfg=[];
cfg.channel=keepchan;
data_reject2=ft_selectdata(cfg, data_reject)
save('data_reject2.mat', 'data_reject2');

% visualize
databrowser_nirs(data_reject2)


%% 6. Concatenate trials and time data 
% this step is necessary in order for short channel regression to work
% properly (if not using short channel regression (for example when
% plotting only short channel data) then skip this step
data_all = data_reject2;
data_all.trial = {cat(2,data_reject2.trial{:})};
data_all.time = {cat(2,data_reject2.time{:})};


%% 7. Short channel regression
% this function finds the short channels closely located to long channels
% and substracts this data from the long channels
cfg = [];
cfg.method = 'OLS'
cfg.verbose = true
cfg.nearest = false
outdata = shortchannel_regression(cfg, data_all);
save('outdata.mat', 'outdata');


%% 8. Redefine the trial and time data
% in order for proper preprocessing the trial and time data has to be
% redefined as previously estimated again 
cfg = [];
cfg.length = length(outdata.time{1})/length(outdata.trialinfo)/outdata.fsample;
data_new = ft_redefinetrial(cfg, outdata);
data_new.trialinfo = data_reject2.trialinfo;
data_new.sampleinfo = outdata.sampleinfo;

databrowser_nirs(data_new)
%% 9. data preprocessing 
% detrend + low pass filter
cfg=[];
cfg.detrend='yes';
cfg.lpfilter          = 'yes';
cfg.lpfreq            = 0.5; 
data_preproc=ft_preprocessing(cfg, data_new);
%data_preproc=ft_preprocessing(cfg, data_reject2); %without short channel regression
save('data_preproc.mat', 'data_preproc')

% visualize
databrowser_nirs(data_preproc);


%% 10. Transform optical densities to concentration changes
% a) Transformation
cfg                 = [];
cfg.target          = {'O2Hb', 'HHb'};
cfg.channel         = 'nirs'; % e.g. one channel incl. wildcards, you can also use ?all? to select all nirs channels
data_conc           = ft_nirs_transform_ODs(cfg, data_preproc);
data_conc.sampleinfo = data_reject2.sampleinfo;
data_conc.hdr = data_reject2.hdr;
data_conc.time = data_reject2.time;
%data_conc = removefields(data_conc, 'data_conc.sampleinfo')
save('data_conc.mat','data_conc'); 

databrowser_nirs(data_conc)

%% 11. Visualize the data
h=multiplot_condition(data_conc, layout, [1:4], {'finger auto', 'finger nonauto', 'foot auto','foot nonauto'}, 'baseline', [-10 0], 'trials', false, 'topoplot', 'yes', 'ylim', [-0.2 0.2]);
% savefig(h{1}, 'hand_auto_multi'); savefig(h{3}, 'hand_nonauto_mulitplot');
% savefig(h{5}, 'foot_auto_multi'); savefig(h{7}, 'foot_nonauto_mulitplot');


%% 12. Statistical testing
con_names={'finger auto',  'finger nonauto', 'foot auto', 'foot nonauto'};
for i=1:4 % loop over the 4 conditions
  [stat_O2Hb, stat_HHb] = statistics_withinsubjects(data_conc, 'data_conc', layout, i, con_names{i});
end

%% 13. Timelock analysis + baselinecorrection
% these steps are necessary for the multisubject_averaging script!
% a) timelock
for task=1:4 %There are 4 tasks (handauto, handnonauto, footauto, footnonauto)
    cfg               = [];
    cfg.trials        = find(data_conc.trialinfo(:,1)==task); % Average the data for given task
    data_TL{task} = ft_timelockanalysis(cfg, data_conc);
end
save('data_TL.mat', 'data_TL');

% b) apply a baseline correction
for task=1:4
    cfg                 = [];
    cfg.baseline        = [-10 0]; % define the amount of seconds you want to use for the baseline
    data_TL_blc{task}     = ft_timelockbaseline(cfg, data_TL{task});
end
save('data_TL_blc.mat','data_TL_blc');



%% HELPER FUNCTIONS
function databrowser_nirs(data, varargin)
bad_chan=ft_getopt(varargin, 'bad_chan', {}); % bad channels are indicated in black
  cfg                = [];
  cfg.preproc.demean = 'yes'; % substracts the mean value (only in the plot)
  cfg.viewmode       = 'vertical';
  cfg.channel = {'Rx*'};
  cfg.fontsize=5;
  cfg.blocksize  = 30;
  cfg.nirsscale =150;
  cfg.ylim = [-1 1];
  cfg.linecolor = 'brkk';
  cfg.colorgroups = repmat([1 2],1, length(data.label)/2)+2*ismember(data.label, bad_chan)';
  ft_databrowser(cfg, data);
end