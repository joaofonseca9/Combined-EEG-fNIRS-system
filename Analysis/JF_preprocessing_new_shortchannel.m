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
load(['sub-',sub,'_rec-',rec_nirs,'_nirs_down.mat'])
load(['sub-',sub,'_rec-',rec_nirs,'_nirs.mat'])
load(['sub-',sub,'_rec-',rec_nirs,'_nirs_events.mat'])
%% NIRS: Short channel regression
% Substact the closely located short channel data from the long channels 
% and thus canceling out some noise.
cfg = [];
cfg.method = 'OLS'; %ordinary least square
cfg.verbose = true;
cfg.nearest = false;
nirs_reg = shortchannel_regression(cfg, nirs_down);
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
save(['sub-',sub,'_rec-',rec_nirs,'_nirs_lpf.mat'],'nirs_lpf'); 

% Visualize
databrowser_nirs(nirs_lpf);

% % Plot Low-pass filtered hemoglobin concentrations 
% idx = find(nirs_lpf.trialinfo==2, 1, 'first'); % check trials
% cfg          = [];
% cfg.channel  = 'Rx*';
% cfg.trials   = idx;
% cfg.baseline = 'yes';
% ft_singleplotER(cfg, nirs_lpf)

%% NIRS: Extract task data (epoch)
% Because of the resampling of the data we can't use the current events for 
% epoching so first need to convert them
clear pre post offset trl sel smp

cd(nirspre_path);
% Extract the data
[nirs_epoch] = extractTaskData_NIRS(nirs_raw, nirs_lpf, nirs_events, marker_table, sub, rec_nirs);
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
save(['sub-',sub,'_rec-',rec_nirs,'_nirs_reject.mat'], 'nirs_reject');

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
save(['sub-',sub,'_rec-',rec_nirs,'_nirs_preprocessed.mat'],'nirs_preprocessed'); 

% Visualize
% databrowser_nirs(nirs_preprocessed)

% % Plot hemoglobin concentration over time averaged over all channels for the epoch around the first deviant
% idx = find(nirs_preprocessed.trialinfo==2, 1, 'first'); % check trials
% cfg          = [];
% cfg.channel  = 'Rx*';
% cfg.trials   = idx;
% cfg.baseline = 'yes';
% ft_singleplotER(cfg, nirs_preprocessed)

end