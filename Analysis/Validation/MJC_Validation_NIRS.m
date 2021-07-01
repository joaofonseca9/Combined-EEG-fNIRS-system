clear, close all
%% Settings

addpath (fullfile(pwd,'..')) %add the path with analysis scripts
% add path to correct folders and open eeglab
laptop = 'laptopJoao';
[mainpath_in, ~, eeglab_path] = addFolders(laptop);
mainpath_in_processed=fullfile(mainpath_in,'processed');
mainpath_in_preprocessed=fullfile(mainpath_in,'pre-processed');
mainpath_out = 'C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\Experiment\Combined-EEG-fNIRS-system\Analysis\Validation';

mainpath_in_nirsonly='C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\Experiment\Data\Validation Data\NIRS';

fig_dir=fullfile(mainpath_out,'Fig_Validation_NIRS');
if ~isfolder(fig_dir)
    mkdir(fig_dir);
end

val_dir=fullfile(mainpath_out,'NIRS');
if ~isfolder(val_dir)
    mkdir(val_dir);
end
%% select ID number and cap
subjects_comb=[{'02','28','64'}];
subject_nirs_only=[{'19','21','43'}];%,'43','69','84'



%% LOAD COMBINED CAP DATA
% CONDITION 1 - NIRS ONLY
% CONDITION 2 - COMBINED CAP
for iSub = 1:size(subjects_comb,2)
    %% 1. Info
    sub = char(subjects_comb(iSub));
    
    layout=load(fullfile(mainpath_in_processed,'..','pre-processed',['sub-',sub],'3d','layout.mat'));
    layout_combined=layout.layout;
    cfg=[];
    cfg.layout=layout_combined;
    ft_layoutplot(cfg); 
    sub_dir=fullfile(mainpath_out,['sub-',sub]);
    if ~isfolder(sub_dir)
        mkdir(sub_dir);
    end
    
    
    switch sub
        case '28'
            rec='02';
        case '64'
            rec='01';
        case '02'
            rec='02';
    end
    
    
    task_label={'AutoDualCue','AutoSingleCue','NonAutoDualCue',...
        'NonAutoSingleCue','AutoDualNoCue','AutoSingleNoCue',...
        'NonAutoDualNoCue','NonAutoSingleNoCue'};
    %% 2. Load pre-processed data
    
    load(fullfile(mainpath_in_preprocessed,['sub-',sub],'nirs',['sub-',sub,'_rec-',rec,'_nirs_preprocessed.mat']));
    
    %% 3. Timelock analysis + baselinecorrection
    % these steps are necessary for the multisubject_averaging script!
    % a) timelock
    
    for task=1:8 %There are 8 tasks
        cfg               = [];
        cfg.trials        = find(nirs_preprocessed.trialinfo(:,1)==task); % Average the data for given task
        data_TL{task}     = ft_timelockanalysis(cfg, nirs_preprocessed);
    end
    save(fullfile(mainpath_in_processed,['sub-',sub],'nirs','data_TL.mat'), 'data_TL');

    % b) apply a baseline correction
    for task=1:8
        cfg                 = [];
        cfg.baseline        = [-10 0]; % define the amount of seconds you want to use for the baseline
        data_TL_blc{task}     = ft_timelockbaseline(cfg, data_TL{task});
        data_all{task}{iSub}=data_TL_blc{task};
    end
    save(fullfile(mainpath_in_processed,['sub-',sub],'nirs','data_TL_blc.mat'),'data_TL_blc');
    
    %% 4. Get the average for uncued tasks
    cap=2;
    load(fullfile(mainpath_in_processed,['sub-',sub],'nirs',['data_TL_blc.mat']));
    uncued_tasks=[6,8];
    % Use only finger auto & non auto single no cue - task 6 and 8
    for task=uncued_tasks
        cfg=[];
        data_all{cap}{iSub} = ft_timelockgrandaverage(cfg, data_TL_blc{task});
    end
        %     for task=uncued_tasks
%         data_all{cap}{iSub}{task}=data_TL_blc{task};
%     end
end
%% 5. Run time-lock and baseline if it hasn't beend one for the NIRS only cap
% these steps are necessary for the multisubject_averaging script!
% a) timelock
layout_nirs_only=load(fullfile(mainpath_in_nirsonly,'layout.mat'));
layout_nirs_only=layout_nirs_only.layout;
for s = 1:length(subject_nirs_only)
    cap=1;
    load(fullfile(mainpath_in_nirsonly, ['sub-',subject_nirs_only{s}], 'data_conc'));
    
    %There are 4 tasks (handauto, handnonauto, footauto, footnonauto)
    % We only want two tasks - 2: handnonauto and handauto
    
    for task=1:2 
        cfg               = [];
        cfg.trials        = find(data_conc.trialinfo(:,1)==task); % Average the data for given task
        data_TL_nirs{task} = ft_timelockanalysis(cfg, data_conc);
    end
    save(fullfile(mainpath_in_nirsonly, ['sub-',subject_nirs_only{s}], 'data_TL.mat'), 'data_TL_nirs');

    % b) apply a baseline correction
    for task=1:2
        cfg                 = [];
        cfg.baseline        = [-10 0]; % define the amount of seconds you want to use for the baseline
        data_TL_blc_nirs{task}     = ft_timelockbaseline(cfg, data_TL_nirs{task});
    end
    save(fullfile(mainpath_in_nirsonly, ['sub-',subject_nirs_only{s}], 'data_TL_blc.mat'), 'data_TL_blc_nirs');
    
    % c) get average of the two tasks
    for task=1:2
        cfg=[];
        data_all{cap}{s}=ft_timelockgrandaverage(cfg, data_TL_blc_nirs{task});
        %     load(fullfile(mainpath_in_nirsonly, ['sub-',subject_nirs_only{s}], 'data_TL_blc.mat'));
    end
    
end
 
%% 6. Average over all subjects -> for each cap seperately!
% cap 1: nirs only cap, cap 2: combined cap
for cap=1:2
cfg=[];
%For each cap, average all tasks
%Combined cap - 2 tasks => Auto/NonAuto, Single, Uncued
grandavg{cap}= ft_timelockgrandaverage(cfg, data_all{cap}{:});
end
save(fullfile(val_dir,'grandavg.mat'),'grandavg');



%% 7. Plot the data

% a) Separate O2Hb and HHb channels (only for cap 1, cap 2 as been done)
for cap=1:2
    cfg=[];
    cfg.channel='* [O2Hb]';
    data_TL_O2Hb{cap}=ft_selectdata(cfg, grandavg{cap});
    % and rename labels such that they have the same name as HHb channels
    for i=1:length(data_TL_O2Hb{cap}.label)
        tmp = strsplit(data_TL_O2Hb{cap}.label{i});
        data_TL_O2Hb{cap}.label{i}=tmp{1};
    end
    save(fullfile(val_dir,'data_TL_O2Hb.mat'),'data_TL_O2Hb');
    
    % The same for HHb channels
    cfg=[];
    cfg.channel='* [HHb]';
    data_TL_HHb{cap}=ft_preprocessing(cfg, grandavg{cap});
    for i=1:length(data_TL_HHb{cap}.label)
        tmp = strsplit(data_TL_HHb{cap}.label{i});
        data_TL_HHb{cap}.label{i}=tmp{1};
    end
    save(fullfile(val_dir,'data_TL_HHb.mat'),'data_TL_HHb');
end

% b) Topoplots for each cap
capname={'NIRS', 'COMBINED'};
cfg          = [];

cfg.marker   = 'labels';
cfg.ylim     = [-0.2 0.2];
cfg.xlim     = [5 10];
cfg.zlim     = cfg.ylim/2;
 % Choose the time window over which you want to average
for cap=1:2
    figure;
    cfg.title=capname{cap};
    if cap==1
        cfg.layout   = layout_nirs_only;
    else
        cfg.layout   = layout_combined;
    end
    ft_topoplotER(cfg, data_TL_O2Hb{cap});
    saveas(gcf,fullfile(fig_dir,['topoplot_',capname{cap},'.jpg']))
end

% c) Plot both on the lay-out

%Labels gotten from data
comb_channels=data_TL_O2Hb{2}.label;
% to_remove={'Rx1-Tx3', 'Rx4-Tx4', 'Rx3-Tx5', 'Rx12-Tx14'};
% to_remove_idx=find(contains(comb_channels,to_remove));
% comb_channels(to_remove_idx)=[];

nirs_channels=data_TL_O2Hb{1}.label;

%Check which channels are closer to each other
distances=zeros(length(nirs_channels),length(comb_channels));
for ii=1:length(nirs_channels)
    chan_nirs_idx=find(contains(layout_nirs_only.label,nirs_channels(ii)));
    chan_nirs_pos=layout_nirs_only.pos(chan_nirs_idx);
    for jj=1:1:length(comb_channels)
        chan_comb_idx=find(contains(layout_combined.label,comb_channels(jj)));
        chan_comb_pos=layout_combined.pos(chan_comb_idx);
        distances(ii,jj)=pdist([chan_nirs_pos;...
            chan_comb_pos],'euclidean');
    end
end
[min_dist,min_idx]=min(distances,[],2);
idx_channels=find(min_dist==distances);


% d. Plot most important channels

for ii=1:size(distances,1)
    nirs_only_channel = nirs_channels(ii);
    comb_channel      = comb_channels(min_idx(ii));
    
    cfg                   = [];
    cfg.title             = ['Combined - ',comb_channel...
        ,' NIRS - ', nirs_only_channel];
    cfg.showlabels        = 'yes';
    cfg.interactive       = 'yes'; % this allows to select a subplot and interact with it
    cfg.linecolor        = 'rbrbmcmc'; % O2Hb is showed in red (finger) and magenta (foot), HHb in blue (finger) and cyan (foot)
    cfg.linestyle = {'--', '--', '-', '-'}; % fingerauto is dashed line, fingernonauto is solid line, footauto is dotted line and footnonauto is a dotted stars line
    cfg.comment = 'NIRS is dashed line, COMBINED is solid line';
    %cfg.ylim = [-0.440 0.540];
    cfg.channel = [nirs_only_channel, nirs_only_channel, comb_channel, comb_channel];
%     cfg.figure = 'gcf';
%     subplot(4,3,ii)
    ft_singleplotER(cfg, data_TL_O2Hb{1}, data_TL_HHb{1}, data_TL_O2Hb{2}, data_TL_HHb{2});
    saveas(gcf, ['Fig_Validation_NIRS/Combined_',comb_channel{1},'NIRS_', nirs_only_channel{1},'.jpg']);
end

% for ii=12:24
%     nirs_only_channel = layout_nirs_only.label(ii);
%     comb_channel      = layout_combined.label(min_idx(ii));
%     
%     cfg                   = [];
%     cfg.title             = ['Combined - ',comb_channel...
%         ,'    NIRS - ', nirs_only_channel];
%     cfg.showlabels        = 'yes';
%     cfg.interactive       = 'yes'; % this allows to select a subplot and interact with it
%     cfg.linecolor        = 'rbrbmcmc'; % O2Hb is showed in red (finger) and magenta (foot), HHb in blue (finger) and cyan (foot)
%     cfg.linestyle = {'--', '--', '-', '-'}; % fingerauto is dashed line, fingernonauto is solid line, footauto is dotted line and footnonauto is a dotted stars line
%     cfg.comment = 'NIRS is dashed line, COMBINED is solid line';
%     %cfg.ylim = [-0.440 0.540];
%     cfg.channel = [ii, ii, min_idx(ii), min_idx(ii)];
%     cfg.figure = 'gcf';
%     subplot(4,3,ii/12)
%     ft_singleplotER(cfg, data_TL_O2Hb{1}, data_TL_HHb{1}, data_TL_O2Hb{2}, data_TL_HHb{2});
% end

%e) Plot all channels
cfg                   = [];
cfg.showlabels        = 'yes';
cfg.layout            = layout_combined;
cfg.interactive       = 'yes'; % this allows to select a subplot and interact with it
cfg.linecolor        = 'rbrbmcmc'; % O2Hb is showed in red (finger) and magenta (foot), HHb in blue (finger) and cyan (foot)
cfg.linestyle = {'--', '--', '-', '-'}; % fingerauto is dashed line, fingernonauto is solid line, footauto is dotted line and footnonauto is a dotted stars line
cfg.comment = 'NIRS is dashed line, COMBINED is solid line';
%cfg.ylim = [-0.440 0.540];
ft_multiplotER(cfg, data_TL_O2Hb{1}, data_TL_HHb{1}, data_TL_O2Hb{2}, data_TL_HHb{2});

% f) Plot for each cap seperately

%taskshort={'complex', 'stroop'};
for cap=1:2
    cfg                   = [];
    cfg.showlabels        = 'yes';
    if cap==1,  cfg.layout            = layout_nirs_only;
    else,       cfg.layout            = layout_combined; 
    end
    cfg.showoutline       = 'yes';
    cfg.interactive       = 'yes'; % this allows to select a subplot and interact with it
    cfg.linecolor         = 'rb';% O2Hb is showed in red, HHb in blue
    cfg.showcomment       = 'yes';
    cfg.title           = capname{cap};
    % cfg.ylim = [-0.440 0.540];
%     cfg.colorgroups=contains(data_TL_blc{task}.label, '[O2Hb]')+2*contains(data_TL_blc{task}.label, '[HHb]');
    ft_multiplotER(cfg, data_TL_O2Hb{cap}, data_TL_HHb{cap})
    saveas(gcf, ['Fig_Validation_NIRS/',char(capname(cap)),'_timelock.jpg']);
end


%% 7. Statistical testing
cap_names={'nirs only',  'combined'};
for i=1:2 % loop over the 4 conditions
  [stat_O2Hb, stat_HHb] = statistics_withinsubjects(grandavg{i}, 'grandavg', layout, i, cap_names{i});
end
%% 8. Statistical analysis (not finished yet! still trying things out) 


% T-TEST
% define the parameters for the statistical comparison
cfg = [];
cfg.channel     = 'Rx1-Tx2 [HHb]', 'Rx1-Tx2 [O2Hb]';
cfg.latency     = 'all';
cfg.avgovertime = 'yes';
cfg.parameter   = 'avg';
cfg.method      = 'analytic';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'no';

Nsub = 2;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

stat = ft_timelockstatistics(cfg, grandavg{1}(:), grandavg{2}(:));   % don't forget the {:}!

% t-test with matlab function
chan = 7;
time = [-10 20];

% find the time points for the effect of interest in the grand average data
timesel_nirs_only = find(grandavg{1}.time >= time(1) & grandavg{1}.time <= time(2));
timesel_combined  = find(grandavg{2}.time >= time(1) & grandavg{2}.time <= time(2));


% select the individual subject data from the time points and calculate the mean
for isub = 1:17
    values_nirs_only(isub) = mean(data_all{1}{isub}.avg(chan,timesel_nirs_only));
    values_combined(isub)  = mean(data_all{2}{isub}.avg(chan,timesel_combined));
end

% plot to see the effect in each subject
M = [values_nirs_only(:) values_combined(:)];
figure; plot(M', 'o-'); xlim([0.5 2.5])
%legend({'subj1', 'subj2', 'subj3', 'subj4', 'subj5', 'subj6', ...
        %'subj7', 'subj8', 'subj9', 'subj10'}, 'location', 'EastOutside');

subtraction = values_nirs_only - values_combined;
[h,p,ci,stats] = ttest(subtraction, 0, 0.05) % H0: mean = 0, alpha 0.05


%loop over channels
time = [-10 20];
timesel_nirs_only = find(grandavg{1}.time >= time(1) & grandavg{1}.time <= time(2));
timesel_combined  = find(grandavg{2}.time >= time(1) & grandavg{2}.time <= time(2));
clear h p

%FAminFNA = zeros(1,17);
FAminFOA = zeros(1,length(subject_comb));

for iChan = 1:2
    for isub = 1:length(subject_comb)
        FAminFOA(isub) = ...
            mean(data_all{1}{isub}.avg(iChan,timesel_nirs_only)) - ...
            mean(data_all{3}{isub}.avg(iChan,timesel_footauto));
    end

    [h(iChan), p(iChan)] = ttest(FAminFOA, 0, 0.05 ) % test each channel separately
end







