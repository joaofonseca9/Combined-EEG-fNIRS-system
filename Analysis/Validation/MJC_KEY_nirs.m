clear, close all
%% Settings

addpath (fullfile(pwd,'..')) %add the path with analysis scripts
% add path to correct folders and open eeglab
laptop = 'laptopJoao';
[mainpath_in, ~, eeglab_path] = addFolders(laptop);
mainpath_in=fullfile(mainpath_in,'pre-processed');
mainpath_out = 'C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\Experiment\Data\Exp\processed';

% select ID number and cap
subjects=[{'02','28','64','76'}];
subject_nirs_only=[{'03','04','10','11'}];

%%
for iSub = 1:size(subjects,2)
    %% 1. Info
    sub = char(subjects(iSub));
    
    load(fullfile(mainpath_in,['sub-',sub],'3d','layout.mat'));
    sub_dir=fullfile(mainpath_out,['sub-',sub]);
    if ~isfolder(sub_dir)
        mkdir(sub_dir);
    end
    
    fig_dir=fullfile(sub_dir,'Fig_Topoplot');
    if ~isfolder(fig_dir)
        mkdir(fig_dir);
    end
    
    switch sub
        case '28'
            rec='02';
        case '64'
            rec='01';
        case '02'
            rec='02';
        case '76'
            rec='01';
    end
    
    
    
    task_label={'AutoDualCue','AutoSingleCue','NonAutoDualCue',...
        'NonAutoSingleCue','AutoDualNoCue','AutoSingleNoCue',...
        'NonAutoDualNoCue','NonAutoSingleNoCue'};
    %% 2. Load pre-processed data
    
%     load(fullfile(mainpath_in,['sub-',sub],'nirs',['sub-',sub,'_rec-',rec,'_nirs_preprocessed.mat']));
    load(fullfile(mainpath_in,['sub-',sub],'nirs',['sub-',sub,'_rec-',rec,'_nirs_epoch.mat']));

%     load(fullfile(mainpath_in,['sub-',sub],'nirs',['sub-',sub,'_rec-',rec,'_nirs_lpf.mat']));
%     nirs_preprocessed.trialinfo=nirs_lpf.trialinfo;
    
    %% 3. Visualize data 
%     h=multiplot_condition(nirs_preprocessed, layout, [1:8], task_label, ....
%         'baseline', [-10 0], 'trials', true, 'topoplot', 'yes');

    %% 4. Timelock analysis + baselinecorrection
    % these steps are necessary for the multisubject_averaging script!
    % a) timelock
    
    for task=1:8 %There are 8 tasks
        cfg               = [];
        cfg.trials        = find(nirs_epoch.trialinfo(:,1)==task); % Average the data for given task
        cfg.nanmean='yes';
        data_TL{task}     = ft_timelockanalysis(cfg, nirs_epoch);
    end
    %save(fullfile(sub_dir,'nirs','data_TL.mat'), 'data_TL');
 %% Add removed channels
labels = {'Rx2-Tx4 [O2Hb]'; 'Rx2-Tx4 [HHb]';...
    'Rx2-Tx3 [O2Hb]'; 'Rx2-Tx3 [HHb]';...
    'Rx1-Tx3 [O2Hb]'; 'Rx1-Tx3 [HHb]';...
    'Rx1-Tx2 [O2Hb]'; 'Rx1-Tx2 [HHb]';...
    'Rx4-Tx4 [O2Hb]'; 'Rx4-Tx4 [HHb]';...
    'Rx4-Tx5 [O2Hb]'; 'Rx4-Tx5 [HHb]';...
    'Rx3-Tx3 [O2Hb]'; 'Rx3-Tx3 [HHb]';...
    'Rx3-Tx2 [O2Hb]'; 'Rx3-Tx2 [HHb]';...
    'Rx3-Tx5 [O2Hb]'; 'Rx3-Tx5 [HHb]';...
    'Rx6-Tx9 [O2Hb]'; 'Rx6-Tx9 [HHb]';...
    'Rx5-Tx8 [O2Hb]'; 'Rx5-Tx8 [HHb]';...
    'Rx5-Tx7 [O2Hb]'; 'Rx5-Tx7 [HHb]';...
    'Rx8-Tx9 [O2Hb]'; 'Rx8-Tx9 [HHb]';...
    'Rx8-Tx10 [O2Hb]'; 'Rx8-Tx10 [HHb]';...
    'Rx7-Tx8 [O2Hb]'; 'Rx7-Tx8 [HHb]';...
    'Rx7-Tx7 [O2Hb]'; 'Rx7-Tx7 [HHb]';...
    'Rx9-Tx12 [O2Hb]'; 'Rx9-Tx12 [HHb]';...
    'Rx9-Tx13 [O2Hb]'; 'Rx9-Tx13 [HHb]';...
    'Rx11-Tx12 [O2Hb]'; 'Rx11-Tx12 [HHb]';...
    'Rx11-Tx13 [O2Hb]'; 'Rx11-Tx13 [HHb]';...
    'Rx10-Tx14 [O2Hb]'; 'Rx10-Tx14 [HHb]';...
    'Rx12-Tx14 [O2Hb]'; 'Rx12-Tx14 [HHb]'};
    % b) apply a baseline correction
    for task=1:8
        cfg                 = [];
        cfg.baseline        = [-10 0]; % define the amount of seconds you want to use for the baseline
        cfg.nanmean='yes';
        data_TL_blc{task}     = ft_timelockbaseline(cfg, data_TL{task});
        data_all{task}{iSub}=data_TL_blc{task};
    end
    %save(fullfile(sub_dir,'nirs','data_TL_blc.mat'),'data_TL_blc');
    
%     %% 5. Separate HbO activity only
%     data_labels=nirs_preprocessed.label;
%     for task=1:8
%         cfg=[];
%         cfg.channel='*[O2Hb]';
%         data_O2Hb{task}=ft_selectdata(cfg, data_TL_blc{task});
%         data_O2Hb{task}.label=data_labels(tasktains(data_TL_blc{task}.label, '[O2Hb]'));
%         for ii=1:length(data_O2Hb{task}.label)
%             label=data_O2Hb{task}.label{ii};
%             label=label(1:end-7);
%             data_O2Hb{task}.label{ii}=label;
%         end
%     end
    

%     %% 6. Topoplot - per subject
%     cfg          = [];
%     cfg.layout   = layout;
%     cfg.marker   = 'labels';
%      % Choose the time window over which you want to average
%     for task=1:8
%         cfg.xlim     = [5 10];
%         cfg.zlim     = [-0.1 0.1];
%         figure;
%         title(task_label{task})
%         ft_topoplotER(cfg, data_O2Hb{task});
%         %saveas(gcf,fullfile(fig_dir,['sub-',sub,'_topoplot_',task_label{task},'.jpg']))
%     end
    
    
end


%% TOPOPLOT FOR ALL SUBJECT
for task=1:8
    cfg=[];
    cfg.nanmean='yes';
    grandavg{task}= ft_timelockgrandaverage(cfg, data_all{task}{:});
%     grandavg{task}.dimord='chan_time';
end

%% 5. Plot the data


% d) Separate O2Hb and HHb channels 
for task=1:8
    cfg=[];
    cfg.channel='* [O2Hb]';
    data_TL_O2Hb{task}=ft_selectdata(cfg, grandavg{task});
    % and rename labels such that they have the same name as HHb channels
    for i=1:length(data_TL_O2Hb{task}.label)
        tmp = strsplit(data_TL_O2Hb{task}.label{i});
        data_TL_O2Hb{task}.label{i}=tmp{1};
    end
    save(fullfile('NIRS/','data_TL_O2Hb.mat'),'data_TL_O2Hb');
    
    % The same for HHb channels
    cfg=[];
    cfg.channel='* [HHb]';
    data_TL_HHb{task}=ft_preprocessing(cfg, grandavg{task});
    for i=1:length(data_TL_HHb{task}.label)
        tmp = strsplit(data_TL_HHb{task}.label{i});
        data_TL_HHb{task}.label{i}=tmp{1};
    end
    save(fullfile('NIRS/','data_TL_HHb.mat'),'data_TL_HHb');
end

%% e) Topoplots for each task
task_label={'AutoDualCue','AutoSingleCue','NonAutoDualCue',...
        'NonAutoSingleCue','AutoDualNoCue','AutoSingleNoCue',...
        'NonAutoDualNoCue','NonAutoSingleNoCue'};
cfg          = []; cfg.layout   = layout; cfg.marker   = 'labels';
cfg.ylim     = [-0.2 0.2]; cfg.xlim     = [5 10]; 
 % Choose the time window over which you want to average
for task=1:8
    cfg.title=task_label{task}; 
    ft_topoplotER(cfg,data_TL_O2Hb{task});
    saveas(gcf,fullfile(pwd,'Fig_NIRS_topo_all_taskd/',['topoplot_',task_label{task},'.jpg']))
end

%% f) Topoplots for total avg

cfg=[];
grand_avg_HbO2= ft_timelockgrandaverage(cfg, data_TL_O2Hb{:});
grand_avg_HHb= ft_timelockgrandaverage(cfg, data_TL_HHb{:});

cfg          = [];
cfg.layout   = layout;
cfg.marker   = 'labels';
cfg.ylim     = [-0.2 0.2];
cfg.xlim     = [5 10];
cfg.zlim     = cfg.ylim/2;
figure;
cfg.title='Grand Average - HbO2 ';
ft_topoplotER(cfg, grand_avg_HbO2);

saveas(gcf,fullfile(pwd,'Fig_NIRS_topo_all_taskd/',['topoplot_alltask_.jpg']))



%% Hemoglobin activity
cfg                   = [];
cfg.showlabels        = 'yes';
cfg.layout   = layout;
cfg.interactive       = 'yes'; % this allows to select a subplot and interact with it
cfg.linecolor        = 'rb'; % O2Hb is showed in red (finger) and magenta (foot), HHb in blue (finger) and cyan (foot)
%cfg.ylim = [-0.440 0.540];
ft_multiplotER(cfg, grand_avg_HbO2, grand_avg_HHb);

%% Single plot

cfg                   = [];
cfg.title             = 'Grand Average';
cfg.showlabels        = 'yes';
cfg.layout   = layout;
cfg.plotstderr = 'yes';
cfg.interactive       = 'yes'; % this allows to select a subplot and interact with it
cfg.linecolor        = 'rb'; % O2Hb is showed in red (finger) and magenta (foot), HHb in blue (finger) and cyan (foot)
%cfg.ylim = [-0.440 0.540];
ft_singleplotER(cfg, grand_avg_HbO2, grand_avg_HHb);
saveas(gcf,fullfile(pwd,'Fig_NIRS_topo_all_taskd/',['hemoglobin_alltask_.jpg']))

%%
% stdplot(grand_avg_HbO2,grand_avg_HHb, cfg.title)
%% Plot per region
cfg                   = [];
cfg.showlabels        = 'yes';
cfg.layout   = layout;
cfg.channel = {'Rx5-Tx7', 'Rx5-Tx8', 'Rx7-Tx7', 'Rx7-Tx8'};
cfg.interactive       = 'yes'; % this allows to select a subplot and interact with it
cfg.linecolor        = 'rb'; % O2Hb is showed in red (finger) and magenta (foot), HHb in blue (finger) and cyan (foot)
% cfg.ylim = [-0.0440 0.0540];

cfg.title             = 'Left dlPFC';
ft_singleplotER(cfg, grand_avg_HbO2, grand_avg_HHb);
saveas(gcf,fullfile(pwd,'Fig_NIRS_topo_all_taskd/',[cfg.title ,'_hemoglobin_alltask_.jpg']))

cfg.title             = 'Right dlPFC';
cfg.channel = {'Rx9-Tx12', 'Rx9-Tx13', 'Rx11-Tx12', 'Rx11-Tx13'};
ft_singleplotER(cfg, grand_avg_HbO2, grand_avg_HHb);
saveas(gcf,fullfile(pwd,'Fig_NIRS_topo_all_taskd/',[cfg.title ,'_hemoglobin_alltask_.jpg']))

cfg.title             = 'Left PPC';
cfg.channel = {'Rx6-Tx9', 'Rx8-Tx9', 'Rx8-Tx10'};
ft_singleplotER(cfg, grand_avg_HbO2, grand_avg_HHb);
saveas(gcf,fullfile(pwd,'Fig_NIRS_topo_all_taskd/',[cfg.title ,'_hemoglobin_alltask_.jpg']))

% 
cfg.title             = 'Right PPC';
cfg.channel = {'Rx10-Tx14', 'Rx12-Tx14'};
ft_singleplotER(cfg, grand_avg_HbO2, grand_avg_HHb);
saveas(gcf,fullfile(pwd,'Fig_NIRS_topo_all_taskd/',[cfg.title ,'_hemoglobin_alltask_.jpg']))

cfg.title             = 'SMA_M1';
cfg.channel = {'Rx4-Tx4', 'Rx4-Tx5', 'Rx1-Tx2', 'Rx1-Tx3', 'Rx3-Tx2'...
    , 'Rx3-Tx3', 'Rx3-Tx5', 'Rx2-Tx4', 'Rx2-Tx3'};
ft_singleplotER(cfg, grand_avg_HbO2, grand_avg_HHb);
saveas(gcf,fullfile(pwd,'Fig_NIRS_topo_all_taskd/',[cfg.title ,'_hemoglobin_alltask_.jpg']))


cfg.title             = 'SMA';
cfg.channel = {'Rx4-Tx4', 'Rx3-Tx5'};
ft_singleplotER(cfg, grand_avg_HbO2, grand_avg_HHb);
saveas(gcf,fullfile(pwd,'Fig_NIRS_topo_all_taskd/',[cfg.title ,'_hemoglobin_alltask_.jpg']))

cfg.title             = 'M1';
cfg.channel = {'Rx4-Tx5', 'Rx1-Tx2', 'Rx1-Tx3', 'Rx3-Tx2'...
    , 'Rx3-Tx3', 'Rx2-Tx4', 'Rx2-Tx3'};
ft_singleplotER(cfg, grand_avg_HbO2, grand_avg_HHb);
saveas(gcf,fullfile(pwd,'Fig_NIRS_topo_all_taskd/',[cfg.title ,'_hemoglobin_alltask_.jpg']))

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function stdplot(grand_avg_HbO2,grand_avg_HHb,figtitle)
time=grand_avg_HbO2.time;
avg1=nanmean(grand_avg_HbO2.avg,1) + 1*std(grand_avg_HbO2.avg,[],1);
avg2=nanmean(grand_avg_HbO2.avg,1) - 1*std(grand_avg_HbO2.avg,[],1);
figure,title(figtitle)
plot(time,nanmean(grand_avg_HbO2.avg,1),'r');hold on; h1 = fill([time,fliplr(time)], [avg1,fliplr(avg2)],'r','LineStyle','none');
set(h1,'FaceAlpha',0.4);
avg1=nanmean(grand_avg_HHb.avg,1) + 1*std(grand_avg_HHb.avg,[],1);
avg2=nanmean(grand_avg_HHb.avg,1) - 1*std(grand_avg_HHb.avg,[],1);
plot(time,nanmean(grand_avg_HHb.avg,1),'b');hold on; h2 = fill([time,fliplr(time)], [avg1,fliplr(avg2)],'b','LineStyle','none');
set(h2,'FaceAlpha',0.4);
end
