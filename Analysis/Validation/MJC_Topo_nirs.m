clear, close all
%% Settings

addpath (fullfile(pwd,'..')) %add the path with analysis scripts
% add path to correct folders and open eeglab
laptop = 'laptopJoao';
[mainpath_in, ~, eeglab_path] = addFolders(laptop);
mainpath_in=fullfile(mainpath_in,'pre-processed');
mainpath_out = 'C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\Experiment\Data\Exp\processed';

% select ID number and cap
subjects=[{'02','64','28','76'}];
subject_nirs_only=[{'03','04','10','11'}];

%%
for iSub = 3:size(subjects,2)
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
    
    load(fullfile(mainpath_in,['sub-',sub],'nirs',['sub-',sub,'_rec-',rec,'_nirs_preprocessed.mat']));

    
    
    %% 3. Visualize data 
    h=multiplot_condition(nirs_preprocessed, layout, [1:8], task_label, ....
        'baseline', [-10 0], 'trials', true, 'topoplot', 'yes');

    %% 4. Timelock analysis + baselinecorrection
    % these steps are necessary for the multisubject_averaging script!
    % a) timelock
    
    for task=1:8 %There are 8 tasks
        cfg               = [];
        cfg.trials        = find(nirs_preprocessed.trialinfo(:,1)==task); % Average the data for given task
        data_TL{task}     = ft_timelockanalysis(cfg, nirs_preprocessed);
    end
    save(fullfile(sub_dir,'nirs','data_TL.mat'), 'data_TL');

    % b) apply a baseline correction
    for task=1:8
        cfg                 = [];
        cfg.baseline        = [-10 0]; % define the amount of seconds you want to use for the baseline
        data_TL_blc{task}     = ft_timelockbaseline(cfg, data_TL{task});
        data_all{task}{iSub}=data_TL_blc{task};
    end
    save(fullfile(sub_dir,'nirs','data_TL_blc.mat'),'data_TL_blc');
    
    %% 5. Separate HbO activity only
    data_labels=nirs_preprocessed.label;
    for task=1:8
        cfg=[];
        cfg.channel='*[O2Hb]';
        data_O2Hb{task}=ft_selectdata(cfg, data_TL_blc{task});
        data_O2Hb{task}.label=data_labels(contains(data_TL_blc{task}.label, '[O2Hb]'));
        for ii=1:length(data_O2Hb{task}.label)
            label=data_O2Hb{task}.label{ii};
            label=label(1:end-7);
            data_O2Hb{task}.label{ii}=label;
        end
    end
    

    %% 6. Topoplot - per subject
    cfg          = [];
    cfg.layout   = layout;
    cfg.marker   = 'labels';
     % Choose the time window over which you want to average
    for task=1:8
        cfg.xlim     = [5 10];
        cfg.zlim     = [-0.1 0.1];
        figure;
        title(task_label{task})
        ft_topoplotER(cfg, data_O2Hb{task});
        saveas(gcf,fullfile(fig_dir,['sub-',sub,'_topoplot_',task_label{task},'.jpg']))
    end
    
    
end


%% TOPOPLOT FOR ALL SUBJECT
for task=1:8
    cfg=[];
    grandavg{task}= ft_timelockgrandaverage(cfg, data_all{task}{:});
end

%% 5. Plot the data


% d) Separate O2Hb and HHb channels (only for task 1, task 2 as been done)
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
cfg          = [];
cfg.layout   = layout;
cfg.marker   = 'labels';
cfg.ylim     = [-0.2 0.2];
cfg.xlim     = [5 10];
cfg.zlim     = cfg.ylim/2;
 % Choose the time window over which you want to average
for task=1:8
    figure;
    title(task_label{task})
    ft_topoplotER(cfg, data_TL_O2Hb{task});
    saveas(gcf,fullfile(pwd,'Fig_NIRS_topo_all_cond/',['topoplot_',task_label{task},'.jpg']))
end

%% f) Topoplots for total avg
for task=1:8
    cfg=[];
    data_all= ft_timelockgrandaverage(cfg, data_TL_O2Hb{task});
end

cfg          = [];
cfg.layout   = layout;
cfg.marker   = 'labels';
cfg.ylim     = [-0.2 0.2];
cfg.xlim     = [5 10];
cfg.zlim     = cfg.ylim/2;
figure;
title('All task averaged')
ft_topoplotER(cfg, data_TL_O2Hb{task});

saveas(gcf,fullfile(pwd,'Fig_NIRS_topo_all_cond/',['topoplot_alltask_.jpg']))


