clear, close all
%% Settings

addpath (fullfile(pwd,'..')) %add the path with analysis scripts
% add path to correct folders and open eeglab
laptop = 'laptopJoao';
[mainpath_in, ~, eeglab_path] = addFolders(laptop);
mainpath_in=fullfile(mainpath_in,'pre-processed');
mainpath_out = 'C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\Experiment\Data\Exp\processed';

% select ID number and cap
subject=[{'28'}];
subject_nirs_only=[{'03','04','10','11','12',}];


for iSub = 1:size(subject,2)
    %% 1. Info
    sub = char(subject(iSub));
    
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
    end
    
    
    task_label={'AutoDualCue','AutoSingleCue','NonAutoDualCue',...
        'NonAutoSingleCue','AutoDualNoCue','AutoSingleNoCue',...
        'NonAutoDualNoCue','NonAutoSingleNoCue'};
    %% 2. Load pre-processed data
    
    load(fullfile(mainpath_in,['sub-',sub],'nirs',['sub-',sub,'_rec-',rec,'_nirs_preprocessed.mat']));
    
%     %% 3. Visualize data
%     h=multiplot_condition(nirs_preprocessed, layout, [1:8], task_label, 'baseline', ....
%         [-10 0], 'trials', false, 'topoplot', 'yes');

    %% 4. Timelock analysis + baselinecorrection
    % these steps are necessary for the multisubject_averaging script!
    % a) timelock
    for task=1:8 %There are 4 tasks (handauto, handnonauto, footauto, footnonauto)
        cfg               = [];
        cfg.trials        = find(nirs_preprocessed.trialinfo(:,1)==task); % Average the data for given task
        data_TL{task} = ft_timelockanalysis(cfg, nirs_preprocessed);
    end
    save(fullfile(sub_dir,'nirs','data_TL.mat'), 'data_TL');

    % b) apply a baseline correction
    for task=1:8
        cfg                 = [];
        cfg.baseline        = [-10 0]; % define the amount of seconds you want to use for the baseline
        data_TL_blc{task}     = ft_timelockbaseline(cfg, data_TL{task});
    end
    save(fullfile(sub_dir,'nirs','data_TL_blc.mat'),'data_TL_blc');
    
    
%     % See which labels coincide in the layout
% 
%     validChannel=contains(data_TL_blc{1}.label,layout.label);
% 
%     validLabels=data_TL_blc{1}.label(validChannel);
%     for i=1:size(validLabels,1)
%         x=strfind(validLabels{i},' ');
%         validLabels{i}=validLabels{i}(1:x-1);
%     end
%     validLabels=unique(validLabels);
% 
%     channels_to_keep=validLabels(ismember(validLabels,layout.label));
% 
% 
% 
%     % Change layout channels to the ones in the data
%     cfg=[];
%     cfg.channel=channels_to_keep;
%     cfg.layout=layout;
%     cfg.opto=layout.cfg.opto;
%     layoutNIRS_new=ft_prepare_layout(cfg);
%     layoutNIRS_new.cfg.opto.label=layoutNIRS_new.cfg.opto.label(ismember(validLabels,layoutNIRS_new.cfg.opto.label));
%     layoutNIRS_new.cfg.opto.chanpos=layoutNIRS_new.cfg.opto.chanpos(ismember(validLabels,layoutNIRS_new.cfg.opto.label));
    %% 5. Topoplot
    for task=1:8
        cfg          = [];
        cfg.layout   = layout;
        cfg.marker   = 'labels';
        cfg.xlim     = [5 10]; % Choose the time window over which you want to average
        figure;
        ft_topoplotER(cfg, data_TL_blc{task});
        saveas(gcf,fullfile(fig_dir,['sub-',sub,'_topoplot_',task_label(task)]))
    end
    
end