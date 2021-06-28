function [stat_O2Hb, stat_HHb, h] = statistics_withinsubjects(data, data_name, layout, condition, condition_name, sub, rec, varargin)
% This function calculates the statistics of the fNIRS activity vs baseline
% for O2Hb and HHb and shows the significant channels in a figure. 
% based on http://www.fieldtriptoolbox.org/tutorial/eventrelatedstatistics/
%
% Use as
%   [stat_O2Hb, stat_HHb] = statistics_withinsubjects(data, data_name, layout, condition, varargin)
%
% INPUT:
%       data         = fnirs data
%       data_name    = string with name of the data input (e.g. 'data_raw')
%       layout       = layout on which the graphs are plotted
%       condition    = number that corresponds to the task which you would
%       like to test (e.g. 1). Only one condition is allowed for this
%       function.
%       condition_name= name of the condition (string)
% Additional options can be specified in key-value pairs and can be:
%       'channel'    = string with channel of which you would like to plot
%       the individual trials. (default = {} --> no trials are plotted)
%       'baseline'   = vector with timepoints that should be used as baseline (default = [-10 0])
%       'xlimit'       = [xmin xmax]; (default = [-10 20])  ! does not work
%       yet
%       'ylimit'       = [ymin ymax]; (default = [-3 5]) ! does not work yet
%
% OUTPUT
%       stat_O2Hb    = fieldtrip stat structure including the statistics
%       for all the O2Hb channels
%       stat_HHb     = fieldtrip stat structure including the statistics
%       for all the HHb channels

%% Get the options
chan = ft_getopt(varargin, 'channel', {});
baseline=ft_getopt(varargin, 'baseline', [-10 0]);
xlimit=ft_getopt(varargin, 'xlimit', [-10 20]);
ylimit=ft_getopt(varargin, 'ylim', [-3 5]);
taskstatisticsO2Hb = {['sub-',sub,'_rec-',rec,'_dualcued_statisticsO2Hb'], ['sub-',sub,'_rec-',rec,'_singlecued_statisticsO2Hb'], ['sub-',sub,'_rec-',rec,'_dualuncued_statisticsO2Hb'], ['sub-',sub,'_rec-',rec,'_singleuncued_statisticsO2Hb']};
taskstatisticsHHb = {['sub-',sub,'_rec-',rec,'_dualcued_statisticsHHb'], ['sub-',sub,'_rec-',rec,'_singlecued_statisticsHHb'], ['sub-',sub,'_rec-',rec,'_dualuncued_statisticsHHb'], ['sub-',sub,'_rec-',rec,'_singleuncued_statisticsHHb']};
  

%% apply baseline correction 
cfg                 = [];
cfg.baseline        = baseline; % define the amount of seconds you want to use for the baseline
cfg.parameter = 'trial';
data_blc = ft_timelockbaseline(cfg, data);

% Collect data for the O2Hb and HHb channel during each trial of the
% complex condition
cfg=[];
cfg.channel='*[O2Hb]';
cfg.trials=[data_blc.trialinfo(:,1)]'==condition;
data_complex_O2Hb=ft_selectdata(cfg, data_blc);

cfg=[];
cfg.trials=[data_blc.trialinfo(:,1)]'==condition;
cfg.channel='*[HHb]';
data_complex_HHb=ft_selectdata(cfg, data_blc);

%% Plot all trials for the given condition (complex in this case)
if ~isempty(chan)
    time = [0 10];
    channel=find(startsWith(data_complex_O2Hb.label, chan));
    
    % Plot trials for complex task
    figure;
    for i = 1:length(data_complex_O2Hb.trialinfo)
        subplot(3,4,i)
        % use the rectangle to indicate the time range used later
        rectangle('Position',[time(1) ylimit(1) (time(2)-time(1)) ylimit(2)-ylimit(1)],'FaceColor',[0.7 0.7 0.7]);
        hold on;
        % plot the lines in front of the rectangle
        plot(data_complex_O2Hb.time, permute(data_complex_O2Hb.trial(i, channel,:), [3 1 2]), 'r');
        plot(data_complex_HHb.time, permute(data_complex_HHb.trial(i, channel,:), [3 1 2]), 'b');
        title(strcat('trial ',num2str(i)))
        ylim(ylimit)
        xlim(xlimit)
    end
    subplot(3,4,11);
    text(0.5, 0.7, sprintf('%s in %s', condition_name, chan), 'color', 'k'); text(0.5,0.5,'oxyHb','color','r') ;text(0.5,0.3,'deoxyHb','color','b')
    axis off;
end
%% t-tests with fieldtrip
    n_trials=length(data_complex_O2Hb.trialinfo);
    % for oxy channels
    cfg = [];
    cfg.latency     = [5 10];
    cfg.avgovertime = 'yes';
    cfg.parameter   = 'trial';
    cfg.method      = 'stats';
    cfg.statistic   = 'ttest';
    cfg.alpha       = 0.05;
    % cfg.alpha       = 0.05/32; % bonferroni correction for mulitple comparison (short channels are not included)
    cfg.tail        = 1;
    cfg.correctm = 'yes';
    cfg.design(1, 1:n_trials) = ones(1,n_trials);
    cfg.ivar =1;
    stat_O2Hb = ft_timelockstatistics(cfg,data_complex_O2Hb); 
    % plot
    data_complex_O2Hb_lay=data_complex_O2Hb;
    for i=1:length(data_complex_O2Hb.label)
        split_label=strsplit(data_complex_O2Hb.label{i});
        data_complex_O2Hb_lay.label{i}=split_label{1};
    end
    cfg = [];
    cfg.style     = 'blank';
    cfg.layout    = layout;
    cfg.highlight = 'on';
    cfg.highlightchannel = find(stat_O2Hb.mask);
    cfg.comment   = 'no';
    h=figure;
    ft_topoplotER(cfg, data_complex_O2Hb_lay); title(sprintf('significant channels (p=0.05) for %s vs baseline [O2Hb]', condition_name))
    saveas(h,taskstatisticsO2Hb{condition},'png'); 
    
    % for deoxy channels
    cfg = [];
    cfg.latency     = [5 10];
    cfg.avgovertime = 'yes';
    cfg.parameter   = 'trial';
    cfg.method      = 'stats';
    cfg.statistic   = 'ttest';
    cfg.alpha       = 0.05;
    % cfg.alpha       = 0.05/32; % bonferroni correction for mulitple comparison (short channels are not included)
    cfg.tail        = -1;
    cfg.correctm = 'yes';
    cfg.design(1, 1:n_trials) = ones(1,n_trials);
    cfg.ivar =1;
    stat_HHb = ft_timelockstatistics(cfg,data_complex_HHb); 
    % plot
    data_complex_HHb_lay=data_complex_HHb;
    for i=1:length(data_complex_HHb.label)
        split_label=strsplit(data_complex_HHb.label{i});
        data_complex_HHb_lay.label{i}=split_label{1};
    end
    cfg = [];
    cfg.style     = 'blank';
    cfg.layout    = layout;
    cfg.highlight = 'on';
    cfg.highlightchannel = find(stat_HHb.mask);
    cfg.comment   = 'no';
    h=figure;
    ft_topoplotER(cfg, data_complex_HHb_lay); title(sprintf('significant channels (p=0.05) for %s vs baseline [HHb]', condition_name))
	saveas(h,taskstatisticsHHb{condition},'png'); 

end



