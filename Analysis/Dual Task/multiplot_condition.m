function [h]=multiplot_condition(data, layout, conditions, condition_names, varargin)
% This function is meant to create multiplots with the average for each
% conidtion plotted in the layout. You can opt to plot all individual
% trials into the plot ('trials'=true)
%
% Use as
%   [h] = multiplot_condition(data, layout, conditions, varargin)
%
% INPUT:
%       data         = fnirs data
%       layout       = layout on which the graphs are plotted
%       conditions   = vector with trlinfo of the conditions you would like
%       to plot (e.g. [1 3] to plot condition 1 and 3)
%       condition_names = cell array with the names of the specified conditions
% Additional options can be specified in key-value pairs and can be:
%       'data_x'       = model response of fnirs data (default = []);
%       'baseline'   = vector with timepoints that should be used as baseline (default = [-10 0])
%       'trials'     = when true, all individual trials will also be
%       plotted (default = false)
%       'xlim'       = [xmin xmax];
%       'ylim'       = [ymin ymax];
%       'topoplot'   = 'yes' or 'no' (default = 'no')
%
% OUTPUT
%       h           = figure handle of the created plots 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% get options
data_x=ft_getopt(varargin, 'data_x', []);
baseline=ft_getopt(varargin, 'baseline', [-10 0]);
trials=ft_getopt(varargin, 'trials', false);
xlim=ft_getopt(varargin, 'xlim', [-20 20]);
ylim=ft_getopt(varargin, 'ylim', [-2 3]);
topoplot=ft_getopt(varargin, 'topoplot', 'no');

%% Perform baseline correction
cfg                 = [];
cfg.baseline        = baseline; 
cfg.parameter = 'trial';
data_blc = ft_timelockbaseline(cfg, data);

f=1;
for con=1:length(conditions)
    % average the data for the selected condition
    trl       = find(data.trialinfo(:,1)==conditions(con)); % index of trials of the selected condition
    cfg=[];
    cfg.trials=trl';
    cfg.avgoverrpt='yes';
    data_avg=ft_selectdata(cfg, data_blc);
    
    % select channels that are both in data and layout
    data_labels=data_blc.label;
    for i=1:length(data_labels)
        tmp = strsplit(data_labels{i});
        data_labels{i}=tmp{1};
    end
    [selchan, sellay] = match_str(data_labels, layout.label);% Take the subselection of channels that is contained in the layout, this is the same in all datasets
    
    % plot the data for the selected condition
    chanWidth=layout.width(sellay)*2;
    chanHeight=layout.height(sellay)*2;
    chanX=layout.pos(sellay,1);
    chanY=layout.pos(sellay,2);
    chanLabel  = layout.label(sellay);
    xmin=xlim(1);
    xmax=xlim(2);
    ymin=ylim(1);
    ymax=ylim(2);
    h{f}=figure; f=f+1;
    for c=1:2:length(selchan)
        xval=data_blc.time;
        for i=1:length(trl)
            tmpdat=data_blc.trial(trl(i),selchan([c c+1]), :);
            yval=permute(tmpdat, [2 3 1]);
            if trials==true
                hold on;
                ft_plot_vector(xval, yval, 'width', chanWidth(c), 'height', chanHeight(c), 'hpos', chanX(c), 'vpos', chanY(c), 'hlim', [xmin xmax], 'vlim', [ymin ymax], 'color', [0.5 0 0; 0 0 0.5],'linewidth', 0.05, 'axis', 'yes');
            end
        end
        yval=data_avg.trial(selchan([c c+1]), :);
        hold on;
        ft_plot_vector(xval, yval, 'width', chanWidth(c), 'height', chanHeight(c), 'hpos', chanX(c), 'vpos', chanY(c), 'hlim', [xmin xmax], 'vlim', [ymin ymax], 'color', 'rb','linewidth', 1, 'axis', 'yes');
    end
    if ~isempty(data_x)
      yval=data_x.trial{1};
      c=find(strcmp('model', layout.label));
      hold on;
      ft_plot_vector(xval, yval, 'width', layout.width(c), 'height', layout.height(c), 'hpos', layout.pos(c,1), 'vpos', layout.pos(c,2), 'hlim', [xmin xmax], 'vlim', [ymin ymax], 'color', 'r','linewidth', 0.1, 'axis', 'yes');
      ft_plot_layout(layout, 'box', 'yes', 'label', 'yes', 'outline', 'yes', 'point', 'no', 'mask', 'no', 'labelyoffset', 1.3*(layout.height(1)/2), 'labelalignh', 'center', 'chanindx', find(~ismember(layout.label, {'COMNT', 'SCALE'})) );
    end
    title(char(condition_names{con}))
    
    % plot the layout, labels and outline
    ft_plot_layout(layout, 'box', 'no', 'label', 'yes', 'outline', 'yes', 'point', 'no', 'mask', 'yes', 'fontsize', 8, 'labelyoffset', 1.4*median(layout.height/2), 'labelalignh', 'center', 'chanindx', find(~ismember(layout.label, {'COMNT', 'SCALE'})), 'interpreter', 'none');

    
    % Plot scale    
    l = find(strcmp(layout.label, 'SCALE'));
    x = layout.pos(l,1);
    y = layout.pos(l,2);
    plotScales([xmin xmax], [ymin ymax], x, y, chanWidth(1), chanHeight(1))
    
    
    
    %add the cfg/data/channel information to the figure under identifier linked to this axis
    ident                 = ['axh' num2str(round(sum(clock.*1e6)))]; % unique identifier for this axis
    set(gca, 'tag', ident);
    info                  = guidata(gcf);
    info.(ident).x        = layout.pos(:, 1);
    info.(ident).y        = layout.pos(:, 2);
    info.(ident).label    = layout.label;
    info.(ident).dataname = 'data_blc';
    info.(ident).cfg=layout.cfg;
    info.(ident).condition=conditions(con);
    info.(ident).varargin = data_blc;
    guidata(gcf, info);

%     info.(ident).cfg.layout  = layout_prep;
%     info.(ident).cfg.showlabels  = 'yes';
%     info.(ident).cfg.interactive = 'yes';
%     info.(ident).varargin = data_blc;
    
    set(gcf, 'WindowButtonUpFcn',  {@ft_select_channel, 'multiple', true, 'callback', {@select_singleplotER}, 'event', 'WindowButtonUpFcn'});
    set(gcf, 'WindowButtonDownFcn', {@ft_select_channel, 'multiple', true, 'callback', {@select_singleplotER}, 'event', 'WindowButtonDownFcn'});
    set(gcf, 'WindowButtonMotionFcn', {@ft_select_channel, 'multiple', true, 'callback', {@select_singleplotER}, 'event', 'WindowButtonMotionFcn'});

    if strcmp(topoplot, 'yes')
        % select O2Hb channels and remove [O2Hb] label
        cfg=[];
        cfg.channel='*[O2Hb]';
        data_O2Hb=ft_selectdata(cfg, data_avg);
        data_O2Hb.label=data_labels(contains(data_avg.label, '[O2Hb]'));
        % create topoplot
        cfg          = [];
        cfg.layout   = layout;
        cfg.marker   = 'labels';
        cfg.xlim     = [5 10]; % Choose the time window over which you want to average
        cfg.zlim     = ylim/2; % Choose the window over which you want to scale based on the plots
        h{f}=figure; f=f+1;
        ft_topoplotER(cfg, data_O2Hb);
        title(sprintf('Topoplot average O2Hb during %s [%d to %d sec]', condition_names{con}, cfg.xlim(1), cfg.xlim(2)));
        
        % select HHb channels and remove [HHb] label
        cfg=[];
        cfg.channel='*[HHb]';
        data_HHb=ft_selectdata(cfg, data_avg);
        data_HHb.label=data_labels(contains(data_avg.label, '[HHb]'));
        % create topoplot
        cfg          = [];
        cfg.layout   = layout;
        cfg.marker   = 'labels';
        cfg.xlim     = [5 10]; % Choose the time window over which you want to average
        cfg.zlim     = ylim/2; % Choose the window over which you want to scale based on the plots
        h{f}=figure; f=f+1;
        ft_topoplotER(cfg, data_HHb);
        title(sprintf('Topoplot average HHb during %s [%d to %d sec]', condition_names{con}, cfg.xlim(1), cfg.xlim(2)));
    
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotScales(hlim, vlim, hpos, vpos, width, height)

% the placement of all elements is identical
placement = {'hpos', hpos, 'vpos', vpos, 'width', width, 'height', height, 'hlim', hlim, 'vlim', vlim};

ft_plot_box([hlim vlim], placement{:}, 'edgecolor', 'k');

if hlim(1)<=0 && hlim(2)>=0
  ft_plot_line([0 0], vlim, placement{:}, 'color', 'k');
  ft_plot_text(0, vlim(1), '0  ', placement{:}, 'rotation', 90, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 8);
end

if vlim(1)<=0 && vlim(2)>=0
  ft_plot_line(hlim, [0 0], placement{:}, 'color', 'k');
  ft_plot_text(hlim(1), 0, '0  ', placement{:}, 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'middle', 'FontSize', 8);
end

ft_plot_text(hlim(1), vlim(1), [num2str(hlim(1), 3) ' '], placement{:}, 'rotation', 90, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top',    'FontSize', 8);
ft_plot_text(hlim(2), vlim(1), [num2str(hlim(2), 3) ' '], placement{:}, 'rotation', 90, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 8);
ft_plot_text(hlim(1), vlim(1), [num2str(vlim(1), 3) ' '], placement{:}, 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'bottom', 'FontSize', 8);
ft_plot_text(hlim(1), vlim(2), [num2str(vlim(2), 3) ' '], placement{:}, 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'top',    'FontSize', 8);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which is called after selecting channels in case of cfg.interactive='yes'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_singleplotER(label, data_blc)
% fetch cfg/data based on axis indentifier given as tag
ident       = get(gca,'tag');
info        = guidata(gcf);
cfg         = info.(ident).cfg;
datvarargin = info.(ident).varargin;
if ~isempty(label)
  cfg = removefields(cfg, 'inputfile');   % the reading has already been done and varargin contains the data
  cfg.baseline = 'no';                    % make sure the next function does not apply a baseline correction again
  cfg.channel = label;
  cfg.dataname = info.(ident).dataname;   % put data name in here, this cannot be resolved by other means
  cfg.trials = 'all';
  cfg.condition=info.(ident).condition;
  % trial selection has already been taken care of
%   fprintf('selected cfg.channel = {%s}\n', join_str(', ', cfg.channel));
  % ensure that the new figure appears at the same position
%   f = figure('position', get(gcf, 'Position'));
  singleplot_condition(cfg, datvarargin);
end
