function singleplot_condition(tmpcfg, data)

% Define the channel
data_labels=data.label;
for i=1:length(data_labels)
    tmp = strsplit(data_labels{i});
    data_labels{i}=tmp{1};
end
c = match_str(data_labels, tmpcfg.channel);

% average the data for the selected condition
trl       = find(data.trialinfo(:,1)==tmpcfg.condition); % index of trials of the selected condition
cfg=[];
cfg.trials=trl';
cfg.avgoverrpt='yes';
data_avg=ft_selectdata(cfg, data);

% Plot the data
xmin=-25;
xmax=30;
ymin=-2;
ymax=3;
figure;
for i=1:length(trl)
    xval=data.time;
    tmpdat=data.trial(trl(i), c, :);
    ytemp=permute(tmpdat, [2 3 1]);
    yval(1,:)=mean(ytemp(1:2:length(c)-1,:),1);
    yval(2,:)=mean(ytemp(2:2:length(c),:),1);
    hold on;
    ft_plot_vector(xval, yval, 'color', [0.5 0 0; 0 0 0.5], 'linewidth', 0.1, 'hlim', [xmin xmax], 'vlim', [ymin ymax], 'axis', 'yes');
end
ytemp=data_avg.trial(c,:);
yval(1,:)=mean(ytemp(1:2:length(c)-1,:),1);
yval(2,:)=mean(ytemp(2:2:length(c),:),1);
hold on;
ft_plot_vector(xval, yval, 'color', 'rb', 'linewidth', 2, 'hlim', [xmin xmax], 'vlim', [ymin ymax], 'axis', 'yes');
title(char(tmpcfg.channel))

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

    