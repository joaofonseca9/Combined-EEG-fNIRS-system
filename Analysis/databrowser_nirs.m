function databrowser_nirs(data, varargin)
bad_chan = ft_getopt(varargin, 'bad_chan', {}); % bad channels are indicated in black
cfg = [];
cfg.preproc.demean = 'yes'; % substracts the mean value (only in the plot)
cfg.viewmode = 'vertical';
cfg.channel = {'Rx*'};
cfg.fontsize=5;
cfg.blocksize  = 30;
cfg.nirsscale =150;
cfg.ylim = [-1 1];
cfg.linecolor = 'brkk';
cfg.colorgroups = repmat([1 2],1, length(data.label)/2)+2*ismember(data.label, bad_chan)';
ft_databrowser(cfg, data);
end
