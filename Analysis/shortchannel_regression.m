function data_rcs=shortchannel_regression(cfg, datain);

% shortchannelregression performs reference channel subtraction for NIRS
% data. It is an adaptation from ft_nirs_referencechannelsubstraction 
% (http://www.fieldtriptoolbox.org/reference/ft_nirs_referencechannelsubtraction/)
%
% Use as
%   outdata = shortchannel_regression(cfg, indata)
% where indata is nirs data and cfg is a configuration structure that should contain
%
%   cfg.method        = string, 'regstat2', 'QR' or 'OLS' (default = 'QR')
%   cfg.verbose       = boolean, whether text output is desired (default = false)
%

%%
% get the options
cfg.method        = ft_getopt(cfg, 'method', 'QR');
cfg.verbose       = ft_getopt(cfg, 'verbose', false);

% find short and long channel indexes
SC=find(contains(datain.label, {'a ', 'b ', 'c ', 'd '})); %index of all short channels
LC=find(~contains(datain.label, {'a ', 'b ', 'c ', 'd '})); %index of all long channels


data_rcs			 = datain;
data_rcs.label = datain.label(LC);
for tr=1:numel(datain.trial)
  %data = datain.trial{tr};
  shallow		= datain.trial{tr}(SC,:);
  shallow		= bsxfun(@minus,shallow,mean(shallow,2)); % mean detrend
  shallowlabel = datain.label(SC);
  
  deep		= datain.trial{tr}(LC,:);
  deep		= bsxfun(@minus,deep,mean(deep,2)); % mean detrend
  deeplabel = datain.label(LC);
  
  time		= datain.time{tr};
  
  % Reference channel subtraction
  ndeep		= size(deep,1);
  signal		= NaN(size(deep));
  x			= shallow';
  for dpIdx	= 1:ndeep
    y				= deep(dpIdx,:)';
    switch (cfg.method)
      case 'regstat2'
        b				= regstats2(y,x,'linear',{'beta','r'});
        beta    = b.beta;
        res			= b.r;
        
      case 'QR'
        [Q,R] = qr(x,0);
        beta = R\(Q'*y);
        yhat = x*beta;
        res = y - yhat;
        
      case 'OLS'
        x2 = [repmat(1, size(x, 1), 1) x];
        beta = x2\y;
        yhat = x2*beta;
        res  = y - yhat;
        cfg.verbose;
        beta(1) = [];
        
      otherwise % it should never come here as we use ft_checkopt
        error('unrecognized method');
    end
    
    signal(dpIdx,:) = res';
    
    % sanity check of results
    if cfg.verbose
      if fprintf('Found the following meaningful shallow channels for deep channel %s:', deeplabel{dpIdx});
        shIdx = find(beta>0.5);
        for s=1:numel(shIdx)
          fprintf('\n\t%s', shallowlabel{shIdx(s)})
        end
      end
      fprintf('\n\n')
    end
  end
  
  % overwrite
  data_rcs.time{tr}	= time;
  data_rcs.trial{tr}	= signal;
end