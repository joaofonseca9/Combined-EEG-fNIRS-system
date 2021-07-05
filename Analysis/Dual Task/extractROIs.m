function [nirs_HbO2_DLPFC, nirs_HbO2_SMA, nirs_HbO2_M1, nirs_HbO2_PPC] =...
    extractROIs(nirs_TLO2Hb)
% DLPFC: Rx5-Tx7, Rx5-Tx8, Rx7-Tx7, Rx7-Tx8, Rx9-Tx13, Rx9-Tx12, Rx11-Tx12,
% Rx11-Tx13
cfg = [];
cfg.channel = {'Rx5-Tx7', 'Rx5-Tx8', 'Rx7-Tx7', 'Rx7-Tx8', 'Rx9-Tx13',...
    'Rx9-Tx12', 'Rx11-Tx12', 'Rx11-Tx13'};
nirs_HbO2_DLPFC{1} = ft_selectdata(cfg, nirs_TLO2Hb{1});
nirs_HbO2_DLPFC{2} = ft_selectdata(cfg, nirs_TLO2Hb{2});
nirs_HbO2_DLPFC{3} = ft_selectdata(cfg, nirs_TLO2Hb{3});
nirs_HbO2_DLPFC{4} = ft_selectdata(cfg, nirs_TLO2Hb{4});

% PMC/SMA: Rx4-Tx5, Rx3-Tx5, Rx4-Tx4
cfg = [];
cfg.channel = {'Rx4-Tx5', 'Rx3-Tx5', 'Rx4-Tx4'};
nirs_HbO2_SMA{1} = ft_selectdata(cfg, nirs_TLO2Hb{1});
nirs_HbO2_SMA{2} = ft_selectdata(cfg, nirs_TLO2Hb{2});
nirs_HbO2_SMA{3} = ft_selectdata(cfg, nirs_TLO2Hb{3});
nirs_HbO2_SMA{4} = ft_selectdata(cfg, nirs_TLO2Hb{4});

% M1: Rx3-Tx2, Rx1-Tx2, Rx3-Tx3, Rx1-Tx3, Rx2-Tx4, Rx2-Tx3
cfg = [];
cfg.channel = {'Rx3-Tx2', 'Rx1-Tx2', 'Rx3-Tx3', 'Rx1-Tx3', 'Rx2-Tx4',...
    'Rx2-Tx3'};
nirs_HbO2_M1{1} = ft_selectdata(cfg, nirs_TLO2Hb{1});
nirs_HbO2_M1{2} = ft_selectdata(cfg, nirs_TLO2Hb{2});
nirs_HbO2_M1{3} = ft_selectdata(cfg, nirs_TLO2Hb{3});
nirs_HbO2_M1{4} = ft_selectdata(cfg, nirs_TLO2Hb{4});

% PPC: Rx8-Tx10, Rx6-Tx9, Rx8-Tx9, Rx12-Tx15, Rx10-Tx14, Rx12-Tx14
cfg = [];
cfg.channel = {'Rx8-Tx10', 'Rx6-Tx9', 'Rx8-Tx9', 'Rx12-Tx15',...
    'Rx10-Tx14', 'Rx12-Tx14'};
nirs_HbO2_PPC{1} = ft_selectdata(cfg, nirs_TLO2Hb{1});
nirs_HbO2_PPC{2} = ft_selectdata(cfg, nirs_TLO2Hb{2});
nirs_HbO2_PPC{3} = ft_selectdata(cfg, nirs_TLO2Hb{3});
nirs_HbO2_PPC{4} = ft_selectdata(cfg, nirs_TLO2Hb{4});

end