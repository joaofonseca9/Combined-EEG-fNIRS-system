function [nirs_TLO2Hb, nirs_TLHHb] = separateHbO2FromHb(conditions, nirs_TLblc, subject)
% Separate O2Hb and HHb channels
for con = 1:length(conditions)
    cfg = [];
    cfg.channel = '* [O2Hb]';
    nirs_TLO2Hb{con} = ft_selectdata(cfg, nirs_TLblc{con}{subject});
    
    % Rename labels such that they have the same name as HHb channels
    for i = 1:length(nirs_TLO2Hb{con}.label)
        tmp = strsplit(nirs_TLO2Hb{con}.label{i});
        nirs_TLO2Hb{con}.label{i}=tmp{1};
    end
    
    % The same for HHb channels
    cfg = [];
    cfg.channel = '* [HHb]';
    nirs_TLHHb{con} = ft_preprocessing(cfg, nirs_TLblc{con}{subject});
    for i=1:length(nirs_TLHHb{con}.label)
        tmp = strsplit(nirs_TLHHb{con}.label{i});
        nirs_TLHHb{con}.label{i}=tmp{1};
    end
end
end