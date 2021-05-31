function [NIRS_divided, file]=extractTaskData_NIRS(data_raw, data_down, event, marker_table,results, file,mainpath_out)
%% 
    TEST_end=TEST_end(end);
    smp.TEST_end=[eeg_fnirs_events(TEST_end).sample]';
    trl=[TEST_end size(data_raw.trial{1},2)];
    pre    =  round(20*data_raw.fsample); % seconds pre-stimulus (100=10 seconden)
    post   =  round(20*data_raw.fsample); % seconds post-stimulus 
    offset = -pre; % see ft_definetrial
    trl(:,3) = offset;
    cfg     = [];
    cfg.inputfile = 'sub-03_rec-02_nirs.mat';
    cfg.trl = trl;
    data_notest = ft_redefinetrial(cfg);
    save('data_notest.mat', 'data_notest');
end