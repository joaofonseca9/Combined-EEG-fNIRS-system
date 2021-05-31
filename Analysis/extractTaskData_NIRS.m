function [NIRS_divided, file]=extractTaskData_NIRS(data_raw, data_down, event, marker_table,results, file,mainpath_out)
%% 
    TEST_end=find(strcmp({event.value}, 'LSL 1600'));
    TEST_end=TEST_end(end);
    smp.TEST_end=[eeg_fnirs_events(TEST_end).sample]';
    trl.TEST_end=[TEST_end size(data_raw.trial{1},2)];
    cfg     = [];
    cfg.trl = trl;
    data_notest = ft_redefinetrial(cfg,data_raw);
    save('nirs/data_notest.mat', 'data_notest');
end