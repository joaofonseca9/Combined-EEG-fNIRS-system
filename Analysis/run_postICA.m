function [EEG] = run_postICA (EEG)
% use the Independent Components, defined using the ICA
% check independent components and look for eyeblink IC, using topoplots,
% power spectra and raw signal, then eliminate IC of eyeblinks 

    % plot raw data, to observe eyeblinks
    figure;
    subplot(2,2,1); plot(EEG.data(1,1000:10000)); title('preICA - Fp1')
    subplot(2,2,2); plot(EEG.data(2,1000:10000)); title('preICA - Fpz')
    % plot all IC headplots and timeseries to select eye component
    pop_selectcomps(EEG);                   
    pop_eegplot(EEG, 0,1,1);                
    % plot power spectrum of separate IC, defined in command window
    idx = 1;
    while idx == 1
        IC = input('check for IC eyeblink [nr/no]: ','s');  % select CHAN of interest
        if ~strcmp(IC, 'no')
            pop_prop(EEG,0,str2double(IC),NaN,{'freqrange' [2 55]}); % plot eyeblink IC
        else
            idx = 0;
        end
    end
    % remove eyeblink IC from original data, defined in command window
    idx = 1;
    while idx == 1
        remove = input('remove eyeblink component [nr/no]: ','s');
        if ~strcmp(remove, 'no')
            [EEGcheck] = pop_subcomp(EEG,str2double(remove),1);
            figure(10);
            subplot(2,2,3); plot(EEGcheck.data(1,1000:10000)); title('postICA - Fp1')
            subplot(2,2,4); plot(EEGcheck.data(2,1000:10000)); title('postICA - Fpz')
            check = input('sure to remove eyeblink component? [yes/no]: ','s');
            switch check
                case 'yes';     EEG = EEGcheck;
                case 'no';      disp('select other eyeblink IC\n')
            end
        else
            idx = 0;
        end
    end
end