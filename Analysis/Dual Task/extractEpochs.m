function [data, time] = extractEpochs(EEG)
% extract epochs (trials) using the event-markers [channels x time x trial]
% ERP-trial: -0.1s to 0.4s;     
% create time vector in seconds

    % task-events
    event_samp  = [EEG.event.latency];
    event_assr  = event_samp((strcmp({EEG.event.type}, 's1777') | (strcmp({EEG.event.type}, 's1797')))==1);
    
    % extract trials ERP
    for iTrial = 1:length(event_assr)
        data(:,:,iTrial) = EEG.data(:, [event_assr(iTrial)-floor(0.1*EEG.srate) : event_assr(iTrial)+ceil(0.7*EEG.srate)]);
    end
    % create time vector [s]
    time = -0.1 : 1/EEG.srate : 0.7;
end