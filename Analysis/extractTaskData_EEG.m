function [EEG, starts, stops]=extractTaskData_EEG(EEG, marker_table,results, file)
% extract events of tasks and check/add correct number of events
event_samp  = [EEG.event.latency];

n_cued_trials = sum([results.events_autodual.trial.cue])+...
    sum([results.events_nonautodual.trial.cue])+...
    sum([results.events_autosingle.trial.cue])+...
    sum([results.events_nonautosingle.trial.cue]);

N_trials = length([results.events_autodual.trial.cue])+...
    length([results.events_nonautodual.trial.cue])+...
    length([results.events_autosingle.trial.cue])+...
    length([results.events_nonautosingle.trial.cue]);

%% CHECK COUNT - only check key, letters and flips
n_letters   = 8*N_trials/2;
n_movie     = N_trials/2;
n_key       = N_trials*12;
n_flips     = 300;
%Letters
if marker_table.Letter(2)==n_letters
    disp('Correct number of letters for EEG only');
else
    disp('Incorrect number of letters for EEG only');
end
%Keys
if marker_table.Key(2)==n_key
    disp('Correct number of keys for EEG only');
else
    disp('Incorrect number of keys for EEG only');
end
%Flips
if marker_table.CheckFlip(2)==n_flips
    disp('Correct number of checkboard flips for EEG only');
else
    disp('Incorrect number of checkboard flips for EEG only');
end

% CHECKERBOARD_________________________________________________________
CHECK_start = event_samp((strcmp({EEG.event.type}, 's1555'))==1);
CHECK_stop  = event_samp((strcmp({EEG.event.type}, 's1500'))==1);
CHECK_flip  = event_samp((strcmp({EEG.event.type}, 's1255'))==1);

if size(CHECK_start,2) == 1 && size(CHECK_stop,2) == 1 && size(CHECK_flip,2) == 300 
    disp('CHECK : correct number of events'); CHECK = true;
elseif size(CHECK_start,2) < 1 || isempty(CHECK_start)  % start-event is missing
    CHECK_start = CHECK_flip(1)-10*EEG.srate;
    disp('CHECK : correct number of events (corrected)'); CHECK = true;
else
    error('CHECK : INCORRECT NUMBER OF EVENTS');
end

% CUED TRIALS - the cued trials are easily separable (including the
% baseline which is from Start_Cue -> Start_Block
CUED_start=event_samp((strcmp({EEG.event.type}, sprintf('s%d',marker_table.StartCue(1))))==1);
CUED_stop=event_samp((strcmp({EEG.event.type}, sprintf('s%d',marker_table.StopCue(1))))==1);

% UNCUED TRIALS - the uncued trials have a cue that only lasts 8 seconds;
%to separate them we get the Task Blocks and then add the baseline from
%Start_Metronome -> Start_Block

% AUTO DUAL
Start_AutoDual  = event_samp((strcmp({EEG.event.type}, sprintf('s%d',marker_table.StartAutoDualNoCue(1))))==1);
Stop_AutoDual   = event_samp((strcmp({EEG.event.type}, sprintf('s%d',marker_table.StopAutoDualNoCue(1))))==1);
Start_Metronome_AutoDual=strcmp({EEG.event.type}, sprintf('s%d',marker_table.StartAutoDualNoCue(1)))==1;
%index of the metronome events in the event vector
idx=find(Start_Metronome_AutoDual);
UNCUED_start=event_samp(idx-2)-5*EEG.srate; %the Start Metronome markers are 2 before the start block marker(Start Metronome-StopMetronome-StartBlock)
UNCUED_stop=Stop_AutoDual+5*EEG.srate;

% NON AUTO DUAL
Start_NonAutoDual  = event_samp((strcmp({EEG.event.type}, sprintf('s%d',marker_table.StartNonAutoDualNoCue(1))))==1);
Stop_NonAutoDual   = event_samp((strcmp({EEG.event.type}, sprintf('s%d',marker_table.StopNonAutoDualNoCue(1))))==1);
Start_Metronome_NonAutoDual=strcmp({EEG.event.type}, sprintf('s%d',marker_table.StartNonAutoDualNoCue(1)))==1;
%index of the metronome events in the event vector
idx=find(Start_Metronome_NonAutoDual);
UNCUED_start    = [UNCUED_start event_samp(idx-2)-5*EEG.srate];
UNCUED_stop     = [UNCUED_stop Stop_NonAutoDual+5*EEG.srate];

% AUTO SINGLE
Start_AutoSingle  = event_samp((strcmp({EEG.event.type}, sprintf('s%d',marker_table.StartAutoNoCue(1))))==1);
Stop_AutoSingle   = event_samp((strcmp({EEG.event.type}, sprintf('s%d',marker_table.StopAutoNoCue(1))))==1);
Start_Metronome_AutoSingle=strcmp({EEG.event.type}, sprintf('s%d',marker_table.StartAutoNoCue(1)))==1;
%index of the metronome events in the event vector
idx=find(Start_Metronome_AutoSingle);
UNCUED_start=[UNCUED_start event_samp(idx-2)-5*EEG.srate];
UNCUED_stop     = [UNCUED_stop Stop_AutoSingle+5*EEG.srate];

% NON AUTO SINGLE
Start_NonAutoSingle  = event_samp((strcmp({EEG.event.type}, sprintf('s%d',marker_table.StartNonAutoNoCue(1))))==1);
Stop_NonAutoSingle   = event_samp((strcmp({EEG.event.type}, sprintf('s%d',marker_table.StopNonAutoNoCue(1))))==1);
Start_Metronome_NonAutoSingle=strcmp({EEG.event.type}, sprintf('s%d',marker_table.StartNonAutoNoCue(1)))==1;
%index of the metronome events in the event vector
idx=find(Start_Metronome_NonAutoSingle);
UNCUED_start=[UNCUED_start event_samp(idx-2)-5*EEG.srate];
UNCUED_stop     = [UNCUED_stop Stop_NonAutoSingle+5*EEG.srate];

starts=sort([CUED_start-5*EEG.srate UNCUED_start]);
stops=sort([CUED_stop+5*EEG.srate UNCUED_stop]);

rej=[2 starts(1)];
for ii=1:length(stops)-1
    rej=[rej;stops(ii) starts(ii+1)];
    if ii==length(stops)-1
        rej=[rej;stops(ii+1) event_samp(end)];
    end
end

[EEG] = eeg_eegrej(EEG, rej);



end