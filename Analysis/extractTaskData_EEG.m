%% FUNCTION
% The usable data from the EEG is extracted into EEG(output), including
% baseline data (10 seconds before every start block)
% File will also get the different epochs: CHECK, autodual, autosingle,
% nonautodual, nonautosingle (cued and uncued)

function [EEG_divided, file]=extractTaskData_EEG(EEG, marker_table,results, file,mainpath_out)
%% Setting
pre=10;
post=10;

%% extract events of tasks and check/add correct number of events
event_samp  = [EEG.event.latency];

% n_cued_trials = sum([results.events_autodual.trial.cue])+...
%     sum([results.events_nonautodual.trial.cue])+...
%     sum([results.events_autosingle.trial.cue])+...
%     sum([results.events_nonautosingle.trial.cue]);

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

%% Remove data before the last test marker
% TEST_end = event_samp((strcmp({EEG.event.type}, 's1600'))==1);
% TEST_end = TEST_end(end);
% EEG=eeg_eegrej(EEG, [1 TEST_end+1]);

%% CHECKERBOARD_________________________________________________________
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
starts=CHECK_start-10*EEG.srate;
stops=CHECK_stop+10*EEG.srate;
[EEG_CHECK] = reject(starts, stops, EEG);
%save(file.CHECK, 'EEG_CHECK');

%% CUED TRIALS - the cued trials are easily separable (including the
% % baseline which is from Start_Cue -> Start_Block
% CUED_start=event_samp((strcmp({EEG.event.type}, sprintf('s%d',marker_table.StartCue(1))))==1);
% CUED_stop=event_samp((strcmp({EEG.event.type}, sprintf('s%d',marker_table.StopCue(1))))==1);
% 
% %Every cued trial (and 3 seconds before and after)
% starts=sort(CUED_start-10*EEG.srate);
% stops=sort(CUED_stop+10*EEG.srate);
% 
% EEG_cued=reject(starts,stops,EEG);
% %save(file.EEG_cued,'EEG_cued');

% AUTO DUAL CUE
Start_AutoDualCue  = event_samp((strcmp({EEG.event.type}, sprintf('s%d',marker_table.StartAutoDualCue(1))))==1);
Stop_AutoDualCue   = event_samp((strcmp({EEG.event.type}, sprintf('s%d',marker_table.StopAutoDualCue(1))))==1);

[EEG_AutoDualCue] = reject(Start_AutoDualCue-pre*EEG.srate, Stop_AutoDualCue+post*EEG.srate,EEG);
%save(file.EEG_AutoDualCue,'EEG_AutoDualCue');

CUED_start=Start_AutoDualCue;
CUED_stop=Stop_AutoDualCue;


% AUTO SINGLE CUE
Start_AutoCue  = event_samp((strcmp({EEG.event.type}, sprintf('s%d',marker_table.StartAutoCue(1))))==1);
Stop_AutoCue   = event_samp((strcmp({EEG.event.type}, sprintf('s%d',marker_table.StopAutoCue(1))))==1);
[EEG_AutoCue] = reject(Start_AutoCue-pre*EEG.srate, Stop_AutoCue+post*EEG.srate,EEG);
%save(file.EEG_AutoCue,'EEG_AutoCue');
CUED_start=[CUED_start Start_AutoCue];
CUED_stop=[CUED_stop Stop_AutoCue];

% NON AUTO DUAL CUE
Start_NonAutoDualCue  = event_samp((strcmp({EEG.event.type}, sprintf('s%d',marker_table.StartNonAutoDualCue(1))))==1);
Stop_NonAutoDualCue   = event_samp((strcmp({EEG.event.type}, sprintf('s%d',marker_table.StopNonAutoDualCue(1))))==1);
[EEG_NonAutoDualCue] = reject(Start_NonAutoDualCue-pre*EEG.srate, Stop_NonAutoDualCue+post*EEG.srate,EEG);
%save(file.EEG_NonAutoDualCue,'EEG_NonAutoDualCue');
CUED_start=[CUED_start Start_NonAutoDualCue];
CUED_stop=[CUED_stop Stop_NonAutoDualCue];

% NON AUTO SINGLE CUE
Start_NonAutoCue  = event_samp((strcmp({EEG.event.type}, sprintf('s%d',marker_table.StartNonAutoCue(1))))==1);
Stop_NonAutoCue   = event_samp((strcmp({EEG.event.type}, sprintf('s%d',marker_table.StopNonAutoCue(1))))==1);
[EEG_NonAutoCue] = reject(Start_NonAutoCue-pre*EEG.srate, Stop_NonAutoCue+post*EEG.srate,EEG);
%save(file.EEG_NonAutoCue,'EEG_NonAutoCue');
CUED_start=[CUED_start Start_NonAutoCue];
CUED_stop=[CUED_stop Stop_NonAutoCue];

EEG_cued=reject(CUED_start-pre*EEG.srate,CUED_stop+post*EEG.srate,EEG);

%% UNCUED TRIALS - the uncued trials have a cue that only lasts 8 seconds;
%to separate them we get the Task Blocks and then add the baseline from
%Start_Metronome -> Start_Block

%% AUTO DUAL NO CUE
Start_AutoDualNoCue  = event_samp((strcmp({EEG.event.type}, sprintf('s%d',marker_table.StartAutoDualNoCue(1))))==1);
Stop_AutoDualNoCue   = event_samp((strcmp({EEG.event.type}, sprintf('s%d',marker_table.StopAutoDualNoCue(1))))==1);
[EEG_AutoDualNoCue] = reject(Start_AutoDualNoCue-pre*EEG.srate, Stop_AutoDualNoCue+post*EEG.srate,EEG);
%save(file.EEG_AutoDualNoCue,'EEG_AutoDualNoCue');

% %Get the baseline
% Start_Metronome_AutoDual=strcmp({EEG.event.type}, sprintf('s%d',marker_table.StartAutoDualNoCue(1)))==1;
% %index of the metronome events in the event vector
% idx=find(Start_Metronome_AutoDual);
% UNCUED_start=event_samp(idx-2)-pre*EEG.srate; %the Start Metronome markers are 2 before the start block marker(Start Metronome-StopMetronome-StartBlock)
% UNCUED_stop=Stop_AutoDualNoCue+post*EEG.srate;

%% NON AUTO DUAL
Start_NonAutoDualNoCue  = event_samp((strcmp({EEG.event.type}, sprintf('s%d',marker_table.StartNonAutoDualNoCue(1))))==1);
Stop_NonAutoDualNoCue   = event_samp((strcmp({EEG.event.type}, sprintf('s%d',marker_table.StopNonAutoDualNoCue(1))))==1);
[EEG_NonAutoDualNoCue] = reject(Start_NonAutoDualNoCue-pre*EEG.srate, Stop_NonAutoDualNoCue+post*EEG.srate,EEG);
%save(file.EEG_NonAutoDualNoCue,'EEG_NonAutoDualNoCue');

% %Get baseline
% Start_Metronome_NonAutoDual=strcmp({EEG.event.type}, sprintf('s%d',marker_table.StartNonAutoDualNoCue(1)))==1;
% %index of the metronome events in the event vector
% idx=find(Start_Metronome_NonAutoDual);
UNCUED_start    = [Start_NonAutoDualNoCue];
UNCUED_stop     = [Stop_NonAutoDualNoCue];

%% AUTO SINGLE
Start_AutoNoCue  = event_samp((strcmp({EEG.event.type}, sprintf('s%d',marker_table.StartAutoNoCue(1))))==1);
Stop_AutoNoCue   = event_samp((strcmp({EEG.event.type}, sprintf('s%d',marker_table.StopAutoNoCue(1))))==1);
[EEG_AutoNoCue] = reject(Start_AutoNoCue-pre*EEG.srate, Stop_AutoNoCue+post*EEG.srate,EEG);
%save(file.EEG_AutoSingleNoCue,'EEG_AutoSingleNoCue');

% %Get baseline
% Start_Metronome_AutoSingle=strcmp({EEG.event.type}, sprintf('s%d',marker_table.StartAutoNoCue(1)))==1;
% %index of the metronome events in the event vector
% idx=find(Start_Metronome_AutoSingle);
UNCUED_start=[UNCUED_start Start_AutoNoCue];
UNCUED_stop     = [UNCUED_stop Stop_AutoNoCue];

%% NON AUTO SINGLE
Start_NonAutoNoCue  = event_samp((strcmp({EEG.event.type}, sprintf('s%d',marker_table.StartNonAutoNoCue(1))))==1);
Stop_NonAutoNoCue   = event_samp((strcmp({EEG.event.type}, sprintf('s%d',marker_table.StopNonAutoNoCue(1))))==1);
[EEG_NonAutoNoCue] = reject(Start_NonAutoNoCue-pre*EEG.srate, Stop_NonAutoNoCue+post*EEG.srate,EEG);
%save(file.EEG_NonAutoSingleNoCue,'EEG_NonAutoSingleNoCue');

% %Get baseline
% Start_Metronome_NonAutoSingle=strcmp({EEG.event.type}, sprintf('s%d',marker_table.StartNonAutoNoCue(1)))==1;
% %index of the metronome events in the event vector
% idx=find(Start_Metronome_NonAutoSingle);
UNCUED_start=[UNCUED_start Start_NonAutoNoCue];
UNCUED_stop     = [UNCUED_stop Stop_NonAutoNoCue];


[EEG_uncued]    =reject(UNCUED_start-pre*EEG.srate, UNCUED_stop+post*EEG.srate,EEG);

[EEG_task]=reject([CUED_start UNCUED_start]-pre*EEG.srate,...
        [CUED_stop UNCUED_stop]+post*EEG.srate,EEG);
%save(file.EEG_task,'EEG_task');


EEG_divided=table(EEG_task,EEG_cued,EEG_uncued,EEG_AutoDualCue,EEG_AutoCue,...
    EEG_NonAutoDualCue,EEG_NonAutoCue,EEG_AutoDualNoCue,EEG_AutoNoCue,...
    EEG_NonAutoDualNoCue,EEG_NonAutoNoCue);
end


function subset=reject(starts,stops,EEG)
event_samp  = [EEG.event.latency];
starts=sort(starts);
stops=sort(stops);

rej=[2 starts(1)];
for ii=1:length(stops)-1
    rej=[rej;stops(ii) starts(ii+1)];
    if ii==length(stops)-1
        rej=[rej;stops(ii+1) event_samp(end)];
    end
end

[subset] = eeg_eegrej(EEG, rej);

end