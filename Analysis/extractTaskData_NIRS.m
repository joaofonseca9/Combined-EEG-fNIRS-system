function [NIRS_divided, file]=extractTaskData_NIRS(data_raw, data_down, event, marker_table,results, file,mainpath_out)
%% 
%% Define events + epochs

%% Starts
%Cued
start_autodual_cue = find(strcmp({event.value}, sprintf('LSL %d',marker_table.StartAutoDualCue(1))));
start_autosingle_cue =  find(strcmp({event.value}, sprintf('LSL %d',marker_table.StartAutoCue(1))));
start_nonautodual_cue= find(strcmp({event.value}, sprintf('LSL %d',marker_table.StartNonAutoDualCue(1))));
start_nonautosingle_cue= find(strcmp({event.value}, sprintf('LSL %d',marker_table.StartNonAutoCue(1))));

%Uncued
start_autodual_nocue = find(strcmp({event.value}, sprintf('LSL %d',marker_table.StartAutoDualNoCue(1))));
start_autosingle_nocue =  find(strcmp({event.value}, sprintf('LSL %d',marker_table.StartAutoNoCue(1))));
start_nonautodual_nocue= find(strcmp({event.value}, sprintf('LSL %d',marker_table.StartNonAutoDualNoCue(1))));
start_nonautosingle_nocue= find(strcmp({event.value}, sprintf('LSL %d',marker_table.StartNonAutoNoCue(1))));


%% Stops
%Cued
stop_autodual_cue = find(strcmp({event.value}, sprintf('LSL %d',marker_table.StopAutoDualCue(1))));
stop_autosingle_cue =  find(strcmp({event.value}, sprintf('LSL %d',marker_table.StopAutoCue(1))));
stop_nonautodual_cue= find(strcmp({event.value}, sprintf('LSL %d',marker_table.StopNonAutoDualCue(1))));
stop_nonautosingle_cue= find(strcmp({event.value}, sprintf('LSL %d',marker_table.StopNonAutoCue(1))));

%Uncued
stop_autodual_nocue = find(strcmp({event.value}, sprintf('LSL %d',marker_table.StopAutoDualNoCue(1))));
stop_autosingle_nocue =  find(strcmp({event.value}, sprintf('LSL %d',marker_table.StopAutoNoCue(1))));
stop_nonautodual_nocue= find(strcmp({event.value}, sprintf('LSL %d',marker_table.StopNonAutoDualNoCue(1))));
stop_nonautosingle_nocue= find(strcmp({event.value}, sprintf('LSL %d',marker_table.StopNonAutoNoCue(1))));

% get the sample number in the original dataset
% note that we transpose them to get columns
smp.start_autodual_cue = [event(start_autodual_cue).sample]';
smp.start_autosingle_cue = [event(start_autosingle_cue).sample]';
smp.start_nonautodual_cue = [event(start_nonautodual_cue).sample]';
smp.start_nonautosingle_cue = [event(start_nonautosingle_cue).sample]';
smp.start_autodual_nocue = [event(start_autodual_nocue).sample]';
smp.start_autosingle_nocue = [event(start_autosingle_nocue).sample]';
smp.start_nonautodual_nocue = [event(start_nonautodual_nocue).sample]';
smp.start_nonautosingle_nocue = [event(start_nonautosingle_nocue).sample]';


smp.stop_autodual_cue = [event(stop_autodual_cue).sample]';
smp.stop_autosingle_cue = [event(stop_autosingle_cue).sample]';
smp.stop_nonautodual_cue = [event(stop_nonautodual_cue).sample]';
smp.stop_nonautosingle_cue = [event(stop_nonautosingle_cue).sample]';
smp.stop_autodual_nocue = [event(stop_autodual_nocue).sample]';
smp.stop_autosingle_nocue = [event(stop_autosingle_nocue).sample]';
smp.stop_nonautodual_nocue = [event(stop_nonautodual_nocue).sample]';
smp.stop_nonautosingle_nocue = [event(stop_nonautosingle_nocue).sample]';

%% Get the sample number after downsampling
factor = data_raw.fsample / data_down.fsample;


% ! maybe convert this structure-array into a matrix?
smp.start_autodual_cue = round((smp.start_autodual_cue-1)/factor+1);
smp.start_autosingle_cue= round((smp.start_autosingle_cue-1)/factor+1);
smp.start_nonautodual_cue = round((smp.start_nonautodual_cue-1)/factor+1);
smp.start_nonautosingle_cue= round((smp.start_nonautosingle_cue-1)/factor+1);
smp.start_autodual_nocue = round((smp.start_autodual_nocue-1)/factor+1);
smp.start_autosingle_nocue= round((smp.start_autosingle_nocue-1)/factor+1);
smp.start_nonautodual_nocue = round((smp.start_nonautodual_nocue-1)/factor+1);
smp.start_nonautosingle_nocue= round((smp.start_nonautosingle_nocue-1)/factor+1);



smp.stop_autodual_cue = round((smp.stop_autodual_cue-1)/factor+1);
smp.stop_autosingle_cue= round((smp.stop_autosingle_cue-1)/factor+1);
smp.stop_nonautodual_cue = round((smp.stop_nonautodual_cue-1)/factor+1);
smp.stop_nonautosingle_cue= round((smp.stop_nonautosingle_cue-1)/factor+1);
smp.stop_autodual_nocue = round((smp.stop_autodual_nocue-1)/factor+1);
smp.stop_autosingle_nocue= round((smp.stop_autosingle_nocue-1)/factor+1);
smp.stop_nonautodual_nocue = round((smp.stop_nonautodual_nocue-1)/factor+1);
smp.stop_nonautosingle_nocue= round((smp.stop_nonautosingle_nocue-1)/factor+1);

%% Get offset and pre and post
pre    =  round(2*data_down.fsample); % seconds pre-stimulus - baseline (20s-25s)
post   =  round(2*data_down.fsample); % seconds post-stimulus 
offset = -pre; % see ft_definetrial

%% Setting the trials
trl.autodual_cue = [smp.start_autodual_cue-pre smp.stop_autodual_cue+post];
trl.autosingle_cue = [smp.start_autosingle_cue-pre smp.stop_autosingle_cue+post];
trl.nonautodual_cue = [smp.start_nonautodual_cue-pre smp.stop_nonautodual_cue+post];
trl.nonautosingle_cue = [smp.start_nonautosingle_cue-pre smp.stop_nonautosingle_cue+post];

trl.autodual_nocue = [smp.start_autodual_nocue-pre smp.stop_autodual_nocue+post];
trl.autosingle_nocue = [smp.start_autosingle_nocue-pre smp.stop_autosingle_nocue+post];
trl.nonautodual_nocue = [smp.start_nonautodual_nocue-pre smp.stop_nonautodual_nocue+post];
trl.nonautosingle_nocue = [smp.start_nonautosingle_nocue-pre smp.stop_nonautosingle_nocue+post];


%% Add the offset
trl.autodual_cue(:,3) = offset;
trl.autosingle_cue(:,3) = offset;
trl.nonautodual_cue(:,3) = offset;
trl.nonautosingle_cue(:,3) = offset;

trl.autodual_nocue(:,3) = offset;
trl.autosingle_nocue(:,3) = offset;
trl.nonautodual_nocue(:,3) = offset;
trl.nonautosingle_nocue(:,3) = offset;

%% trialinfo
trl.autodual_cue(:,4) = 1;
trl.autosingle_cue(:,4) = 2;
trl.nonautodual_cue(:,4) = 3;
trl.nonautosingle_cue(:,4) = 4;

trl.autodual_nocue(:,4) = 5;
trl.autosingle_nocue(:,4) = 6;
trl.nonautodual_nocue(:,4) = 7;
trl.nonautosingle_nocue(:,4) = 8;

end