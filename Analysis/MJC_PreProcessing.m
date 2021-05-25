clear;

%% Initialize FieldTrip & EEGLAB
% laptop='laptopJoao';
laptop='laptopMariana';
% laptop='laptopCatarina';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);
eeglab;
ft_defaults;

sub='02';
rec='02';

file = getFileNames(mainpath_out, sub, rec);

%% Select Data Folder (if pilots, select pilots folder)
sub_path    = fullfile(mainpath_in,'incoming',['sub-',sub]);
eeg_path = fullfile(sub_path,'eeg');
nirseeg_path = fullfile(sub_path,'nirseeg');
sub_vhdr    = fullfile(['sub-',sub,'_rec-',rec,'_eeg.vhdr']);
cd(sub_path);

oxyfile=fullfile(sub_path,'nirseeg',['sub-',sub,'_rec-',rec,'_nirseeg.edf']);

%% Upload files
% fprintf('\n \n Select the Oxysoft file \n \n')
% [baseName, folder] = uigetfile({'*.oxy4;*.oxy3;*.edf', 'All Oxysoft file(*.oxy4, *.oxy3, *.edf)';},'Select the Oxysoft file');
% oxyfile = fullfile(folder, baseName);
% oxyfile = 'C:\Users\maria\Universidade do Porto\João Pedro Barbosa Fonseca - Internship\Experiment\Data\Pilots\incoming\sub-02\nirseeg\sub-02_rec-02_nirseeg.edf';
% oxyfile = 'C:\Users\maria\Universidade do Porto\João Pedro Barbosa Fonseca - Internship\After Experiment\pilots\pilot 1\sub-01_rec-01_nirseeg.edf';
% oxyfile = 'C:\Users\catar\OneDrive - Universidade do Porto\Internship\After Experiment\pilots\pilot 1\sub-01_rec-01_nirseeg.edf'

% fprintf('\n \n Select the EEGo file \n \n')
% [baseName, folder, idx] = uigetfile({'*.eeg;*.vhdr;*.vmrk', 'All EEGo file(*.eeg, *.vhdr, *.vmrk)';},'Select the EEGo file');
% eegfile = fullfile(folder, baseName);
% eegfile = 'C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\After Experiment\pilots\pilot 1\pilotfnirs_01.vhdr';
% eegfile = 'C:\Users\maria\Universidade do Porto\João Pedro Barbosa Fonseca - Internship\After Experiment\pilots\pilot 1\pilotfnirs_01.vhdr';
% eegfile = 'C:\Users\catar\OneDrive - Universidade do Porto\Internship\After Experiment\pilots\pilot 1\pilotfnirs_01.vhdr'


correct=input('Load data (Y/N)? If not, its assumed the .set files have been generated \n', 's');
if strcmpi(correct, 'y')
    done=1;
    data_loaded=0;
elseif strcmpi(correct, 'n')
    data_loaded=1;
end

%% Read data 
if data_loaded==0
    %FIELDTRIP - load the eeg&nirs data
    cfg = [];
    cfg.dataset = oxyfile;
    [data_raw] = ft_preprocessing(cfg);
    eeg_fnirs_events=ft_read_event(cfg.dataset);
    save(['sub-',sub,'_rec-',rec,'_nirseeg.mat'], 'data_raw');
    save(['sub-',sub,'_rec-',rec,'_nirseeg_events.mat'], 'eeg_fnirs_events');
    
    %EEGLAB load eeg only data
    [EEG,~]         = pop_loadbv(fullfile(sub_path,'eeg'), sub_vhdr);
    [ALLEEG,EEG,~]  = pop_newset(ALLEEG, EEG, 1,'setname','rawData','gui','off');

else %if data has been loaded and the datasets created, load the structs
    load(['sub-',sub,'_rec-',rec,'_nirseeg.mat']);    %Avoids call to ft_preprocessing
    load(['sub-',sub,'_rec-',rec,'_nirseeg_events.mat']);%Avoids call to ft_readevents
    
    [EEG]  = pop_loadset(['sub-',sub,'_rec-',rec,'_eeg.set'],fullfile(sub_path,'eeg'));
end

%% Read stimuli results
results=load(fullfile(sub_path,'stim',['results_sub-',sub,'_rec-',rec]));

%% VERIFY MARKERS

%Get nirseeg numerical markers
time=data_raw.time{1,1};
eeg_fnirs_markers=zeros(1,data_raw.sampleinfo(2));
for i=1:length(eeg_fnirs_events)
    if ~isempty(eeg_fnirs_events(i).value)
        s=eeg_fnirs_events(i).value;
        s=s(end-4:end);
        eeg_fnirs_markers(eeg_fnirs_events(i).sample)=str2double(s);
    end
end

%Get eeg only numerical markers
eeg_events=EEG.event;
eeg_markers=zeros(1,EEG.pnts);
for i=2:length(EEG.event)
    a=EEG.event(i).type;
    eeg_markers(EEG.event(i).latency)=str2double(a(end-3:end));
end


nonzero=find(eeg_fnirs_markers~=0);
start_eeg_fnirs=nonzero(1);
stop_eeg_fnirs=nonzero(end);
eeg_fnirs_markers=eeg_fnirs_markers(start_eeg_fnirs:stop_eeg_fnirs);

nonzero=find(eeg_markers~=0);
start_eeg=nonzero(1);
stop_eeg=nonzero(end);
eeg_markers=eeg_markers(start_eeg:stop_eeg);

%% GET COLUMN OF MARKERS ONLY
tmp=struct2table(eeg_fnirs_events);
value=tmp.value;
eeg_fnirs_markers_only=value(~cellfun('isempty',value));

tmp=struct2table(eeg_events);
eeg_markers_only=tmp.type;

%% Get labels only
eeg_fnirs_lbl=zeros(size(eeg_fnirs_markers_only));

for i=1:length(eeg_fnirs_markers_only)
    a=eeg_fnirs_markers_only(i);
    a=a{:};
    eeg_fnirs_lbl(i)=str2double(a(end-3:end));
end

eeg_lbl=zeros(size(eeg_markers_only));
for i=1:length(eeg_markers_only)
    a=eeg_markers_only(i);
    a=a{:};
    eeg_lbl(i)=str2double(a(end-3:end));
end

%% Get count of each marker
clear marker_summary
lbl=1698:1717;
lbl(end+1)=1777;
lbl(end+1)=1255;
lbl(end+1)=1500;
lbl(end+1)=1555;
lbl(end+1:end+3)=1797:1799;
lbl_numeric=lbl;
n=0;
for i=lbl_numeric
    n=n+1;
    eeg_no=sum(eeg_lbl==i);
    eeg_fnirs_no=sum(eeg_fnirs_lbl==i);
    marker_summary(1,n)=i;
    marker_summary(2,n)=eeg_no;
    marker_summary(3,n)=eeg_fnirs_no;
end

lbl = sprintfc('%d',lbl_numeric);
marker_table=array2table(marker_summary);
marker_table.Properties.VariableNames={'StartMetronome','StopMetronome','StartCue','StopCue','StartAutoCue','StartAutoNoCue','StartNonAutoCue','StartNonAutoNoCue','StartAutoDualCue','StartAutoDualNoCue','StartNonAutoDualCue','StartNonAutoDualNoCue','StopAutoCue','StopAutoNoCue','StopNonAutoCue','StopNonAutoNoCue','StopAutoDualCue','StopAutoDualNoCue','StopNonAutoDualCue','StopNonAutoDualNoCue','Key','CheckFlip','CheckStart','CheckStop','Letter','StartMovie','StopMovie'};

%% Extract Task Data
[EEG, starts, stops]=extractTaskData_EEG(EEG,marker_table, results, file);
[ALLEEG,EEG,~]  = pop_newset(ALLEEG, EEG, 1,'setname','taskData','gui','off');
pop_saveset( EEG, ['sub-',sub,'_rec-',rec,'_eeg_task.set'],fullfile(sub_path,'eeg'));


%% Filter EEG - 50 Hz noise and harmonics

% Determine the power spectrum of the raw data.
raw_data = EEG.data;
[P_raw, f_raw] = periodogram(raw_data', [], [] , EEG.srate);

% Filter the signal.
if ~isfile(file.filtered) 
    EEG_filtered = EEG;
    filtered_data = filterNoise(double(raw_data), EEG, 4);
    EEG_filtered.data = filtered_data;
    [ALLEEG, EEG_filtered, ~] = pop_newset(ALLEEG, EEG_filtered, 1,...
        'setname','filtData','gui','off');
    save(file.filtered, 'EEG_filtered');
else
    load(file.filtered, 'EEG_filtered');
    filtered_data = EEG_filtered.data;
    [ALLEEG, EEG_filtered, ~] = pop_newset(ALLEEG, EEG_filtered, 1,...
        'setname', 'filtData', 'gui', 'off');
end

% Determine the power spectrum of the filtered data.
[P_filt, f_filt] = periodogram(filtered_data', [], [] , EEG.srate);

% Plot the power spectrums.
figure;
subplot(1, 3, 1); plot(f_raw, P_raw); 
xlim([0 200]); ylim([0 7e5]); title('Raw data');
subplot(1, 3, 2); plot(f_filt, P_filt); 
xlim([0 200]); ylim([0 7e5]); title('Filtered data - same scale');
subplot(1, 3, 3); plot(f_filt, P_filt); 
xlim([0 50]); title('Filtered data - different scale');

%% Removal of eye blinks - preICA

if ~isfile(file.preICA)  
    [EEG_preICA] = pop_runica(EEG_filtered,'icatype', 'runica',...
        'extended', 1, 'interrupt', 'on');
    [ALLEEG, EEG_preICA, ~] = pop_newset(ALLEEG, EEG_preICA, 1,...
        'setname', 'preICA', 'gui', 'off');
    save(file.preICA, 'EEG_preICA');
else
    load(file.preICA, 'EEG_preICA');
    [ALLEEG, EEG_preICA, ~] = pop_newset(ALLEEG, EEG_preICA, 1,...
        'setname', 'preICA', 'gui', 'off');
end

%% Removal of eye blinks - pstICA

if ~isfile(file.pstICA)
    [EEG_preICA] = pop_chanedit(EEG_preICA, 'lookup', join([eeglab_path, '\\plugins\\dipfit\\standard_BESA\\standard-10-5-cap385.elp']), 'lookup', join([eeglab_path, '\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc']));
    [EEG_pstICA] = run_postICA(EEG_preICA);
    [ALLEEG, EEG_pstICA, ~] = pop_newset(ALLEEG, EEG_pstICA, 1,...
        'setname', 'pstICA','gui','off');
    save(file.pstICA, 'EEG_pstICA');
else                          
    load(file.pstICA, 'EEG_pstICA');
    [ALLEEG, EEG_pstICA, ~] = pop_newset(ALLEEG, EEG_pstICA, 1,...
        'setname', 'pstICA','gui','off');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
