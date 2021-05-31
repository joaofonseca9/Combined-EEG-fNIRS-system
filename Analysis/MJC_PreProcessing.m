clear;
close all;

%% Initialize FieldTrip & EEGLAB
laptop='laptopMariana';
% laptop='laptopJoao';
% laptop='laptopCatarina';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);

eeglab;
ft_defaults;

sub='03';
rec='02';

file = getFileNames(mainpath_out, sub, rec);

%% Select data folder 
sub_path = fullfile(mainpath_in,'incoming',['sub-',sub]);
eeg_path = fullfile(sub_path,'eeg');
nirseeg_path = fullfile(sub_path,'nirseeg');
sub_vhdr = fullfile(['sub-',sub,'_rec-',rec,'_eeg.vhdr']);

%Before changing directory to the subpath, add current directory to access
%the function files
addpath(pwd)
cd(sub_path);

oxyfile = fullfile(sub_path,'nirseeg',['sub-',sub,'_rec-',rec,'_nirseeg.edf']);

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

%% Verify markers
% Get nirseeg numerical markers
time=data_raw.time{1,1};
eeg_fnirs_markers=zeros(1,data_raw.sampleinfo(2));
for i=1:length(eeg_fnirs_events)
    if ~isempty(eeg_fnirs_events(i).value)
        s=eeg_fnirs_events(i).value;
        s=s(end-4:end);
        eeg_fnirs_markers(eeg_fnirs_events(i).sample)=str2double(s);
    end
end

% Get eeg only numerical markers
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

%% Get column of markers
tmp=struct2table(eeg_fnirs_events);
value=tmp.value;
eeg_fnirs_markers_only=value(~cellfun('isempty',value));

tmp=struct2table(eeg_events);
eeg_markers_only=tmp.type;

%% Get labels 
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
marker_table.Properties.VariableNames={'StartMetronome','StopMetronome','StartCue','StopCue','StartAutoCue','StartAutoNoCue','StartNonAutoCue','StartNonAutoNoCue','StartAutoDualCue','StartAutoDualNoCue','StartNonAutoDualCue','StartNonAutoDualNoCue','StopAutoCue','StopAutoNoCue','StopNonAutoCue','StopNonAutoNoCue','StopAutoDualCue','StopAutoDualNoCue','StopNonAutoDualCue','StopNonAutoDualNoCue','Key','CheckFlip','CheckStop','CheckStart','Letter','StartMovie','StopMovie'};



%% Load MNI coordinates
% Load channel coordinates/positions of the standard MNI model of eeglab: 
% MNI dipfit channel positions
[EEG] = pop_chanedit(EEG, 'lookup', join([eeglab_path,...
        '\\plugins\\dipfit\\standard_BESA\\standard-10-5-cap385.elp']),...
        'lookup', join([eeglab_path,...
        '\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc']));

%% Filter EEG - 50 Hz noise and harmonics
% Determine the power spectrum of the raw data
% raw_data = EEG.data;
% [P_raw, f_raw] = periodogram(raw_data', [], [] , EEG.srate);

% Filter the signal
if ~isfile(file.filtered) 
    filtered_data = filterNoise(double(raw_data), EEG, 4);
    EEG.data = filtered_data;
    [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1, 'setname', 'filtData',...
        'gui', 'off');
    save(file.filtered, 'EEG');
else
    load(file.filtered, 'EEG');
    filtered_data = EEG.data;
    [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1, 'setname', 'filtData',...
        'gui', 'off');
end

% Determine the power spectrum of the filtered data
% [P_filt, f_filt] = periodogram(filtered_data', [], [] , EEG.srate);

% Plot the power spectrums
% figure;
% subplot(1, 3, 1); plot(f_raw, P_raw); 
% xlim([0 200]); ylim([0 7e5]); title('Raw data');
% subplot(1, 3, 2); plot(f_filt, P_filt); 
% xlim([0 200]); ylim([0 7e5]); title('Filtered data - same scale');
% subplot(1, 3, 3); plot(f_filt, P_filt); 
% xlim([0 50]); title('Filtered data - different scale');

%% Remove bad channels 
% Visually inspect the signals and choose if a signals is too bad that it
% needs to be removed.
% First see the power spectrum and then check if the signal is actually bad
% on the plot.
if ~isfile(file.removedBadChannels) 
    figure; 
    pop_spectopo(EEG, 1, [0 EEG.pnts], 'EEG', 'percent', 50, 'freqrange',...
        [2 75], 'electrodes', 'off');
    pop_eegplot(EEG);
    RC = input('Remove channel [nr/no]: ','s');
    while ~strcmp(RC, 'no')
        [EEG] = pop_select(EEG, 'nochannel', eval(RC));
        figure;
        pop_spectopo(EEG, 1, [0 EEG.pnts], 'EEG', 'percent', 50,...
            'freqrange', [2 75], 'electrodes', 'off');
        pop_eegplot(EEG);
        RC = input('Remove channel [nr/no]: ','s');
    end
    [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1, 'setname',...
        'removedBadChannels', 'gui', 'off');
    save(file.removedBadChannels, 'EEG');
else
    load(file.removedBadChannels, 'EEG');
    [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1, 'setname',...
        'removedBadChannels', 'gui', 'off');
end

%% Removal of eye blinks - preICA
if ~isfile(file.preICA)  
    [EEG] = pop_runica(EEG,'icatype', 'runica', 'extended', 1,...
        'interrupt', 'on');
    [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1, 'setname', 'preICA',...
        'gui', 'off');
    save(file.preICA, 'EEG');
else
    load(file.preICA, 'EEG');
    [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1, 'setname', 'preICA',...
        'gui', 'off');
end

%% Removal of eye blinks - pstICA
if ~isfile(file.pstICA)
    [EEG] = run_postICA(EEG);
    [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1, 'setname', 'pstICA',...
        'gui', 'off');
    save(file.pstICA, 'EEG');
else                          
    load(file.pstICA, 'EEG');
    [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1, 'setname', 'pstICA',...
        'gui', 'off');
end

%% Set reference
% Re-reference the system to Cz 

if ~isfile(file.preprocessed)
    [EEG] = pop_reref(EEG, 'Cz');
    [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1, 'setname',...
        'preprocessed', 'gui', 'off');
    save(file.preprocessed, 'EEG');
else                          
    load(file.preprocessed, 'EEG');
    [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1, 'setname',...
        'preprocessed', 'gui', 'off');
end

%% Extract task data
[EEG_divided, file]=extractTaskData_EEG(EEG,marker_table, results, file, mainpath_out);
save(file.EEG_divided,'EEG_divided');
%[ALLEEG,EEG,~]  = pop_newset(ALLEEG, EEG_task, 1,'setname','taskData','gui','off');

