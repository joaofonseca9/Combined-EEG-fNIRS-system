%% MAIN

% This script is used to set up the experiment, randomize the order of the
% tasks

%% Settings:

% Clear the workspace and the screen
%sca; close all; clear all; clc


%Before starting the automaticity test, clear the workspace.
clear all

%Synch test skip => comment when actually testing patient
Screen('Preference', 'SkipSyncTests', 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SET UP PARAMETERS

%The sequences used for this study (automatic and non-automatic sequences
%randomized between participants)

%Sequences used in order to be able to print in the command window if
%to generate a new sequence use randi([1 4], 1, 12)
sequencesprint = {('4 3 4 1 4 1 2 4 3 2 1 2'),('2 1 2 3 2 1 3 2 4 2 4 1')};

sequences = {split(sequencesprint(1))',split(sequencesprint(2))'} ;

%Parameters for the resting period in between the trials
t1 = 20; %Resting period in seconds
t2 = 5;  %Random interval around the resting period time
t3 = 9.5; %Duration of a trial (tapping the sequence 1 time)

%Amount of letters presented during test for automaticity for one trial.
%Should be adjusted when letter presenting speed is changed!
N_letters=8; % 8 letters presented during a trial
N_trials=2; % number of trials per block

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LSL SETUP
% LSL outlet sending events

% Create and load the lab streaming layer library


addpath(genpath('C:\Users\joaop\Downloads\liblsl-Matlab'));
% addpath(genpath('C:\Users\catar\Downloads\liblsl-Matlab-master'));
% addpath(genpath('C:\Users\maria\OneDrive\Documentos\GitHub\liblsl-Matlab'));
lib = lsl_loadlib(); version = lsl_library_version(lib);
lib = lsl_loadlib();

% Make a new stream outlet.
% info = lsl_streaminfo([lib handle],[name],[type],[channelcount],[fs],[channelformat],[sourceid])
% > name = name of stream; describes device/product
% > type = content type of stream (EEG, Markers)
% > channelcount = nr of channels per sample
% > fs = samplking rate (Hz) as advertized by data source
% > channelformat = cf_float32, cf__double64, cf_string, cf_int32, cf_int16
% > sourceid = unique identifier for source or device, if available
info    = lsl_streaminfo(lib, 'Dual Task', 'Markers', 1, 0.0, 'cf_int32', 'Automaticity_DualTask'); 

% Open an outlet for the data to run through.
outlet = lsl_outlet(info);

%% MARKER SETUP
% Block related
% instructions = 'instructions'; %NEVER USED (?)
% finger_test='finger_test';

Marker_StartBlock_Cue       = 1700;         
Marker_EndBlock_Cue         = 1701;

Marker_StartBlock_AutomaticSequence     = 1702;
Marker_StartBlock_NonAutomaticSequence  = 1703;

Marker_StartBlock_AutomaticSequence_Dual     = 1704;
Marker_StartBlock_NonAutomaticSequence_Dual  = 1705;

Marker_EndBlock_AutomaticSequence       = 1712;
Marker_EndBlock_NonAutomaticSequence       = 1713;

Marker_EndBlock_AutomaticSequence_Dual      = 1714;
Marker_EndBlock_NonAutomaticSequence_Dual      = 1715;

Marker_CHECK = 1255;        % checkerboard flip
Marker_start = 1555;        % start signal 
Marker_stop = 1500;         % stop signal 

%Open Pshychtoolbox.
PsychDefaultSetup(2);
KbName('UnifyKeyNames'); %Links the key presses to the key board names
KbQueueCreate;
KbQueueStart; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD METRONOME SOUNDS (PsychToolbox)
audio_dir='.\metronomesounds';
addpath(audio_dir)
[WAVMetronome8.wave,WAVMetronome8.fs]       = audioread('Metronome8.wav');
[WAVMetronome600.wave,WAVMetronome600.fs]       = audioread('Metronome600.wav');
[WAVMetronome300.wave,WAVMetronome300.fs]       = audioread('Metronome300.wav');

% change rows<>columns
WAVMetronome8.wave = WAVMetronome8.wave';         WAVMetronome8.nrChan=2;
WAVMetronome600.wave = WAVMetronome600.wave';         WAVMetronome600.nrChan=2;
WAVMetronome300.wave = WAVMetronome300.wave';         WAVMetronome300.nrChan=2;

% Get Cueing Files
Cue1_25Hz       = 'Metronome120.wav';
[Cue1_25Hz]     = CreateWAVstruct(Cue1_25Hz);
Cue1_25HzLength = length(Cue1_25Hz.wavedata)/Cue1_25Hz.fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CREATE AND FILL AUDIO BUFFER
% Initialize Sounddriver
% This routine loads the PsychPortAudio sound driver for high precision, low latency,
% multichannel sound playback and recording
% Call it at the beginning of your experiment script, optionally providing the
% 'reallyneedlowlatency'-flag set to 1 to push really hard for low latency
InitializePsychSound(1);

priority = 0;                       % 0 = better quality, increased latency; 1 = minimum latency
duration = 1;                       % number of repetitions of the wav-file
PsychPortAudio('Verbosity',1);      % verbosity = "wordiness" -> 1= print errors

% Get audio device
h_device = PsychPortAudio ('GetDevices');

% Open handle
h_Metronome8   = PsychPortAudio('Open', [], [], priority, WAVMetronome8.fs, WAVMetronome8.nrChan);
h_Metronome600   = PsychPortAudio('Open', [], [], priority, WAVMetronome600.fs, WAVMetronome600.nrChan);
h_Metronome300   = PsychPortAudio('Open', [], [], priority, WAVMetronome300.fs, WAVMetronome300.nrChan);
PPA_cue1_25Hz = PsychPortAudio('Open', [], [], priority, Cue1_25Hz.fs, Cue1_25Hz.nrChan);

% Fill buffer
PsychPortAudio('FillBuffer', h_Metronome8, WAVMetronome8.wave);
PsychPortAudio('FillBuffer', h_Metronome600, WAVMetronome600.wave);
PsychPortAudio('FillBuffer', h_Metronome300, WAVMetronome300.wave);
PsychPortAudio('FillBuffer', PPA_cue1_25Hz, Cue1_25Hz.wavedata);

%AudioFile
file = [h_Metronome8; PPA_cue1_25Hz; h_Metronome300; h_Metronome600];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE FILES IN FOLDER

fprintf('Select the project directory \n')
root_dir=uigetdir('C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\Combined-EEG-fNIRS-system', 'Select the project directory');
addpath(root_dir);
complete=0;
while complete==0
    sub_ID=input('What is the subject ID (2 digit number) \n', 's');
    sub=sprintf('sub-%s', sub_ID);
        rec_n=input('What is the number of the recording? \n');
        rec=sprintf('rec-%.2d', rec_n);

    inf=fprintf('\n root_dir = %s \n sub = %s \n rec = %s \n', root_dir, sub, rec);
    correct=input('Is the above information correct? (y/n) \n', 's');
    if strcmp(correct, 'y')
        complete=1;
    else
        continue
    end
end

% go to subject folder
sub_dir=fullfile(root_dir, sub);
if ~exist(sub_dir)
    mkdir(sub_dir)
end
cd(sub_dir)

logname=sprintf('%s_%s_triggers.log', sub, rec); diary(logname);
% save current script in subject directory
script=mfilename('fullpath');
script_name=mfilename;
copyfile(sprintf('%s.m', script), fullfile(sub_dir, sprintf('%s_%s.m', sub, script_name)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SCREEN PREPARATION

% Get the screen numbers.
screens = Screen('Screens');

% Select the external screen if it is present, else revert to the native
% screen
screenNumber = max(screens);

% Define black, white and grey
black = BlackIndex(screenNumber);
white = WhiteIndex(screenNumber);

% Open an on screen window and color it grey
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, black);

% Get the size of the on screen window in pixels
% For help see: Screen WindowSize?
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Get the centre coordinate of the window in pixels
% For help see: help RectCenter
[xCenter, yCenter] = RectCenter(windowRect);

% Query the frame duration
ifi = Screen('GetFlipInterval', window);

% Set up alpha-blending for smooth (anti-aliased) lines
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% preparations for the fixation cross so we only need to do this once
fixCrossDimPix = 40; % Here we set the size of the arms of our fixation cross
xCoords = [-fixCrossDimPix fixCrossDimPix 0 0]; % Set the coordinates (these are all relative to zero we will let the drawing routine center the cross in the center of our monitor for us)
yCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
allCoords = [xCoords; yCoords];
lineWidthPix = 4;% Set the line width for the fixation cross
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save result structs for later
% str=['events_autosingle',sub,'_',rec];
% save(str,'events_autosingle');
% 
% str=['events_autodual',sub,'_',rec];
% save(str,'events_autodual');
% 
% str=['events_nonautosingle',sub,'_',rec];
% save(str,'events_nonautosingle');
% 
% str=['events_nonautodual',sub,'_',rec];
% save(str,'events_nonautodual');

%% PSEUDORANDOMIZATION

%Randomize the first sequence to be tested
%1=auto sequence, 2 non-auto sequence
order_sequence=randperm(2,2);

%Pseudorandomize which trials are cued and uncued (must be 50/50 split)
events_autosingle=randCuedTrials(N_trials);
events_autodual=randCuedTrials(N_trials);
events_nonautosingle=randCuedTrials(N_trials);
events_nonautodual=randCuedTrials(N_trials);

%% WELCOME SCREEN
%Instruction automaticity test
Screen('TextSize',window,45);
DrawFormattedText(window, 'Welcome to the experiment! \n \n If you have any questions please ask now. \n \n Thank you for participating!','center', 'center', white);
vbl = Screen('Flip', window);
WaitSecs(10);


%% Demonstrate how cueing works
Screen('TextSize',window,25);
DrawFormattedText(window,'While you perform the tasks, \n you may hear a rhytmic sound like this one. \n Or you may just hear nothing.','center', 'center', white);
vbl = Screen('Flip', window);
PsychPortAudio('Start', file(2), 1, [], []);
WaitSecs(5);
PsychPortAudio('Stop', file(2));

Exp_Script_Automaticity;

for sequence_idx=order_sequence
    if sequence_idx==1
        Exp_Script_AutoSingle
    else
        Exp_Script_NonAutoSingle
    end
end

Exp_Script_Checkerboard;

%% HELPER FUNCTIONS
function events=randCuedTrials(n)
    %Generate a vector where half is 0's and half is 1's
    %1=cued, 0=uncued
    isCued = zeros(n, 1);
    isCued(randperm(numel(isCued), round(n/2))) = 1;
    for j=1:n
        events.trial(j).cue=isCued(j);
    end
end

% To Play Back Sound
function [WAVstruct] = CreateWAVstruct(WAVfilename)
% This function creates a struct with the information from the wav-files.

    wav = WAVfilename;                                          
    WAVstruct = struct('wavedata',[],'fs',[],'nrChan',[]);      
    [WAVstruct.wavedata, WAVstruct.fs] = psychwavread(wav);     
    WAVstruct.wavedata = WAVstruct.wavedata';                   
    WAVstruct.nrChan = size(WAVstruct.wavedata,1);  
end