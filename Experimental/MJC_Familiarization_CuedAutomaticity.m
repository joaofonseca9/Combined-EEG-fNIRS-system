% FAMILIARIZATION SCRIPT
% Presenting the studied sequence on the screen for 2 minutes, while a
% metronome sound of 75bpm (1.25Hz) is playing. 2 min practicing for 
% the finger tapping, for the participant to get used to the experiment 
% set up. 

%% SETTINGS:

% Clear the workspace and the screen
clear;

% Skip screen synchronization to prevent Pyshtoolbox for freezing
Screen('Preference', 'SkipSyncTests', 1);

% root_dir='C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\Experiment\Combined-EEG-fNIRS-system\Experimental';
% root_dir='C:\Users\maria\OneDrive\Documentos\GitHub\Combined-EEG-fNIRS-system\Experimental';
root_dir='C:\Users\catar\OneDrive - Universidade do Porto\Internship\Experiment\Combined-EEG-fNIRS-system\Experimental';

%% SET RIGHT AUTOMATIC SEQUENCE
sequenceA = '2 4 3 4 1 3 4 1 2 1 3 2';
sequenceB = '4 1 3 2 4 1 4 2 3 2 1 3';
sequenceauto = sequenceA;

%% Open Phsychtoolbox
PsychDefaultSetup(2);
KbName('UnifyKeyNames'); % Links the key presses to the key board names
KbQueueCreate;
KbQueueStart; 

%% LOAD METRONOME SOUNDS (PsychToolbox)
audio_dir=fullfile(root_dir, 'metronomesounds');
addpath(audio_dir)
[WAVMetronome8.wave,WAVMetronome8.fs] = audioread('Metronome8.wav');
[WAVMetronome300.wave,WAVMetronome300.fs] = audioread('Metronome300.wav');

% Change rows<>columns
WAVMetronome8.wave = WAVMetronome8.wave';         WAVMetronome8.nrChan=2;
WAVMetronome300.wave = WAVMetronome300.wave';         WAVMetronome300.nrChan=2;

% Get Cueing Files
Cue1_25Hz = 'Metronome120.wav';
[Cue1_25Hz] = CreateWAVstruct(Cue1_25Hz);
Cue1_25HzLength = length(Cue1_25Hz.wavedata)/Cue1_25Hz.fs;

%% CREATE AND FILL AUDIO BUFFER
% Initialize Sounddriver
% This routine loads the PsychPortAudio sound driver for high precision, low latency,
% multichannel sound playback and recording
% Call it at the beginning of your experiment script, optionally providing the
% 'reallyneedlowlatency'-flag set to 1 to push really hard for low latency
InitializePsychSound(1);

priority = 0;                       % 0 = better quality, increased latency; 1 = minimum latency
duration = 1;                       % Number of repetitions of the wav-file
PsychPortAudio('Verbosity',1);      % Verbosity = "wordiness" -> 1= print errors

% Get audio device
h_device = PsychPortAudio ('GetDevices');

% Open handle
h_Metronome8   = PsychPortAudio('Open', [], [], priority, WAVMetronome8.fs, WAVMetronome8.nrChan);
h_Metronome300   = PsychPortAudio('Open', [], [], priority, WAVMetronome300.fs, WAVMetronome300.nrChan);
PPA_cue1_25Hz = PsychPortAudio('Open', [], [], priority, Cue1_25Hz.fs, Cue1_25Hz.nrChan);

% Fill buffer
PsychPortAudio('FillBuffer', h_Metronome8, WAVMetronome8.wave);
PsychPortAudio('FillBuffer', h_Metronome300, WAVMetronome300.wave);
PsychPortAudio('FillBuffer', PPA_cue1_25Hz, Cue1_25Hz.wavedata);

% % Set Handle For Audio Capture (delay check)
% CAP_cue1_25Hz = PsychPortAudio('Open', [], 2, priority, Cue1_25Hz.fs, Cue1_25Hz.nrChan);

%AudioFile
file = [h_Metronome8; PPA_cue1_25Hz; h_Metronome300];

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

%% START FAMILIARIZATION
HideCursor;
% Instruction familiarization
Screen('TextSize',window,25);
DrawFormattedText(window,'In the upcoming 5 minutes you have time to get familiarized with the experimental set up. \n\n You have 5 minutes to practice the finger tapping task. \n\n The sequence will be shown on the screen. \n\n You will also hear a metronome sound at times, which is at the same speed as you studied at home. \n\n Press any key to start practicing the sequence for finger tapping.','center', 'center', white);
vbl = Screen('Flip', window);
KbStrokeWait; 

% Presenting the sequence to study on the screen
Screen('TextSize', window, 50);
DrawFormattedText(window, sprintf('%s', sequenceauto), 'center', 'center', white);
vbl= Screen('Flip', window);
WaitSecs(180)
PsychPortAudio('Start', file(2), 1, [], []); %Play metronome sound file (2 minutes)
WaitSecs(120)
PsychPortAudio('Stop', file(2));

%End of finger tapping familiarization
Screen('TextSize',window,25);
DrawFormattedText(window,'The time of the familiarization is over. \n\n Press any key to end the familiarization.','center', 'center', white);
vbl = Screen('Flip', window);
WaitSecs(5)
KbStrokeWait; 
sca

%% HELPER FUNCTIONS 
% To Play Back Sound
function [WAVstruct] = CreateWAVstruct(WAVfilename)
% This function creates a struct with the information from the wav-files.

    wav = WAVfilename;                                          
    WAVstruct = struct('wavedata',[],'fs',[],'nrChan',[]);      
    [WAVstruct.wavedata, WAVstruct.fs] = psychwavread(wav);     
    WAVstruct.wavedata = WAVstruct.wavedata';                   
    WAVstruct.nrChan = size(WAVstruct.wavedata,1);  
end
