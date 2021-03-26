%% Lab Streaming Layer

lib     = lsl_loadlib();                      % Load Lab Streaming Layer Library
info    = lsl_streaminfo(lib, 'ReactionTimeTest', 'Markers', 1, 0.0, 'cf_int32', 'ReactionTime'); 
outlet  = lsl_outlet(info);                   % Open outlet

% % % Marker Setup % % %
% Blocks Related
Marker_StartBlockCue1HzAddition       = 6000;        
Marker_StartBlockCue3HzAddition       = 6001;       

Marker_EndBlockCue1HzAddition         = 6002;
Marker_EndBlockCue3HzAddition         = 6003;

% Sample Related
Marker_GoStimulusAddition             = 6004;      

% Test Marker
Marker_TestTest                      = 5555;     
                                        

% General Screen Preparation

% Default Setup for Psychtoolbox
PsychDefaultSetup(2);

% Find Screen Numbers
screens     = Screen('Screens');
ScreenNr    = max(screens);

% Set Level Of Priority
PriorityLevel = MaxPriority(ScreenNr);
Priority(PriorityLevel);

% Defenition of Colors
white       = WhiteIndex(ScreenNr);
black       = BlackIndex(ScreenNr);
grey        = white/2;

%% Prepare Screen Window

% Open A Window (With Specific Color)
[Window, WindowRect] = PsychImaging('OpenWindow', ScreenNr, black);

% Find and Set the Center Coordinate of Window
[Xpixel, Ypixel]     = Screen('WindowSize', Window);
[Xcenter, Ycenter]   = RectCenter(WindowRect);

% Anti-Aliasing (Enable Alpha Blending)
Screen('BlendFunction', Window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% Find Flip-Interval (Refresh Rate of Screen)
ifi = Screen('GetFlipInterval', Window);

%% Texture Generation

% Text Please Wait
PWLine1     = 'Initializing... \n \n Please wait for welcome screen';
DrawFormattedText(Window, PWLine1, 'Center', 'Center', white);
Screen('Flip', Window);
WaitSecs(3);

% Create Rectangle (Size 60% of the Screen)
BaseRect    = [0 0 .6*Xpixel 0.6*Ypixel];
Centered    = CenterRectOnPointd(BaseRect, Xcenter, Ycenter);

% Focus Dot Settings
DotSize     = 20;                   % Pixels

% Go Stimulus
ColorRect   = [0 1 0];              % Green Rectangle
DotColor    = [1 0 0];              % Red Dot

% Show Texture
Screen('FillRect', Window, ColorRect, Centered);
% Screen('DrawDots', Window, [Xcenter, Ycenter], DotSize, DotColor, [], 2);
Screen('Flip', Window);

% Save Texture as Go Stimulus in Texture(1)
GoStimulus  = Screen('GetImage', Window);
Texture(1)  = Screen('MakeTexture', Window, GoStimulus);

%% SoundFiles

% Get Cueing Files
Cue1Hz       = 'Cue_075bps_350Hz.wav';
[Cue1Hz]     = CreateWAVstruct(Cue1Hz);
Cue1HzLength = length(Cue1Hz.wavedata)/Cue1Hz.fs;

Cue3Hz       = 'Cue_3bps_350Hz.wav';
[Cue3Hz]     = CreateWAVstruct(Cue3Hz);
Cue3HzLength = length(Cue3Hz.wavedata)/Cue3Hz.fs;

%% General Sound Preparation

% PsychPortAudio
InitializePsychSound(1);
Priority = 0;

% Get Audio Device
PPA_device = PsychPortAudio('GetDevices');

% Set Handle For Audio Playback
PPA_cue1Hz = PsychPortAudio('Open', [], [], Priority, Cue1Hz.fs, Cue1Hz.nrChan);
PPA_cue3Hz = PsychPortAudio('Open', [], [], Priority, Cue3Hz.fs, Cue3Hz.nrChan);
file       = [PPA_cue1Hz; PPA_cue3Hz];

% Load Buffer For Audio Playback
PsychPortAudio('FillBuffer', PPA_cue1Hz, Cue1Hz.wavedata);
PsychPortAudio('FillBuffer', PPA_cue3Hz, Cue3Hz.wavedata);

%% Experiment

outlet.push_sample(Marker_TestTest);
WaitSecs(5);
outlet.push_sample(Marker_TestTest);

% Baseline Period
Screen('DrawDots', Window, [Xcenter, Ycenter], DotSize, DotColor, [], 2);
Screen('Flip',Window);
WaitSecs(10);

% AEP 1HZ
PsychPortAudio('Start', file(1), 1, [], []);
outlet.push_sample(Marker_StartBlockCue1HzAddition);

for ii = 1:60
    WaitSecs(1/0.75);
    outlet.push_sample(Marker_GoStimulusAddition);
end

PsychPortAudio('Stop', file(1));
outlet.push_sample(Marker_EndBlockCue1HzAddition);

WaitSecs(5);

% AEP 3HZ
PsychPortAudio('Start', file(2), 1, [], []);
outlet.push_sample(Marker_StartBlockCue3HzAddition);

for ii = 1:120
    WaitSecs(1/3);
    outlet.push_sample(Marker_GoStimulusAddition);
end

PsychPortAudio('Stop', file(2));
outlet.push_sample(Marker_EndBlockCue3HzAddition);

WaitSecs(5);

% VEP
for ii = 1:60
    Screen('DrawTexture',Window,Texture(1));
    Screen('Flip',Window);
    outlet.push_sample(Marker_GoStimulusAddition);
    WaitSecs(0.75);
    
    Screen('DrawDots', Window, [Xcenter, Ycenter], DotSize, DotColor, [], 2);
    Screen('Flip',Window);
    
    WaitSecs(2);
end

sca
        
%% Functions

% To Play Back Sound
function [WAVstruct] = CreateWAVstruct(WAVfilename)
% This function creates a struct with the information from the wav-files.

    wav = WAVfilename;                                          
    WAVstruct = struct('wavedata',[],'fs',[],'nrChan',[]);      
    [WAVstruct.wavedata, WAVstruct.fs] = psychwavread(wav);     
    WAVstruct.wavedata = WAVstruct.wavedata';                   
    WAVstruct.nrChan = size(WAVstruct.wavedata,1);  
end