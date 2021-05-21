%% REACTIONTIMETEST %%
clear all; close all; clc;

Subject     = '05';

Dir         = ['C:\Users\Helena\Documents\dryeeg\Go-No-Go\Results_struct\'];
Filename    = ['s' Subject '_GoCue.mat'];

Test = 0;

%% Lab Streaming Layer

Screen('Preference', 'SkipSyncTests', 1);

% lib     = lsl_loadlib();                      % Load Lab Streaming Layer Library
% info    = lsl_streaminfo(lib, 'ReactionTimeTest', 'Markers', 1, 0.0, 'cf_int32', 'ReactionTime'); 
% outlet  = lsl_outlet(info);                   % Open outlet

% % % Marker Setup % % %
% Test Related
Marker_StartTest              = 1000;         % Start of Experiment
Marker_EndTest                = 1001;         % End of Experiment

% Blocks Related
Marker_StartBlockCue1Hz       = 1100;         % Start of Block with Cues 1Hz
Marker_StartBlockCue3Hz       = 1200;         % Start of Block with Cues 3Hz
Marker_StartBlockNoCue        = 1300;         % Start of Block with no Cues
Marker_StartBlockGoOnly       = 1400;         % Start of Block with Go only

Marker_EndBlockCue1Hz         = 1101;         % End of Block with Cues 1Hz
Marker_EndBlockCue3Hz         = 1201;         % End of Block with Cues 3Hz
Marker_EndBlockNoCue          = 1301;         % End of Block with no Cues
Marker_EndBlockGoOnly         = 1401;         % End of Block with Go Only

% Sample Related
Marker_GoStimulus             = 2000;         % Go Stimulus is selected (Onset)
Marker_NoGoStimulus           = 2001;         % No Go Stimulus is selected (Onset)

Marker_Click                  = 3001;         % Mouse click to Stimulus (Offset)
Marker_NoClick                = 3002;         % No Reaction to Stimulus (Offset)

Marker_ErrorClick             = 3010;         % Clicked while not expected (error)
Marker_ErrorNoClick           = 3011;         % Not Clicked while expected (error)

% Test Marker
Marker_TestTest               = 5555;         % Test marker to see if markers are received and
                                              % to gather all other markers

%% General Screen Preparation

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
grey        = [102 102 102]./255;           % site 1 
grey        = [126 126 126]./255;           % Site 2 (Lum = +/- 50%)
grey        = [138 138 138]./255;

%% Prepare Screen Window

% Open A Window (With Specific Color)
[Window, WindowRect] = PsychImaging('OpenWindow', ScreenNr, grey);
HideCursor;

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

% Color settings
red         = [204 51 102]./255;        % Same Luminance as grey & green (+/- 50%)
% red         = [219 0 86]./255;          % Other option for red
green       = [0 153 102]./255;         % Same Luminance as grey & red (+/- 50%)

red = [250 78 37]./255;
green = [0 170 0]./255;

% Focus Dot
DotSize = 20;

% Create Rectangle (Size 60% of the Screen)
BaseRect    = [0 0 .6*Xpixel 0.6*Ypixel];
Centered    = CenterRectOnPointd(BaseRect, Xcenter, Ycenter);

% Show Texture
Screen('FillRect', Window, green, Centered);
Screen('Flip', Window);

% Save Texture as GoStimulus in Texture(1)
GoStimulus  = Screen('GetImage', Window);
Texture(1)  = Screen('MakeTexture', Window, GoStimulus);

% Wait momentarily before showing No Go Stimulus
WaitSecs(.5);    

% Show Texture
Screen('FillRect', Window, red, Centered);
Screen('Flip', Window);

% Save Texture as NoGoStimulus in Texture(2)
NoGoStimulus = Screen('GetImage', Window);
Texture(2)   = Screen('MakeTexture', Window, NoGoStimulus);

WaitSecs(.5);

% Return to initializing screen
DrawFormattedText(Window, PWLine1, 'Center', 'Center', white);
Screen('Flip', Window);

WaitSecs(3);

%% SoundFiles

% Get Cueing Files
Cue1Hz       = 'Cues_075.wav';
[Cue1Hz]     = CreateWAVstruct(Cue1Hz);
Cue1HzLength = length(Cue1Hz.wavedata)/Cue1Hz.fs;

Cue3Hz       = 'Cues_3.wav';
[Cue3Hz]     = CreateWAVstruct(Cue3Hz);
Cue3HzLength = length(Cue3Hz.wavedata)/Cue3Hz.fs;

%% General Sound Preparation

% PsychPortAudio
InitializePsychSound(1);
Priority = 1;
Cue1Hz.fs = [];
Cue3Hz.fs = [];

% Get Audio Device
PPA_device = PsychPortAudio('GetDevices');

% Set Handle For Audio Playback
PPA_cue1Hz = PsychPortAudio('Open', [], [], Priority, Cue1Hz.fs, Cue1Hz.nrChan);
PPA_cue3Hz = PsychPortAudio('Open', [], [], Priority, Cue3Hz.fs, Cue3Hz.nrChan);
file       = [PPA_cue1Hz; PPA_cue3Hz];

% Load Buffer For Audio Playback
PsychPortAudio('FillBuffer', PPA_cue1Hz, Cue1Hz.wavedata);
PsychPortAudio('FillBuffer', PPA_cue3Hz, Cue3Hz.wavedata);

% Set Handle For Audio Capture
CAP_cue1Hz = PsychPortAudio('Open', [], 2, Priority, 4*Cue1Hz.fs, Cue1Hz.nrChan);
CAP_cue3Hz = PsychPortAudio('Open', [], 2, Priority, 4*Cue3Hz.fs, Cue3Hz.nrChan);

%% Test Settings

% Set Sequence of Stimuli and Cued Blocks
[TestNr, StimSeq] = Sequence;
Results.TestNr    = TestNr;
Results.StimSeq   = StimSeq;
Results.StimShow  = zeros(size(StimSeq));

% Reset Amount of Errors
Results.Errors    = zeros(size(StimSeq));

% Settings to wait for mouse click
NumSecs           = 1.5;
NumFrames         = round(NumSecs / ifi);
WaitFrames        = 1;

%% Experiment

Screen('TextSize',Window,50);

WaitSecs(2);
% outlet.push_sample(Marker_TestTest);

    % Screen 1: Welcome Screen
    Line1 = 'Welcome, thank you for participating in this study! \n \n \n';
    Line2 = 'You will be performing a Go/No-Go task, \n \n while you are hearing repetitive, rhythmic sounds \n \n called cues.';
    
    DrawFormattedText(Window,[Line1 Line2],'Center','Center',white);
    DrawFormattedText(Window,'>>',Xpixel*.9,Ypixel*.9,white);
    Screen('Flip', Window);
    GetClicks(Window);
    
    % Screen 2: Introduction
    Line1 = 'During this experiment, you will be performing 19 blocks, with 20 stimuli each. \n \n \n';
    Line2 = 'The stimuli will be either a green or a red rectangle. \n \n \n';
    Line3 = 'In between the presentation of two stimuli, \n \n you have to focus on the red focus dot in the middle of the screen.';
    
    DrawFormattedText(Window,[Line1 Line2 Line3],'Center','Center',white);
    DrawFormattedText(Window,'>>',Xpixel*.9,Ypixel*.9,white);
    Screen('Flip', Window);
    GetClicks(Window);
    
    % Screen 3: Example Go Stimulus
    Line1 = 'If you see a green stimulus, \n \n you have to click the left mouse button as fast as possible.';
    
    DrawFormattedText(Window,Line1,'Center',.2*Ypixel,white);
    DrawFormattedText(Window,'>>',Xpixel*.9,Ypixel*.9,white);
    BaseRect = [0 0 .4*Xpixel .4*Ypixel];
    Centered = CenterRectOnPointd(BaseRect, Xpixel/2, 3*Ypixel/5);
    Screen('FillRect', Window, green, Centered);
    Screen('Flip', Window);
    GetClicks(Window);
    
    % Screen 4: Example No-Go Stimulus
    Line1 = 'If you see a red stimulus, \n \n you have to withhold your response by NOT clicking the mouse.';
    
    DrawFormattedText(Window,Line1,'Center',.2*Ypixel,white);
    DrawFormattedText(Window,'>>',Xpixel*.9,Ypixel*.9,white);
    BaseRect = [0 0 .4*Xpixel .4*Ypixel];
    Centered = CenterRectOnPointd(BaseRect, Xpixel/2, 3*Ypixel/5);
    Screen('FillRect', Window, red, Centered);
    Screen('Flip', Window);
    GetClicks(Window);
    
    % Screen 5: Go/No-Go vs Go only
    Line1 = 'Go';
    Line2 = '/';
    Line3 = 'No-Go';
    Line4 = 'task';
    Line5 = '\n \n In these blocks you will receive both green and red stimuli.';
    Line6 = '\n \n In these blocks you will receive only green stimuli.';
    Line7 = 'Before each block the screen will tell you the things you should expect.';
        
    DrawFormattedText(Window,Line1,.40*Xpixel,.3*Ypixel,green);
    DrawFormattedText(Window,Line2,.44*Xpixel,.3*Ypixel,white);
    DrawFormattedText(Window,Line3,.45*Xpixel,.3*Ypixel,red);
    DrawFormattedText(Window,Line4,.53*Xpixel,.3*Ypixel,white);
    DrawFormattedText(Window,Line5,'Center',.3*Ypixel,white);
    DrawFormattedText(Window,Line1,.45*Xpixel,.5*Ypixel,green);
    DrawFormattedText(Window,Line4,.49*Xpixel,.5*Ypixel,white);
    DrawFormattedText(Window,Line6,'Center',.5*Ypixel,white);
    DrawFormattedText(Window,Line7,'Center',.7*Ypixel,white);
    DrawFormattedText(Window,'>>',Xpixel*.9,Ypixel*.9,white);
    Screen('Flip', Window);
    GetClicks(Window);
    
    % Screen 6: Introduction Cues
    Line1 = 'During the performance of the Go/No-Go task, \n \n you will hear repetitive, rhythmic sounds called cues.';
    Line2 = '\n \n \n You do not have to do anything with these cues.';
    Line3 = '\n \n \n 3 different cueing conditions will be presented during the experiment:';
    
    DrawFormattedText(Window,[Line1 Line2 Line3],'Center','Center',white);
    DrawFormattedText(Window,'>>',Xpixel*.9,Ypixel*.9,white);
    Screen('Flip', Window);
    GetClicks(Window);    
    
    % Screen 7: Silence
    Line1= '1. Silence \n\n';
    
    DrawFormattedText(Window,Line1,'Center','Center',white);
    DrawFormattedText(Window,'>>',Xpixel*.9,Ypixel*.9,white);
    Screen('Flip',Window);
    GetClicks(Window);
    
    % Screen 8: Slow Cues
    Line1 = '2. Slow cues \n\n';
    Line2 = '[Click for an example]';
    
    DrawFormattedText(Window,[Line1 Line2],'Center','Center',white);
    Screen('Flip',Window);
    GetClicks(Window);
    
    DrawFormattedText(Window,Line1,'Center','Center',white);
    Screen('Flip',Window);
    
    PsychPortAudio('Start', file(1), 1, [], []);
    WaitSecs(5);
    PsychPortAudio('Stop',file(1));
    
    DrawFormattedText(Window,Line1,'Center','Center',white);
    DrawFormattedText(Window,'>>',Xpixel*.9,Ypixel*.9,white);
    Screen('Flip',Window);
    GetClicks(Window);
    
    % Screen 9: Fast Cues
    Line1 = '3. This sound \n\n';
    Line2 = '[Click for an example]';
    
    DrawFormattedText(Window,[Line1 Line2],'Center','Center',white);
    Screen('Flip',Window);
    GetClicks(Window);
    
    DrawFormattedText(Window,Line1,'Center','Center',white);
    Screen('Flip',Window);
    
    PsychPortAudio('Start', file(2), 1, [], []);
    WaitSecs(5.2);
    PsychPortAudio('Stop',file(2));
    
    DrawFormattedText(Window,Line1,'Center','Center',white);
    DrawFormattedText(Window,'>>',Xpixel*.9,Ypixel*.9,white);
    Screen('Flip',Window);
    GetClicks(Window);
    
    % Screen 10: Instruction Pause
    Line1 = 'In between blocks you can take a break as long as you prefer.\n\n';
    Line2 = 'You can look away or stretch a little.\n\n';
    Line3 = 'You can continue with the next block when you are ready and if you see >>. \n\n';
    Line4 = 'By clicking the mouse, the next block will start';
    
    DrawFormattedText(Window,[Line1 Line2 Line3 Line4],'Center','Center',white);
    DrawFormattedText(Window,'>>',Xpixel*.9,Ypixel*.9,white);
    Screen('Flip',Window);
    GetClicks(Window);
    
    % Screen 11: Questions?
    Line1 = 'If you have any questions, you can ask Robin \n\n';
    Line2 = 'Good Luck!';    
    
    DrawFormattedText(Window,[Line1 Line2],'Center','Center',white);
    DrawFormattedText(Window,'>>',Xpixel*.9,Ypixel*.9,white);
    Screen('Flip',Window);
    GetClicks(Window);
    
% Screens during experiment
    % Instruction screen (Go/No-Go)
    Line1 = 'Go';
    Line2 = '/';
    Line3 = 'No-Go';
    Line4 = 'task';
    Line5 = 'In this block you will receive green and red stimuli \n\n';
    Line6 = 'You will hear no cues';
    Line7 = 'You will hear slow cues';
    Line8 = 'You will hear fast cues';
    
    % Instruction screen (Go only)
    Line9 = 'In this block you will receive only green stimuli';
    
    % Baseline
    Line10= 'Sit still and do not speak \n \n';
    Line11= 'Focus on the dot \n \n \n \n \n';
    Line12= 'The next block will start in 10 seconds';
    
    % Break
    Line13= 'Take a short break!';

% outlet.push_sample(Marker_StartTest);

PsychPortAudio('GetAudioData',CAP_cue1Hz,120,[],[],[]);
PsychPortAudio('GetAudioData',CAP_cue3Hz,120,[],[],[]);

for iBlock = 1%:size(Results.StimSeq,2)
    BlockType = Results.TestNr(iBlock);
    BlockType = 1;
    switch BlockType
        case 1          % 1Hz Cues
            % Information Screen
            DrawFormattedText(Window,Line1,.1*Xpixel,.1*Ypixel,green);
            DrawFormattedText(Window,Line2,.14*Xpixel,.1*Ypixel,white);
            DrawFormattedText(Window,Line3,.15*Xpixel,.1*Ypixel,red);
            DrawFormattedText(Window,Line4,.23*Xpixel,.1*Ypixel,white);
            
            DrawFormattedText(Window,[Line5 Line7],'Center','Center',white);
            DrawFormattedText(Window,'>>',Xpixel*.9,Ypixel*.9,white);
            Screen('Flip',Window);
            GetClicks(Window);
            
            startrecord = GetSecs;
            PsychPortAudio('Start',CAP_cue1Hz,1, [], []);
            
            WaitSecs(1)

            
%             outlet.push_sample(Marker_StartBlockCue1Hz);
            startmoment = GetSecs;
            PsychPortAudio('Start', file(1), 1, [], []);

%             
%             [audiodata1Hz] = PsychPortAudio('GetAudioData',CAP_cue1Hz);
%             status1 = PsychPortAudio('GetStatus',CAP_cue1Hz);
              
        case 2          % 3Hz Cues
            % Information Screen
            DrawFormattedText(Window,Line1,.1*Xpixel,.1*Ypixel,green);
            DrawFormattedText(Window,Line2,.14*Xpixel,.1*Ypixel,white);
            DrawFormattedText(Window,Line3,.15*Xpixel,.1*Ypixel,red);
            DrawFormattedText(Window,Line4,.23*Xpixel,.1*Ypixel,white);
            
            DrawFormattedText(Window,[Line5 Line8],'Center','Center',white);
            DrawFormattedText(Window,'>>',Xpixel*.9,Ypixel*.9,white);
            Screen('Flip',Window);
            GetClicks(Window);
                        
%             outlet.push_sample(Marker_StartBlockCue3Hz);
            PsychPortAudio('Start', file(2), 1, [], []);
            
            PsychPortAudio('Start',CAP_cue3Hz,1, [], []);
%             [audiodata3Hz] = PsychPortAudio('GetAudioData',CAP_cue3Hz);
            
        case 3          % No Cues
            % Information Screen
            DrawFormattedText(Window,Line1,.1*Xpixel,.1*Ypixel,green);
            DrawFormattedText(Window,Line2,.14*Xpixel,.1*Ypixel,white);
            DrawFormattedText(Window,Line3,.15*Xpixel,.1*Ypixel,red);
            DrawFormattedText(Window,Line4,.23*Xpixel,.1*Ypixel,white);
            
            DrawFormattedText(Window,[Line5 Line6],'Center','Center',white);
            DrawFormattedText(Window,'>>',Xpixel*.9,Ypixel*.9,white);
            Screen('Flip',Window);
            GetClicks(Window);
            
%             outlet.push_sample(Marker_StartBlockNoCue);
        case 4          % Go Only Task
            % Information Screen
            DrawFormattedText(Window,Line1,.1*Xpixel,.1*Ypixel,green);
            DrawFormattedText(Window,Line4,.14*Xpixel,.1*Ypixel,white);
           
            DrawFormattedText(Window,[Line9 Line6],'Center','Center',white);
            DrawFormattedText(Window,'>>',Xpixel*.9,Ypixel*.9,white);
            Screen('Flip',Window);
            GetClicks(Window);
            
%             outlet.push_sample(Marker_StartBlockGoOnly);
    end
    
    % Baseline
    DrawFormattedText(Window,[Line10 Line11 Line12],'Center','Center',white);
    Screen('DrawDots', Window, [Xcenter, Ycenter], DotSize, red, [], 2);
    Screen('Flip',Window);
    WaitSecs(6);
    
    Screen('DrawDots', Window, [Xcenter, Ycenter], DotSize, red, [], 2);
    Screen('Flip',Window);
    WaitSecs(4);
    
    % % % Start Collection of Samples % % %
    for iRep = 1:size(Results.StimSeq,1)
        clear Button

        % % % Stimulus presentation % % %
        % Select Stimulus from provided sequence.
        Stimulus = Results.StimSeq(iRep, iBlock);

        % Present Stimulus after 2 to 4 seconds.
        WaitSecs(2+1.5*rand(1));                            % Randomly wait 2 to 4 seconds.
        Screen('DrawTexture', Window, Texture(Stimulus));   % Select Texture corresponding to stimulus.
        Results.StimShow(iRep, iBlock) = GetSecs;   
        Screen('Flip', Window);                                                    % Present the stimulus.
        
        % Push marker corresponding with stimulus and start counting for
        % reaction times
        switch Stimulus
            case 1
                Onset = GetSecs;
%                 outlet.push_sample(Marker_GoStimulus);
            case 2
                Onset = GetSecs;
%                 outlet.push_sample(Marker_NoGoStimulus);
        end
        
        % % % Reaction % % %
        % Wait until NumSecs have passed or until left mouse button is
        % clicked. Don't forget to correct for one more WaitSecs in analysis.
        Offset = 0;
        for frame = 1:NumFrames                             % Select how much frames the for loop takes.
            [~, ~, Button] = GetMouse(Window);
            if Button(1)                                    % If left mouse button is pressed. 
                Offset = GetSecs;
%                 outlet.push_sample(Marker_Click);                                            % Push the marker for the reaction time offset
                break                                       % and break out of for loop.
            end
            WaitSecs(WaitFrames * ifi);                     % Select how long one frame takes at most.
        end

        % Change Screen to only focus dot.
        Screen('DrawDots', Window, [Xcenter, Ycenter], DotSize, red, [], 2);
        Screen('Flip', Window);

        % % % Collect every result % % %
        % Calculate Reaction Time.
        ReactionTime = Offset - Onset;
        
        % Save Reaction Times and Errors.
        switch Stimulus                                         % Check which stimulus is presented
            case 1                                              % Go Stimulus
                if ReactionTime < NumSecs                       % Check if the reaction time is within the presentation of the stimulus
                    Results.ReactionTime(iRep, iBlock) = ReactionTime;
                    Results.Errors(iRep, iBlock)       = 0;     % If so, save reaction time and no errors.
                else
%                     outlet.push_sample(Marker_ErrorNoClick);    % If not so, Error No Click is presented
                    Results.ReactionTime(iRep, iBlock) = NaN;   % Reaction time is non existent
                    Results.Errors(iRep, iBlock)       = 1;     % An error is made.
                end
            case 2                                              % No Go Stimulus
                Results.ReactionTime(iRep, iBlock)     = NaN;   % Normally reaction time is non existent
                if any(Button)                                  % If a reaction is given:
                    Results.Errors(iRep, iBlock)       = 1;     % An error is made.
                else
%                     outlet.push_sample(Marker_NoClick);         % If no reaction is given, show correct No Click.
                end
        end

    end
    
    audio1 = PsychPortAudio('GetAudioData',CAP_cue1Hz);
    stoprecord = GetSecs;
    
    audio3 = PsychPortAudio('GetAudioData',CAP_cue3Hz);
    
    % Stop Cues
    switch BlockType
        case 1
            PsychPortAudio('Stop', file(1));
%             outlet.push_sample(Marker_EndBlockCue1Hz);
            stopmoment = GetSecs;
        case 2
            PsychPortAudio('Stop', file(2));
%             outlet.push_sample(Marker_EndBlockCue3Hz);
        case 3
%             outlet.push_sample(Marker_EndBlockNoCue);
        case 4
%             outlet.push_sample(Marker_EndBlockGoOnly);
    end
    
    Line18 = ['End of block ' num2str(iBlock) '\n \n'];
    
    % Pause
    if iBlock == 5 || iBlock == 10 || iBlock == 15
        DrawFormattedText(Window,[Line18 Line13],'Center','Center',white);
        Screen('Flip',Window);
        
        WaitSecs(30);
        
        DrawFormattedText(Window,[Line18 Line13],'Center','Center',white);
        DrawFormattedText(Window,'>>',Xpixel*.9,Ypixel*.9,white);
        Screen('Flip',Window);
        GetClicks(Window);
        
    elseif iBlock == 19
        
    else
        DrawFormattedText(Window,[Line18 Line13],'Center','Center',white);
        Screen('Flip',Window);
        
        WaitSecs(5);
        
        DrawFormattedText(Window,[Line18 Line13],'Center','Center',white);
        DrawFormattedText(Window,'>>',Xpixel*.9,Ypixel*.9,white);
        Screen('Flip',Window);
        GetClicks(Window);
    end
end

    
oldbias1 = PsychPortAudio('LatencyBias', PPA_cue1Hz);
oldbias3 = PsychPortAudio('LatencyBias', PPA_cue3Hz);
oldbias1rec = PsychPortAudio('LatencyBias', CAP_cue1Hz);
oldbias3rec = PsychPortAudio('LatencyBias', CAP_cue3Hz);

% outlet.push_sample(Marker_EndTest);

% End Screen
Line1 = 'End of the experiment \n \n';
Line2 = 'Thank you very much! \n \n';
Line3 = 'Press the left mouse button to exit';

DrawFormattedText(Window,[Line1 Line2 Line3],'Center','Center',white);
Screen('Flip',Window);
GetClicks(Window);
sca;

% cd(Dir);
% save(Filename, '-struct','Results');
 

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

% To Set Sequence of Go and No Go stimuli and Sequence of Cued blocks
function [TestNr, Sequence] = Sequence
% function that randomly determines the sequence of block presentation and the
% stimuli in those blocks.

Samples = 100;
perc = 30;
N = perc/100*Samples;

Exp = ones(Samples,1);
Exp(1:N) = 2;

TrueExp = zeros(21,19);

TrueExp(1,:) = [ones(1,5) 2.*ones(1,5) 3.*ones(1,5) 4.*ones(1,4)];
TrueExp(2:end,1:5) = reshape(Exp(randperm(numel(Exp))),20,5);
TrueExp(2:end,6:10) = reshape(Exp(randperm(numel(Exp))),20,5);
TrueExp(2:end,11:15) = reshape(Exp(randperm(numel(Exp))),20,5);
TrueExp(2:end,16:19) = ones(20,4);

ShuffledTrueExp       = TrueExp(:,randperm(19));
TestNr                = ShuffledTrueExp(1,:);
Sequence              = ShuffledTrueExp(2:end,:);

for icol = 1:19
    a = Sequence(:, icol);                  % Select corresponding column
    v = a == 1;                             % Check for Go Stimuli
    v = v.*(1:length(v))';                  % Get position of said Go Stimuli
    v = v(v~=0);                            % Delete No Go Stimuli
    v(v==1) = [];                           % Remove first two rows
    v(v==2) = [];
    v = v(randperm(numel(v)));              % Randomise the positions
    x = v(1);                               % Select two random positions
    y = v(2);
    a([1 2 x y], :) = a([x y 1 2], :);      % Swap first two rows with random Go Stimuli
    Sequence(:, icol) = a;                  % Replace column
end
end