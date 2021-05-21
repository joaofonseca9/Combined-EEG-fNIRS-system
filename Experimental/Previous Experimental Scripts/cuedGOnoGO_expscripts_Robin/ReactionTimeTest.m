%% ReactionTimeTest

clear all; close all; clc;

Subject     = '05';

Dir         = ['C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\EEG scripts (introductory)'];
Filename    = ['s' Subject '_GoCue.mat'];

%% Lab Streaming Layer

lib     = lsl_loadlib();                      % Load Lab Streaming Layer Library
info    = lsl_streaminfo(lib, 'ReactionTimeTest', 'Markers', 1, 0.0, 'cf_int32', 'ReactionTime'); 
outlet  = lsl_outlet(info);                   % Open outlet

% % % Marker Setup % % %
% Test Related
Marker_StartTest              = 1000;         % Start of Experiment
Marker_EndTest                = 1001;         % End of Experiment

% Blocks Related
Marker_StartBlockCue1Hz       = 1100;         % Start of Block with Cues 1Hz
Marker_StartBlockCue3Hz       = 1200;         % Start of Block with Cues 3Hz
Marker_StartBlockNoCue        = 1300;         % Start of Block with no Cues

Marker_EndBlockCue1Hz         = 1101;         % End of Block with Cues 1Hz
Marker_EndBlockCue3Hz         = 1201;         % End of Block with Cues 3Hz
Marker_EndBlockNoCue          = 1301;         % End of Block with no Cues

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
grey        = white/2;

%% Prepare Screen Window

% Open A Window (With Specific Color)
[Window, WindowRect] = PsychImaging('OpenWindow', ScreenNr, black);
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

% Wait momentarily before showing No Go Stimulus
WaitSecs(.5);

% No Go Stimulus
ColorRect   = [1 0 0];              % Red Rectangle  
DotColor    = black;                % Black Dot

% Show Texture
Screen('FillRect', Window, ColorRect, Centered);
% Screen('DrawDots', Window, [Xcenter, Ycenter], DotSize, DotColor, [], 2);
Screen('Flip', Window);

% Save Texture as NoGoStimulus in Texture(1)
NoGoStimulus = Screen('GetImage', Window);
Texture(2)   = Screen('MakeTexture', Window, NoGoStimulus);

WaitSecs(.5);

% Reset Dot Color to Red
DotColor    = [1 0 0];

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

%% Test Settings

% Set Amount of Repetitions
Blocks            = 5;
Repetitions       = 20;

% Set Sequence of Stimuli and Cued Blocks
[TestNr, StimSeq] = Sequence(Blocks, Repetitions);
Results.TestNr    = TestNr;
Results.StimSeq   = StimSeq;
Results.StimShow  = zeros(size(StimSeq));

% Reset Amount of Errors
Results.Errors    = zeros(size(StimSeq));

% Settings to wait for mouse click
NumSecs           = 1.5;
NumFrames         = round(NumSecs / ifi);
WaitFrames        = 1;

%% Experiment text

[WLine1, WLine2, ILine1, ILine2, ILine3, ILine4, ILine5, ILine6, ILine7, ILine8, ILine9, ...
    ILine10, ILine11, ILine12, TLine1, TLine2, TLine4, TLine5, TLine6, ELine1, ELine2, ELine3, Arrow]...
    = GenerateExpText(Blocks, Repetitions);

%% Experiment

WaitSecs(2);
outlet.push_sample(Marker_TestTest);

% Welcome Screen
Screen('TextSize',Window,50);
WelcomeScreen(Window, Xpixel, Ypixel, white, WLine1, WLine2, Arrow);

% Introduction
IntroductionScreen(Window, Xpixel, Ypixel, white, ILine1, ILine2, ILine3, ILine4, ILine5, ILine6, ...
    ILine7, ILine8, ILine9, ILine10, ILine11, ILine12, Arrow, file, outlet, Marker_TestTest);

% Actual Experiment
outlet.push_sample(Marker_StartTest);
for iBlock = 1:size(Results.StimSeq,2)
    
    % % % Cueing % % %
    % Check if Cueing is needed for this block and if so, start cueing.
    Cueing = Results.TestNr(iBlock);
    switch Cueing
        case 1 % Start Cueing with 1 Hz.
            PsychPortAudio('Start', file(1), 1, [], []);
            outlet.push_sample(Marker_StartBlockCue1Hz);
        case 2 % Start Cueing with 3 Hz.
            PsychPortAudio('Start', file(2), 1, [], []);
            outlet.push_sample(Marker_StartBlockCue3Hz);
        case 3
            outlet.push_sample(Marker_StartBlockNoCue);
    end
    
    % Set the baseline to even out electrodes.
    Baseline(Window, white, TLine1, TLine2, TLine6, Xcenter, Ycenter, DotSize, DotColor);

    % % % Start Collection of Samples % % %
    for iRep = 1:size(Results.StimSeq,1)
        clear Button

        % % % Stimulus presentation % % %
        % Select Stimulus from provided sequence.
        Stimulus = Results.StimSeq(iRep, iBlock);

        % Present Stimulus after 2 to 4 seconds.
        WaitSecs(randi([2, 4]));                            % Randomly wait 2 to 4 seconds.
        Screen('DrawTexture', Window, Texture(Stimulus));   % Select Texture corresponding to stimulus.
        Results.StimShow(iRep, iBlock) = GetSecs;   
        Screen('Flip', Window);                                                    % Present the stimulus.
        
        % Push marker corresponding with stimulus and start counting for
        % reaction times
        switch Stimulus
            case 1
                Onset = GetSecs;
                outlet.push_sample(Marker_GoStimulus);
            case 2
                Onset = GetSecs;
                outlet.push_sample(Marker_NoGoStimulus);
        end
        
        % % % Reaction % % %
        % Wait until NumSecs have passed or until left mouse button is
        % clicked. Don't forget to correct for one more WaitSecs in analysis.
        for frame = 1:NumFrames                             % Select how much frames the for loop takes.
            [~, ~, Button] = GetMouse(Window);
            if Button(1)                                    % If left mouse button is pressed. 
                Offset = GetSecs;
                outlet.push_sample(Marker_Click);                                            % Push the marker for the reaction time offset
                break                                       % and break out of for loop.
            end
            WaitSecs(WaitFrames * ifi);                     % Select how long one frame takes at most.
        end

        % Change Screen to only focus dot.
        Screen('DrawDots', Window, [Xcenter, Ycenter], DotSize, DotColor, [], 2);
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
                    outlet.push_sample(Marker_ErrorNoClick);    % If not so, Error No Click is presented
                    Results.ReactionTime(iRep, iBlock) = NaN;   % Reaction time is non existent
                    Results.Errors(iRep, iBlock)       = 1;     % An error is made.
                end
            case 2                                              % No Go Stimulus
                Results.ReactionTime(iRep, iBlock)     = NaN;   % Normally reaction time is non existent
                if any(Button)                                  % If a reaction is given:
                    Results.Errors(iRep, iBlock)       = 1;     % An error is made.
                else
                    outlet.push_sample(Marker_NoClick);         % If no reaction is given, show correct No Click.
                end
        end

    end

    % % % Stop Cueing % % %
    switch Cueing
        case 1
            PsychPortAudio('Stop', file(1));
            outlet.push_sample(Marker_EndBlockCue1Hz);
        case 2
            PsychPortAudio('Stop', file(2));
            outlet.push_sample(Marker_EndBlockCue3Hz);
        case 3
            outlet.push_sample(Marker_EndBlockNoCue);
    end

    % % % Pause % % %
    if iBlock < Blocks*3                                    % Check if all blocks have been run through.
        PauseScreen(Window, white, TLine4, TLine5, Xcenter, Ycenter, DotSize, DotColor, iBlock);
    end                                                     % If not, show Pause Screen.
    
end
outlet.push_sample(Marker_EndTest);

% End Screen
EndScreen(Window, white, ELine1, ELine2, ELine3)

sca;

cd(Dir);
save(Filename, '-struct','Results');

%% Functions Used

% Text Generation
function [WLine1, WLine2, ILine1, ILine2, ILine3, ILine4, ILine5, ILine6, ILine7, ILine8, ILine9, ...
    ILine10, ILine11, ILine12, TLine1, TLine2, TLine4, TLine5, TLine6, ELine1, ELine2, ELine3, Arrow]...
    = GenerateExpText(Blocks, Repetitions)
% This function generates the text seen on the screen during the
% experiment.

    % Welcome Screen
    WLine1          = 'Welcome to this experiment';
    WLine2          = 'Press the left mouse button for a short introduction';

    % Introduction
    ILine1          = 'Introduction';
    ILine2          = 'During the experiment, you should do only two things \n \n';
    ILine3          = 'The next screens will show you what to do';
    ILine4          = 'Press the left mouse button';
    ILine5          = 'Do not press the left mouse button';
    ILine6          = ['You will receive ' num2str(Repetitions) ' images per block \n \n'];
    ILine7          = ['After each of the ' num2str(Blocks*3) ' blocks you can take a short break \n \n'];
    ILine8          = 'If everything is clear, press the left mouse button to start';
    ILine9          = 'During the experiment, you will encounter 3 types of sounds \n \n';
    ILine10         = '1. You could experience silence';
    ILine11         = '2. You could experience this sound';
    ILine12         = '3. You could experience this sound';

    % Actual Experiment
    TLine1          = 'Sit still and do not speak \n \n';
    TLine2          = 'Focus on the dot \n \n \n \n \n';
%     TLine3          = 'End of block X \n \n';
    TLine4          = 'Take a short break \n \n';
    TLine5          = 'Press the left mouse button to start again';
    TLine6          = 'The next block will start in 10 seconds';

    % End Screen
    ELine1          = 'End of the experiment \n \n';
    ELine2          = 'Thank you very much! \n \n';
    ELine3          = 'Press the left mouse button to exit';

    % Arrow
    Arrow           = '>>';
end

% Screen Types
function WelcomeScreen(Window, X,Y, Color, WLine1, WLine2, Arrow)
% The Welcome Screen of the experiment. It will show a short welcome and
% that an introduction will be given. It keeps waiting on user input (mouse
% click) to continue.

    DrawFormattedText(Window,WLine1,'Center','Center',Color);
    DrawFormattedText(Window,Arrow,X*.9,Y*.9,Color);
    Screen('Flip',Window);
    GetClicks(Window);

    DrawFormattedText(Window,WLine2,'Center','Center',Color);
    DrawFormattedText(Window,Arrow,X*.9,Y*.9,Color);
    Screen('Flip',Window);
    GetClicks(Window);
end

function IntroductionScreen(Window, X, Y, Color, ILine1, ILine2, ILine3, ILine4, ILine5, ILine6, ...
    ILine7, ILine8, ILine9, ILine10, ILine11, ILine12, Arrow, Soundfile, outlet, Marker_TestTest)
% The introduction, it shows how to perform the experiment. That a click
% must be provided on a go stimulus and that the user must not give a
% reaction on a no go stimulus. 

    % Screen 1
    DrawFormattedText(Window,ILine1,X*.1,Y*.1,Color);
    DrawFormattedText(Window,[ILine2 ILine3],'Center','Center',Color);
    DrawFormattedText(Window,Arrow,X*.9,Y*.9,Color);
    Screen('Flip',Window);
    GetClicks(Window);

    % Screen 2 (Example Go Stimulus)
    DrawFormattedText(Window,ILine1,X*.1,Y*.1,Color);
    DrawFormattedText(Window,ILine4,'Center',Y*.75,Color);
    DrawFormattedText(Window,Arrow,X*.9,Y*.9,Color);
    
    BaseRect = [0 0 .4*X 0.4*Y];
    Centered = CenterRectOnPointd(BaseRect, X/2, Y/2);
    Screen('FillRect', Window, [0 1 0], Centered);
    
    Screen('Flip',Window);    GetClicks(Window);

    % Screen 3 (Example No Go Stimulus)
    DrawFormattedText(Window,ILine1,X*.1,Y*.1,Color);
    DrawFormattedText(Window,ILine5,'Center',Y*.75,Color);
    DrawFormattedText(Window,Arrow,X*.9,Y*.9,Color);
    
    Screen('FillRect', Window, [1 0 0], Centered);
    
    Screen('Flip',Window);
    GetClicks(Window);
    
    % Screen 4 (Examples of sound Silence)
    DrawFormattedText(Window,ILine1,X*.1,Y*.1,Color);
    DrawFormattedText(Window,[ILine9 ILine10],'Center','Center',Color);
    DrawFormattedText(Window,Arrow,X*.9,Y*.9,Color);
    
    Screen('Flip',Window);
    GetClicks(Window);
    
    % Screen 5 (Example of sound 1 Hz)
    DrawFormattedText(Window,ILine1,X*.1,Y*.1,Color);
    DrawFormattedText(Window,ILine11,'Center','Center',Color);
    DrawFormattedText(Window,Arrow,X*.9,Y*.9,Color);
        
    Screen('Flip',Window);
    
    PsychPortAudio('Start', Soundfile(1), 1, [], []);
    WaitSecs(5);
    PsychPortAudio('Stop',Soundfile(1));
    GetClicks(Window);
    
    outlet.push_sample(Marker_TestTest);
    
    % Screen 6 (Example of sound 3 Hz)
    DrawFormattedText(Window,ILine1,X*.1,Y*.1,Color);
    DrawFormattedText(Window,ILine12,'Center','Center',Color);
    DrawFormattedText(Window,Arrow,X*.9,Y*.9,Color);
    
    Screen('Flip',Window);
    
    PsychPortAudio('Start', Soundfile(2), 1, [], []);
    WaitSecs(3);
    PsychPortAudio('Stop',Soundfile(2));
    GetClicks(Window);

    % Screen 7 (Amount Of Samples, Questions)
    DrawFormattedText(Window,ILine1,X*.1,Y*.1,Color);
    DrawFormattedText(Window,[ILine6 ILine7 ILine8],'Center','Center',Color);
    
    Screen('Flip',Window);
    GetClicks(Window);
end

function PauseScreen(Window, Color, TLine4, TLine5, Xcenter, Ycenter, DotSize, DotColor, iBlock)
% This function shows the pause screen that is given when not all the
% blocks are completed. 

    TLine3 = ['End of block ' num2str(iBlock) ' \n \n'];
    
    DrawFormattedText(Window, [TLine3, TLine4], 'Center', 'Center', Color);
    Screen('Flip', Window);
    GetClicks(Window);
    
    WaitSecs(3);

    DrawFormattedText(Window, [TLine3, TLine5], 'Center', 'Center', Color);
    Screen('Flip', Window);
    GetClicks(Window);

    Screen('DrawDots', Window, [Xcenter, Ycenter], DotSize, DotColor, [], 2);
    Screen('Flip', Window);
end

function EndScreen(Window, white, ELine1, ELine2, ELine3)
% This screen is shown when the experiment is completed.

    DrawFormattedText(Window,[ELine1 ELine2 ELine3],'Center','Center',white);
    Screen('Flip',Window);
    GetClicks(Window);
end

function Baseline(Window, Color, TLine1, TLine2, TLine6, Xcenter, Ycenter, DotSize, DotColor)
% This function shows the user to sit still and focus on the dot in the
% middle of the screen.

    DrawFormattedText(Window,[TLine1 TLine2 TLine6],'Center','Center',Color);
    Screen('DrawDots', Window, [Xcenter, Ycenter], DotSize, DotColor, [], 2);
    Screen('Flip',Window);
    WaitSecs(6);
    
    Screen('DrawDots', Window, [Xcenter, Ycenter], DotSize, DotColor, [], 2);
    Screen('Flip',Window);
    WaitSecs(4);
end

% To Set Sequence of Go and No Go stimuli and Sequence of Cued blocks
function [TestNr, Sequence] = Sequence(Blocks,Reps)
% This function creates the sequences for both the Blocks and the
% Repetitions in those blocks. Every experiment has a set amount of Go and
% No Go cues. Between experiments the amounts stay the same. 

Samples               = Blocks*Reps;
perc                  = 30;
N                     = perc/100*Samples;

Exp                   = ones(Samples,1);
Exp(1:N)              = 2;

TrueExp               = zeros(Reps+1,Blocks*3);
TrueExp(1,:)          = [1.*ones(1,Blocks) 2.*ones(1,Blocks) 3.*ones(1,Blocks)];

TrueExp(2:end,1:Blocks)             = reshape(Exp(randperm(numel(Exp))), Reps, Blocks);
TrueExp(2:end,Blocks+1:Blocks*2)    = reshape(Exp(randperm(numel(Exp))), Reps, Blocks);
TrueExp(2:end,Blocks*2+1:Blocks*3)  = reshape(Exp(randperm(numel(Exp))), Reps, Blocks);

ShuffledTrueExp       = TrueExp(:,randperm(Blocks*3));
TestNr                = ShuffledTrueExp(1,:);
Sequence              = ShuffledTrueExp(2:end,:);

% Check if Reps are enough to have at least 2 different Go Stimuli to
% switch out with [Can be removed, however you can receive No Go stimuli 
% in the first two reps then]
if Reps >= 10
    % Replace first two samples of each block with Go stimuli
    % Without losing any No Go Stimuli
    for icol = 1:Blocks*3
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
   
%% Further Notes
