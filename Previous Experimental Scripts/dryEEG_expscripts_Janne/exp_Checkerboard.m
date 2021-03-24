












 % EXP CHECKERBOARD

% Clear the workspace and the screen
%sca; close all; clear all; clc

%% SELECT SETTINGS

% screen laptop [cm]
% > laptop Janne [34.4 19.4]
% > laptop Jurjan [31 17.4]
Xlength = 31;           % X length screen [cm]
Ylength = 17.4;         % Y length screen [cm]

% number of trials
nrTrials = 300;         % nr of flips to be made (300 in exp)

%% LSL OUTLET SENDING EVENTS

% load lsl library and make a new stream outlet
lib = lsl_loadlib();            % load labstreaming layer library
info = lsl_streaminfo(lib,'EOEC','Markers',1,0.0,'cf_int32','sdfwerr32432');
outlet = lsl_outlet(info);      % open an outlet

% create marker id's 
Marker_CHECK = 1255;        % checkerboard flip
Marker_start = 1555;        % start signal 
Marker_stop = 1500;         % stop signal 


%% General screen preparation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Call some default settings for Psychtoolbox
PsychDefaultSetup(2);

% Get screen numbers
% > external = max(screens) ; laptop = 1
screens = Screen('Screens');  
screenNumber = 0;
% > set priority of computer to screen
maxPriorityLevel = MaxPriority(screenNumber);
Priority(maxPriorityLevel);

% Start screen of PTB is black (not white)
% Screen('Preference', 'VisualDebugLevel', 1);

% Define white/black/grey
% > returns intensity value to produce white/black at current screen depth
white = WhiteIndex(screenNumber);   
black = BlackIndex(screenNumber);   
grey = white / 2;                   

%% Open/prepare screen window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Open an on screen window
% > get size of the on-screen window (screenX/Ypixels
% > get the centre coordinate of the window (x/yCenter)
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey);
[Xpixels, Ypixels] = Screen('WindowSize', window);
[xCenter, yCenter] = RectCenter(windowRect);

% Enable alpha blending for anti-aliasing
% > Alpha blending is a way to combine color values of pixels already in the window
% > with new color values from drawing commands.
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
ifi = Screen('GetFlipInterval', window);    % Query the frame duration [s]?

%% GENERATE CHECKERBOARD 

XYChecks = [24 18];                                             % nr of checks 4:3
nrChecks = XYChecks(1)*XYChecks(2);                             % total number of checks
dim = Ypixels/XYChecks(2);                                      % dimension of a check [nr of pixels]
if dim*XYChecks(1)>Xpixels
    error('ERROR: dimensions of checks are not correct');
end
baseRect = [0 0 dim*XYChecks(1) dim*XYChecks(2)];               % rect of dim*nrXpixels by dim pixels*nrYpixels
[Xpos,Ypos] = meshgrid(0:1:XYChecks(1)-1, 0:1:XYChecks(2)-1);   % Make coordinates for our grid of squares

% Reshape the matrices of coordinates into a vector (reshape) & 
% Scale grid spacing to the size of our squares and centre (.*dim)
Xgrid = reshape(Xpos, 1, nrChecks) .*dim;    
Ygrid = reshape(Ypos, 1, nrChecks) .*dim;  

% Visual angle
% > visual angle (V) = 2*arctan(size on screen (S) / 2* distance to screen(D))
Vangle = 1;                                                     % visual angle [degree]
dist = 65;                                                      % distance to screen [cm]
Lcheck = (2*dist)*tand(Vangle/2);                               % length of check on screen [cm]

% Set the colors of each of our squares
% > 1= white; 0=black
for iNrChecks = 1:XYChecks(2)
    if mod(iNrChecks,2)==0
        checkColor(iNrChecks,:) = repmat([1 0],1,XYChecks(1)/2); % [1 0 1 0 ...]
    else
        checkColor(iNrChecks,:) = repmat([0 1],1,XYChecks(1)/2); % [0 1 0 1 ...]
    end
end
checkColor = reshape(checkColor, 1, nrChecks);

% pixelnumbers for the checks 
% > [(xpos;ypos;length;heigth) x nrChecks] 
checks = [Xgrid ; Ygrid ; Xgrid+dim ; Ygrid+dim];

% Create Grid for the color with size of 'checks' (GridColor) &
% Overlap the GridColor on the Grid of the Checkerboard (=GridChecks)
for i=1:nrChecks
    GridColor = meshgrid(checks(1,i):1:checks(3,i)-1, checks(2,i):1:checks(4,i)-1) *0 +checkColor(1,i);
    GridChecks(checks(2,i)+1:1:checks(4,i), checks(1,i)+1:1:checks(3,i)) = GridColor;
    % imshow(GridChecks)
end

% Create checkerboard texture (1) and inverse checkerboard (2)
Texture(1)  = Screen('MakeTexture', window, GridChecks);
Texture(2)  = Screen('MakeTexture', window, 1-GridChecks);

% create red focus dot
dotColor = [1 0 0];                 % Set the color of our dot to full red (RBG)
dotSizePix = 20;                    % Dot size in pixels

% Draw the dot to the screen. 
% > Screen('DrawDots', windowPtr, xy [,size] [,color] [,center] [,dot_type]);
% > dot_type: round dots (1,2,3) > 2 tries to use high-quality anti-aliasing
Screen('DrawDots', window, [xCenter yCenter], dotSizePix, dotColor, [], 2);

%% GENERAL VARIABLES

% experimental variables
% > 'start' = time between keypress to start and first checkerboard presentation; 
% > grey screen with focus dot will be presented
startTime = 10;                             % time between pressing start and first checkerboard presentation (=grey screen) [sec]
startFrame = round(startTime / ifi);        % time between pressing start and first checkerboard presentation [nrSamples]
startFrameCount = 0;                        % counter for nr of start frames

% > 'flip' = flipping the contrast of the checkerboard presentations
flipTime = 0.45;                            % duration of checkerboard presentation, before flipping [s]
flipFrame = round(flipTime / ifi);          % duration of checkerboard presentation, before flipping [nrSamples]
flipNr = nrTrials;                          % total nr of flips made (300)
flipNrCount = 0;                            % counter for nr of flips made
flipFrameCount = 0;                         % counting the number of frames                         
flipFrameWait = 0;                          % time to wait in frames for a flip

textureCue = [1 2];                         % cue that determines which texture (checkerboard) will be shown
SaveFrameLog = [];                          % save Screen presentation > (start=0; checkerboard = 1/2)

%% >> EXPERIMENT << %%
%%%%%%%%%%%%%%%%%%%%%%

%%%>> START SCREEN <<%%%
% > grey screen with focus dot with text: 'Focus on red dot; press any key to START'
% > after pressing any key, the text removes
Screen('TextSize',window,70);
DrawFormattedText(window, 'Focus on the red dot \n\n\n\n Press any key to START','center','center', white);
vbl = Screen('Flip', window); 
Screen('DrawDots', window, [xCenter yCenter], dotSizePix, dotColor, [], 2);
KbStrokeWait; outlet.push_sample(Marker_start); % wait for key press and send LSL-marker (start)

% > after pressing any key, the text removes
% > 10sec grey screen with focus dot before checkerboards are presented 
while startFrameCount ~= startFrame
    vb1 = Screen('Flip', window);                   % Sync us to the vertical retrace
    startFrameCount = startFrameCount + 1;          % increment of counter
    SaveFrameLog = [SaveFrameLog; 0];               % save screenpresentation (startstimuli=0) for each frame
end

%%%>> CHECKERBOARD <<%%%
if startFrameCount == startFrameCount
    while flipNrCount ~= flipNr
        flipFrameCount = flipFrameCount + 1;            % increment of counter
        SaveFrameLog = [SaveFrameLog; textureCue(1)];   % save screenpresentation (checkerboard 1/2) for each frame

        % Draw our texture to the screen
        Screen('DrawTexture', window, Texture(textureCue(1)));
        Screen('DrawDots', window, [xCenter yCenter], dotSizePix, dotColor, [], 2);
        vbl = Screen('Flip', window, vbl + flipFrameWait); % ?? (flipFrameWait - 0.5) * ifi) / 0= default; flip on the next possible video retrace
              
        % Reverse texture cue to show inverse checkerboard if the time is up
        if flipFrameCount == flipFrame
            outlet.push_sample(Marker_CHECK);           % send LSL-marker (checkerboard flip)
            textureCue = fliplr(textureCue);            % flip index of texture
            flipNrCount = flipNrCount + 1;              % increment counter
            flipFrameCount = 0;                         % set frame counter on zero
        end
    end
end

%%%>> END SCREEN <<%%%
% > grey screen with focus dot with text: End of experiment; Press any key to EXIT'
if flipNrCount == flipNr
    Screen('DrawDots', window, [xCenter yCenter], dotSizePix, dotColor, [], 2);
    Screen('TextSize',window,70); % Stop screen
    DrawFormattedText(window, 'End of experiment \n\n Press any key to EXIT','center','center', white);
    vbl = Screen('Flip', window); 
    outlet.push_sample(Marker_stop); 
end
KbStrokeWait; % wait for keypress
sca; % close screen

close all;
clear all;


%% NOTES SCRIPT

% info = lsl_streaminfo([lib handle],[name],[type],[channelcount],[fs],[channelformat],[sourceid])
%   > name = name of stream; describes device/product 
%   > type = content type of stream (EEG, Markers)
%   > channelcount = nr of channels per sample
%   > fs = samplking rate (Hz) as advertized by data source 
%   > channelformat = cf_float32, cf__double64, cf_string, cf_int32, cf_int16
%   > sourceid = unique identifier for source or device, if available
