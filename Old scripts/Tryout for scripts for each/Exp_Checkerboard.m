%% CHECKERBOARD TASK

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

%% GENERATE CHECKERBOARD 

XYChecks = [24 18];                                             % nr of checks 4:3
nrChecks = XYChecks(1)*XYChecks(2);                             % total number of checks
dim = screenYpixels/XYChecks(2);                                      % dimension of a check [nr of pixels]
if dim*XYChecks(1)>screenXpixels
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
flipNr = 300;                          % total nr of flips made (300)
flipNrCount = 0;                            % counter for nr of flips made
flipFrameCount = 0;                         % counting the number of frames                         
flipFrameWait = 0;                          % time to wait in frames for a flip

textureCue = [1 2];                         % cue that determines which texture (checkerboard) will be shown
SaveFrameLog = [];                          % save Screen presentation > (start=0; checkerboard = 1/2)

%% >> EXPERIMENT << %%
%%%%%%%%%%%%%%%%%%%%%%
%% Instructions
Screen('TextSize',window,36);
DrawFormattedText(window, 'The only thing you need to do for this task is look at the screen \n \n A red dot will indicate where you should look. \n \n Press any key to continue.','center','center', white);
vbl = Screen('Flip', window); 
KbStrokeWait



%% %>> START SCREEN <<%%%
% > grey screen with focus dot with text: 'Focus on red dot; press any key to START'
% > after pressing any key, the text removes
Screen('TextSize',window,36);
DrawFormattedText(window, 'Focus on the red dot \n\n\n\n Press any key to START','center','center', white);
vbl = Screen('Flip', window); 
Screen('DrawDots', window, [xCenter yCenter], dotSizePix, dotColor, [], 2);
KbStrokeWait; 
outlet.push_sample(Marker_start); % wait for key press and send LSL-marker (start)

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
    Screen('TextSize',window,36); % Stop screen
    DrawFormattedText(window, 'End of experiment. \n \n Thank you for participating! \n\n Press any key to EXIT','center','center', white);
    vbl = Screen('Flip', window); 
    outlet.push_sample(Marker_stop); 
end
KbStrokeWait; % wait for keypress
sca; % close screen

close all;
clear all;