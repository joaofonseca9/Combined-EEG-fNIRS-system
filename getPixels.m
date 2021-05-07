%% SCREEN PREPARATION
Screen('Preference', 'SkipSyncTests', 1);
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

%% Draw crosses and letters
Screen('TextSize', window, 36);
Screen('DrawLines', window, allCoords,...
    lineWidthPix, white, [xCenter yCenter], 2);
Screen('Flip', window);
cross=Screen('GetImage',window);

Screen('TextSize', window, 100);
DrawFormattedText(window, 'G','center','center', white);
Screen('Flip', window);
G=Screen('GetImage',window);

DrawFormattedText(window, 'Q','center','center', white);
Screen('Flip', window);
Q=Screen('GetImage',window);

DrawFormattedText(window, 'O','center','center', white);
Screen('Flip', window);
O=Screen('GetImage',window);

DrawFormattedText(window, 'C','center','center', white);
Screen('Flip', window);
C=Screen('GetImage',window);

DrawFormattedText(window, 'D','center','center', white);
Screen('Flip', window);
D=Screen('GetImage',window);

sca;
%% Number of pixels
pixelsG=sum(sum(G(:,:,1)))
pixelsO=sum(sum(O(:,:,1)))
pixelsQ=sum(sum(Q(:,:,1)))
pixelsC=sum(sum(C(:,:,1)))
pixelsD=sum(sum(D(:,:,1)))
pixelsCross=sum(sum(cross(:,:,1)))

imwrite(G,'G.png')
imwrite(O,'O.png')
imwrite(Q,'Q.png')
imwrite(C,'C.png')
imwrite(D,'D.png')
imwrite(cross,'cross.png')

