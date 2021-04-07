%% CREATE MOVIE FILE FOR DUAL TASKING

%Before starting the automaticity test, clear the workspace.
clear all

%Synch test skip => comment when actually testing patient
Screen('Preference', 'SkipSyncTests', 1);

%Sequences used in order to be able to print in the command window if
%to generate a new sequence use randi([1 4], 1, 12)
sequencesprint = {('4 3 4 1 4 1 2 4 3 2 1 2'),('2 1 2 3 2 1 3 2 4 2 4 1')};

sequences = {split(sequencesprint(1))',split(sequencesprint(2))'} ;

N_letters=8;
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


%% Presentation of the letters on the screen (dual task). -> is random.
%Participant has to count the times that G was presented.
Letterlist='AGOL';
letter_order=randi(length(Letterlist), 1, N_letters);
value={Letterlist(letter_order)};

%% LETTER PRESENTATION
% onset=zeros(1,N_letters);
moviename = 'LetterPresentation.mov';
moviePtr = Screen('CreateMovie', window, moviename,512, 512, 30, ':CodecSettings=AddAudioTrack=2@48000 Videoquality=0.5 Profile=2');
for n=1:N_letters
    %Present random letter on the screen
    Screen('TextSize', window, 100);
    DrawFormattedText(window, value{1}(n),'center','center', white); 
%     onset(n)=GetSecs;
    vbl = Screen('Flip', window);
    Screen('AddFrameToMovie',window,[],[],moviePtr);
    WaitSecs (1/1.25/2); 

    %Between each letter show a red fixation cross
    Screen('DrawLines', window, allCoords,...
        lineWidthPix, [255 0 0], [xCenter yCenter], 2);
    Screen('Flip', window);
    Screen('AddFrameToMovie',window,[],[],moviePtr);
    WaitSecs (1/1.25/2);
end


Screen('FinalizeMovie', moviePtr);

%% Try to play movie

Screen('TextSize', window, 100);
DrawFormattedText(window, 'Playing movie in 2 seconds','center','center', white); 
%     onset(n)=GetSecs;
vbl = Screen('Flip', window);
WaitSecs (2); 

% Settings to open&play movie
% moviefile=[pwd filesep moviename];
% [moviePtr, movieduration, fps, imgw, imgh, ~, ~, hdrStaticMetaData] = Screen('OpenMovie', window, moviename);
% moviePtr = Screen('OpenMovie', window, [pwd filesep moviefile]);
Screen('PlayMovie', moviePtr, 1);

while ~KbCheck
    % Wait for next movie frame, retrieve texture handle to it
    tex = Screen('GetMovieImage', window, moviePtr);

    % Valid texture returned? A negative value means end of movie reached:
    if tex<=0
        % We're done, break out of loop:
        break;
    end

    % Draw the new texture immediately to screen:
    Screen('DrawTexture', window, tex);

    % Update display:
    Screen('Flip', window);

    % Release texture:
    Screen('Close', tex);
end

% Stop playback:
Screen('PlayMovie', movie, 0);

% Close movie:
Screen('CloseMovie', movie);

sca;

%Check letter showing frequency
% avg_wavelength=0;
% for ii=2:N_letters
%     avg_wavelength=avg_wavelength+onset(ii)-onset(ii-1);
% end
% avg_wavelength=avg_wavelength/(N_letters-1);
% freq=1/avg_wavelength


