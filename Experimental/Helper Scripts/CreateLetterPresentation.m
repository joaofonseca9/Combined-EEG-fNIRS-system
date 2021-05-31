%% CREATE MOVIE FILE FOR DUAL TASKING
%Please define filepath before
%It will save a .mov file of the presentation as well as a .txt with the
%number
filepath='C:\Users\joaop\Desktop\LetterPresentation';

%Number of trials per block
N_trials=20;
%Number of letter counting trials = 2 * N_trials (2 blocks of dual tasking)
N=N_trials*2; 


%Create N letter presentations
for ii=1:N
    filename=['LetterPresentation_', num2str(ii,'%d'),'.mov'];
    [letters,isLetterFrame]=letterPresentation (filepath,filename);
    
    % Save letter order
    str=['LetterPresentation_', num2str(ii,'%d')];
    save(str,'letters');
    str=['LetterPresentation_isLetterFrame_', num2str(ii,'%d')];
    save(str,'isLetterFrame');
end
function [letters,isLetterFrame]=letterPresentation (filepath,filename)
cd(filepath);
%Synch test skip => comment when actually testing patient
Screen('Preference', 'SkipSyncTests', 1);

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
fixCrossDimPix = 50; % Here we set the size of the arms of our fixation cross
xCoords = [-fixCrossDimPix fixCrossDimPix 0 0]; % Set the coordinates (these are all relative to zero we will let the drawing routine center the cross in the center of our monitor for us)
yCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
allCoords = [xCoords; yCoords];
lineWidthPix = 10;% Set the line width for the fixation cross


%% Presentation of the letters on the screen (dual task). -> is random.
%Participant has to count the times that G was presented.
Letterlist='DGOQ';
letter_order=randi(length(Letterlist), 1, N_letters);
letters={Letterlist(letter_order)};



%% LETTER PRESENTATION
% onset=zeros(1,N_letters);
framerate=30;
moviePtr = Screen('CreateMovie', window, filename,[], [], framerate);

%We want a frequency of 1.25 Hz for the letter presentation.
%‘frameDuration’ which defaults
% to one. The parameter defines the display duration of that frame as the fraction
% ‘frameDuration’ / ‘frameRate’ seconds, so ‘frameRate’ defines the denominator of
% that term
%So we have to define how many frames the letters and the red cross will
%occupy. so frame_letter= round(time_letter * framerate)

frame_letter=zeros(1,N_letters);
frame_cross=round(0.25*framerate)*ones(1,N_letters);
% frameduration= time_letter/framerate;
for n=1:N_letters
    %Defining number of frames per letter and cross
    frame_letter(n)=round((rand(1)+0.55)*framerate); %Speed with which the letters are presented
    
    %Present random letter on the screen
    Screen('TextSize', window, 100);
    DrawFormattedText(window, letters{1}(n),'center','center', white); 
    onset(n)=GetSecs;
    vbl = Screen('Flip', window);
    Screen('AddFrameToMovie',window,windowRect,'frontBuffer',moviePtr,frame_letter(n));
    
    %Between each letter show a red fixation cross
    Screen('DrawLines', window, allCoords,...
        lineWidthPix, [255 0 0], [xCenter yCenter], 2);
    Screen('Flip', window);
    Screen('AddFrameToMovie',window,windowRect,'frontBuffer',moviePtr,frame_cross(n));
end
Screen('DrawLines', window, allCoords,...
        lineWidthPix, white, [xCenter yCenter], 2);
Screen('Flip', window);
frame_white_cross=round(((8-5).*rand(1) + 5)*framerate);
Screen('AddFrameToMovie',window,windowRect,'frontBuffer',moviePtr,frame_white_cross);

Screen('FinalizeMovie', moviePtr);
sca;

% letter_timestamps=zeros(1,N_letters);
% for k=2:N_letters
%     letter_timestamps(k)=letter_timestamps(k-1)+frame_cross(k-1)+frame_letter(k-1);
% end

isLetterFrame=zeros(1,sum(frame_letter)+sum(frame_cross)+frame_white_cross);
tmp=1;
for k=1:length(frame_letter)
    isLetterFrame(tmp)=1; %mark only the first frame of the letter to 1
    isLetterFrame(tmp+frame_letter(k):tmp+frame_letter(k)+frame_cross(k)-1)=0;
    tmp=tmp+frame_letter(k)+frame_cross(k);
end

%The last tmp will be the sample of the first white cross frame
isLetterFrame(tmp)=2;

end