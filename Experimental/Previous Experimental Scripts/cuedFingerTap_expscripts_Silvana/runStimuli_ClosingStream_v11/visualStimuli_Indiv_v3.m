% Clear the workspace
% close all;
% clearvars;
% sca;
function [markers, timeMarkers]=visualStimuli_Indiv_v3(freq,timeTot,window,ifi,screenYpixels,outlet,black)

%This section is a modification of the code developed by Janne Heijs
% load lab streaming layer library
% lib = lsl_loadlib();
% make a new stream outlet

% info = lsl_streaminfo([lib handle],[name],[type],[channelcount],[fs],[channelformat],[sourceid])

%   > name = name of stream; describes device/product

%   > type = content type of stream (EEG, Markers)

%   > channelcount = nr of channels per sample

%   > fs = samplking rate (Hz) as advertized by data source

%   > channelformat = cf_float32, cf__double64, cf_string, cf_int32, cf_int16

%   > sourceid = unique identifier for source or device, if available
% fs=0;
% id='sdfwerr32432';
% info = lsl_streaminfo(lib,'EOEC','Markers',1,fs,'cf_int32',id);
% outlet = lsl_outlet(info); % thing you push your data through
% info2 = lsl_streaminfo(lib,'TimeVisualIso','Markers',1,fs,'cf_int32',id);
% outlet2 = lsl_outlet(info2); % thing you push your data through
%% Definition markers
V1IIM=1101; %This is the auditory isorhythmic 1 Hz initial marker
V1IFM=1109;%This is the auditory isorhythmic 1 Hz final marker
V3IIM=1103;%This is the auditory isorhythmic 3.2 Hz initial marker
V3IFM=1111;%This is the auditory isorhythmic 3.2 Hz final marker
%%
% Screen('Preference', 'SkipSyncTests', 1)
% % Here we call some default settings for setting up Psychtoolbox
% % Screen('Preference', 'SkipSyncTests', 1)
% PsychDefaultSetup(2);
% 
% % Get the screen numbers
% screens = Screen('Screens');
% 
% % Draw to the external screen if avaliable
% screenNumber = max(screens);
% 
% % Define black and white
% % white = WhiteIndex(screenNumber);
% black = BlackIndex(screenNumber);
% % grey = white / 2;
% % inc = white - grey;
% 
% % Open an on screen window
% [window, windowRect] = PsychImaging('OpenWindow', screenNumber, black);
% 
% % Get the size of the on screen window
% [screenXpixels, screenYpixels] = Screen('WindowSize', window);
% 
% % Query the frame duration
% ifi = Screen('GetFlipInterval', window);
% 
% % Get the centre coordinate of the window
% [xCenter, yCenter] = RectCenter(windowRect);
% 
% % Set up alpha-blending for smooth (anti-aliased) lines
% Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Here we load in an image from file. This one is a image of rabbits that
% is included with PTB
%%
%Two vectors to save the markers and the time of markers
markers=[];
timeMarkers=[];
if freq==3.2
    markerstart=V3IIM;
    startTask=1014;
    endTask=1015;

else
    markerstart=V1IIM;
    startTask=1010;
    endTask=1011;
    
end 
theImageLocation = 'WhiteCircle2.png';
theImage = imread(theImageLocation);

% Get the size of the image
[s1, s2, s3] = size(theImage);

% Here we check if the image is too big to fit on the screen and abort if
% it is. See ImageRescaleDemo to see how to rescale an image.
if s1 > screenYpixels || s2 > screenYpixels
    disp('ERROR! Image is too big to fit on the screen');
    sca;
    return;
end

% Make the image into a texture
imageTexture = Screen('MakeTexture', window, theImage);

% Our will fade in and out with a sine wave function
% See: http://en.wikipedia.org/wiki/Sine_wave
% amplitude = 0.5;
frequency = freq;
angFreq = 2 * pi*frequency;
% startPhase = 0;
time = 0;
if floor(frequency)==frequency
    numImages=timeTot*frequency;
%     disp('Es divisible por freq') 
    
else
    numImages=timeTot*frequency-1;
%     disp('No es divisible')
end 
l=0;
thisContrastSave=1;
% Presentation loop (press any key to exit)

while l<numImages
    %~KbCheck

    % Position of the square on this frame
    thisContrast =  square(angFreq*time);
    
%     plot(thisContrast)
%     disp(thisContrast);
    
%     thisContrastSave;
   if l==0
       outlet.push_sample(startTask)
   end 
    if thisContrast==1 && thisContrastSave==-1
        l=l+1;
       
       outlet.push_sample(markerstart)
    elseif thisContrast==-1 && thisContrastSave==1

    end
    % Draw the image to the screen
    Screen('DrawTexture', window, imageTexture, [], [], 0, [], thisContrast);

    % Increment the time
   
    time = time +1/frequency+ifi;
%     disp(time)
    
    thisContrastSave=thisContrast;
    
    % Flip to the screen
    Screen('Flip', window);
    if KbCheck
        sca;
        return
        
    end
end
outlet.push_sample(endTask)
DrawFormattedText(window, '', 'center', 'center', black);
end
% outlet.delete()
% outlet2.delete()
