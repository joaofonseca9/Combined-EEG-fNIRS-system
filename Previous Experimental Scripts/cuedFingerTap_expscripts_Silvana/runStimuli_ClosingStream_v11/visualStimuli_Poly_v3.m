function [markers, timeMarkers]=visualStimuli_Poly_v3(freq,timeTot,window,ifi,xCenter,yCenter,screenXpixels, screenYpixels,outlet,black)
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
% fs=0; %This value was selected as 4 times the frequency the highest beat.
% %In this case: secondary beat 3/2*main and if main is 3.2 this will be 4.8
% %This multiply by 4 is 19.2. So we use 20.
% id='sdfwerr32432';
% info = lsl_streaminfo(lib,'EOEC','Markers',1,fs,'cf_int32',id);
% outlet = lsl_outlet(info); % thing you push your data through
% info2 = lsl_streaminfo(lib,'TimeVisualPoly','Markers',1,fs,'cf_int32',id);
% outlet2 = lsl_outlet(info2); % thing you push your data through
%% Definition markers
V1PIM=1105; %This is the auditory isorhythmic 1 Hz initial marker
V1PFM=1113;%This is the auditory isorhythmic 1 Hz final marker
V3PIM=1107;%This is the auditory isorhythmic 3.2 Hz initial marker
V3PFM=1115;%This is the auditory isorhythmic 3.2 Hz final marker
V1PIS=1125; %This is the auditory isorhythmic 1 Hz initial marker
V1PFS=1129;%This is the auditory isorhythmic 1 Hz final marker
V3PIS=1123;%This is the auditory isorhythmic 3.2 Hz initial marker
V3PFS=1131;%This is the auditory isorhythmic 3.2 Hz final marker
% screen_size = get(0, 'ScreenSize');
% xCent = screen_size(3)/2;
% yCent = screen_size(4)/2;
% screen_size
% xCent
% yCent

% Set up alpha-blending for smooth (anti-aliased) lines
% Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% We generate two vectors to save the markers and the times
markers=[];
timeMarkers=[];
if freq==3.2
    markerstart=V3PIM;
    startTask=1016;
    endTask=1017;

else
    markerstart=V1PIM;
    startTask=1012;
    endTask=1013;
    
end
theImageLocation = 'WhiteCircle2.png';
theImage = imread(theImageLocation);


position1=[xCenter-400,yCenter-100,xCenter-200,yCenter+100];
position2=[xCenter+200,yCenter-100,xCenter+400,yCenter+100];

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

% 
% imageTexture2=Screen('MakeTexture', window, theImage);

% Our will fade in and out with a sine wave function
% See: http://en.wikipedia.org/wiki/Sine_wave
% amplitude = 0.5;
frequency = freq;
angFreq = 2 * pi*frequency;
angFreq2 = angFreq*2;
% startPhase = 0;
time = 0;
time2=0;
if floor(frequency)==frequency
    numImages=timeTot*frequency;
%     disp('Es divisible por freq') 
    
else
    numImages=timeTot*frequency-1;
%     disp('No es divisible')
end 
l=0;
thisContrastSave=1;
thisContrastSave2=1;
timeIterator=0;
y=0;
time2vec=0;
% Presentation loop (press any key to exit)

%create red dot
Screen('DrawDots', window, [xCenter yCenter], 20, [1 0 0], [], 2);
showTime = 0.01;%the time each dot has to be shown
timeStep = (1/freq)/2;
secondsynch=1;
while l<numImages
    
    %~KbCheck

    % Position of the square on this frame
    thisContrast =  square(angFreq*time);
    thisContrast2 = square(angFreq2*time);
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
        if secondsynch == 1
            secondsynch = 0;
        else
            secondsynch = 1;
        end
        
    end
    if  secondsynch == 1 && thisContrast==1
        thisContrast2 = -1;
    end
    
    
      Screen('DrawTexture', window, imageTexture, [], position1, 0, [], thisContrast);
      Screen('DrawTexture', window, imageTexture, [], position2, 0, [], thisContrast2);
      Screen('DrawDots', window, [xCenter yCenter], 20, [1 0 0], [], 2);
%     
    
    % Increment the time
      time = time + ifi;
      time2vec=[time2vec time];

      thisContrastSave=thisContrast;
      thisContrastSave2=thisContrast2;
     
    


%     
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
