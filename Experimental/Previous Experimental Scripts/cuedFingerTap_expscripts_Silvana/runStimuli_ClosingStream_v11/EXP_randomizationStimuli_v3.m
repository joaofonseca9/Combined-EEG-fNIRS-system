%Randomization of the stimuli. We have 8 cases: 

% 1. Visual, 1Hz, isorhythmic
%2. Visual, 1Hz, polyrhythmics
%3. Visual, 3.2Hz, isorhythmic
%4. Visual, 3.2Hz,polyrhythmic

%5. Auditory, 1Hz, isorhythmic
%6. Auditory, 1Hz, polyrhythmic
%7. Auditory, 3.2Hz, isorhythmic
%8. Auditory, 3.2Hz, polyrhythmic

%% This section generates the .wav files with the Beep sounds 
%This needs to be run only if you want to generate again the tones. 
fs=48000; %This is the sample rate (Frequency at which the audio port can sample)
[beep_main,beep_sec]=generationTones(fs); %This function generates the two
%wav signals at different frequencies for the polyrhythm


%% This section runs the different stimuli in a randomized way. 

randomStimuli=8; %Number of stimuli trials. The explanation of each is done above
lib = lsl_loadlib();
fs=0.0; %This value was selected since the sampling frequency of the markers 
%is not constant. 
%%This value was selected as 4 times the frequency the highest beat.
%%In this case: secondary beat 3/2*main and if main is 3.2 this will be 4.8
%%This multiply by 4 is 19.2. So we use 20.
id='sdfwerr32432'; 
% id='SilvanaMarkers';
info = lsl_streaminfo(lib,'EOEC','Markers',1,fs,'cf_int32',id);
outlet = lsl_outlet(info); % thing you push your data through

% id='SilvanaMarkersTime';
% info2 = lsl_streaminfo(lib,'STime','Markers',1,fs,'cf_int32',id);
% outlet2 = lsl_outlet(info2); % thing you push your data through
% %The auditory files (beeps) are loaded
nameSignal1='Beep880Hz.wav'; 
nameSignal2='Beep146Hz.wav';
nameSignal='Beep146Hz.wav';

timeStimuli=[];
typeStimuli=[];
totalStimuli=1:randomStimuli;
stimuli=totalStimuli;
timeTrials=5; %Pause time between each trial (in sconds) of each type of task
numTrials=5;% Number trials per task
timeTask=15; %This is the pause in time between types of task
Screen('Preference', 'SkipSyncTests', 1)
% Here we call some default settings for setting up Psychtoolbox
% Screen('Preference', 'SkipSyncTests', 1)
PsychDefaultSetup(2);
% Get the screen numbers
screens = Screen('Screens');
% Draw to the external screen if avaliable
screenNumber = max(screens);
% Define black and white
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
% Open an on screen window
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, black);
% Set up alpha-blending for smooth (anti-aliased) lines
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Get the size of the on screen window
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect);
% Initiation variables for visual cues
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
position1=[xCenter-400,yCenter-100,xCenter-200,yCenter+100];
position2=[xCenter+200,yCenter-100,xCenter+400,yCenter+100];
ifi = Screen('GetFlipInterval', window);

%Initiation of the presentation with the instructions

% fs=0;
% id='sdfwerr32432';
% info = lsl_streaminfo(lib,'EOEC','Markers',1,fs,'cf_int32',id);
% outlet = lsl_outlet(info);
%Sounds to follow are show. This will be the baseline
% [time,type]=
StimuliIn(1000,outlet);
% timeStimuli=[timeStimuli time];
% typeStimuli=[typeStimuli type];
line1='Thank you for participating';
line2='\n This first part is';
line3='\n an explanation of the experiment';
Screen('TextSize', window, 80);
Screen('TextFont', window, 'Courier');
DrawFormattedText(window,[line1 line2 line3] , 'center', 'center', white);
Screen('Flip', window);
WaitSecs(5)
line1='Please read the explanation';
line2='\n hear the sounds';
line3='\n see the images';
line4='\n But do not tap during this part';
Screen('TextSize', window, 80);
Screen('TextFont', window, 'Courier');
DrawFormattedText(window,[line1 line2 line3 line4] , 'center', 'center', white);
Screen('Flip', window);
WaitSecs(5)

%The explanation/baseline starts
% outlet.push_sample(1000)
line1='During the experiment';
line2='\n you will tap everytime';
line3='\n you hear this sound';
Screen('TextSize', window, 80);
Screen('TextFont', window, 'Courier');
DrawFormattedText(window,[line1 line2 line3] , 'center', 'center', white);
Screen('Flip', window);
WaitSecs(5)
DrawFormattedText(window, '', 'center', 'center', white);
Screen('Flip', window);
freqSound=0.5;
time=3;
% [markers,timeMarkers]=
auditoryStimuli_v2(freqSound,nameSignal1,time,outlet);
% typeStimuli=[typeStimuli markers];
% timeStimuli=[timeStimuli timeMarkers];
WaitSecs(1)
% [markers,timeMarkers]=
auditoryStimuli_v2(freqSound,nameSignal1,time,outlet);
% typeStimuli=[typeStimuli markers];
% timeStimuli=[timeStimuli timeMarkers];

%We will show, for a few seconds, the presentation of the audio
line1='You will hear this sound';
line2='\n at different moments.';
DrawFormattedText(window,[line1 line2] , 'center', 'center', white);
Screen('Flip', window); 
WaitSecs(5)
line3='An example of this';
line4='\n is the following';
Screen('TextSize', window, 80);
Screen('TextFont', window, 'Courier');
DrawFormattedText(window,[line3 line4] , 'center', 'center', white);
Screen('Flip', window);
WaitSecs(5)
DrawFormattedText(window, '', 'center', 'center', white);
Screen('Flip', window);
freqSound=3.2;
time=5;
% [markers,timeMarkers]=
auditoryStimuli_v2(freqSound,nameSignal1,time,outlet);
% typeStimuli=[typeStimuli markers];
% timeStimuli=[timeStimuli timeMarkers];

line1='You can also hear it';
line2='\n Accompanied by other sounds,';
line3='\n Like in the following example'; 
Screen('TextSize', window, 80);
Screen('TextFont', window, 'Courier');
DrawFormattedText(window,[line1 line2 line3] , 'center', 'center', white);
Screen('Flip', window);
WaitSecs(5)
DrawFormattedText(window, '', 'center', 'center', white);
Screen('Flip', window);
freqSound=3.2;
time=5;
% tic
% [markers,timeMarkers]=
auditoryStimuliPoly_v3(freqSound,nameSignal2,nameSignal1,time,outlet);
% disp(['audi poly 3.2Hz 5s:',num2str(toc)])
% typeStimuli=[typeStimuli markers];
% timeStimuli=[timeStimuli timeMarkers];

% Here he show the images that they need to tap at
line1='You should also tap';
line2='\n everytime the next image appears,';
line3='\n During the experiment';
Screen('TextSize', window, 80);
Screen('TextFont', window, 'Courier');
DrawFormattedText(window,[line1 line2 line3] , 'center', 'center', white);
Screen('Flip', window);
WaitSecs(3)
Screen('DrawTexture', window, imageTexture, [], [], 0);
Screen('Flip', window);
WaitSecs(3)
line1='As with the sound, the image';
line2='\n will also appear at';
line3='\n different moments';
Screen('TextSize', window, 80);
Screen('TextFont', window, 'Courier');
DrawFormattedText(window,[line1 line2 line3] , 'center', 'center', white);
Screen('Flip', window);
WaitSecs(5)
line1='An example of this';
line2='\n is the following';
Screen('TextSize', window, 80);
Screen('TextFont', window, 'Courier');
DrawFormattedText(window,[line1 line2] , 'center', 'center', white);
Screen('Flip', window);
WaitSecs(2)
freqVisual=3.2;
time=5;
% [markers, timeMarkers]=
visualStimuli_Indiv_v3(freqVisual,time,window,ifi,screenYpixels,outlet,black);
% typeStimuli=[typeStimuli markers];
% timeStimuli=[timeStimuli timeMarkers];
WaitSecs(2)
line1='If you see two circles and a red dot,';
line2='\n focus your gaze on the red dot';
line3='\n but tap everytime ';
line4='\n the left circle appears';
Screen('TextSize', window, 80);
Screen('TextFont', window, 'Courier');
DrawFormattedText(window,[line1 line2 line3 line4] , 'center', 'center', white);
Screen('Flip', window);
WaitSecs(10)
Screen('DrawTexture', window, imageTexture, [], position1, 0);
Screen('DrawTexture', window, imageTexture, [], position2, 0);
Screen('DrawDots', window, [xCenter yCenter], 20, [1 0 0], [], 2);
Screen('Flip', window);
WaitSecs(5)
line1='An example of this';
line2='\n is the following';
Screen('TextSize', window, 80);
Screen('TextFont', window, 'Courier');
DrawFormattedText(window,[line1 line2] , 'center', 'center', white);
Screen('Flip', window);
WaitSecs(5)
freqVisual=3.2;
time=5;
visualStimuli_Poly_v3(freqVisual,time,window,ifi,xCenter,yCenter,screenXpixels, screenYpixels,outlet,black);
WaitSecs(2)
line1='The sounds/images will start';
line2='\n in 15 seconds.';
line3='\n Try to not move.';
line4='\n And prepare to start.';
Screen('TextSize', window, 80);
Screen('TextFont', window, 'Courier');
DrawFormattedText(window,[line1 line2 line3 line4] , 'center', 'center', white);
Screen('Flip', window);
WaitSecs(15)
% [time,type]=
StimuliIn(1001,outlet);
% timeStimuli=[timeStimuli time];
% typeStimuli=[typeStimuli type];
% Query the frame duration
ifi = Screen('GetFlipInterval', window);
DrawFormattedText(window, '', 'center', 'center', white);
Screen('Flip', window);
% outlet.push_sample(1000)
% [time,type]=
StimuliIn(1999,outlet);
% timeStimuli=[timeStimuli time];
% typeStimuli=[typeStimuli type];
for i=1:length(totalStimuli)
    
    
    selectedStimuli=stimuli(randperm(length(stimuli),1));
    stimuli(stimuli==selectedStimuli)=[];
%     disp(stimuli)
    
%     selectedStimuli
%     text='Tasks will start in 5 seconds';
%     Screen('DrawText', window, text, white, black);
    
    
    
%     [time,type]=
    StimuliIn(1018,outlet);
%     timeStimuli=[timeStimuli time];
%     typeStimuli=[typeStimuli type];
    switch selectedStimuli
        
%         pause(timeTask)
        %Auditory, 1Hz, isorhythmic
        
        case 1
            freqSound=1;
            time=16;

            for j=1:numTrials
%                 tic
%                 [markers,timeMarkers]=
                auditoryStimuli_v2(freqSound,nameSignal1,time,outlet);
%                 disp("Audi Iso 1: "+num2str(toc))
                pause(timeTrials)
%                 typeStimuli=[typeStimuli markers];
%                 timeStimuli=[timeStimuli timeMarkers];
            end
        DrawFormattedText(window, '', 'center', 'center', white);
        Screen('Flip', window);
        %Auditory, 3.2Hz, isorhythmic
%         [time,type]=
        StimuliIn(1999,outlet);
%         timeStimuli=[timeStimuli time];
%         typeStimuli=[typeStimuli type];
        pause(timeTask)
        case 2

            freqSound=3.2;
            time=16;
            for j=1:numTrials
%                 tic
%                 [markers,timeMarkers]=
                auditoryStimuli_v2(freqSound,nameSignal1,time,outlet);
%                 disp("Audi Iso 3.2: "+num2str(toc))
                
                pause(timeTrials)
%                 typeStimuli=[typeStimuli markers];
%                 timeStimuli=[timeStimuli timeMarkers];
            end
        DrawFormattedText(window, '', 'center', 'center', white);
        Screen('Flip', window);
        %Visual, 1Hz, isorhythmic
%         [time,type]=
        StimuliIn(1999,outlet);
%         timeStimuli=[timeStimuli time];
%         typeStimuli=[typeStimuli type];
        pause(timeTask)
        case 3
            freqVisual=1;
            time=16;
            for j=1:numTrials
%                 tic
%                 [markers,timeMarkers]=
                visualStimuli_Indiv_v3(freqVisual,time,window,ifi,screenYpixels,outlet,black);
%                 disp("visuial iso 1: "+num2str(toc))
                
                pause(timeTrials)
%                 typeStimuli=[typeStimuli markers];
%                 timeStimuli=[timeStimuli timeMarkers];
            end 
            %Visual, 3.2Hz, isorhythmic
        DrawFormattedText(window, '', 'center', 'center', white);
        Screen('Flip', window);    
%         [time,type]=
        StimuliIn(1999,outlet);
%         timeStimuli=[timeStimuli time];
%         typeStimuli=[typeStimuli type];
        pause(timeTask)
        case 4
            freqVisual=3.2;
            time=16;
            for j=1:numTrials
%                 tic
%                 [markers,timeMarkers]=
                visualStimuli_Indiv_v3(freqVisual,time,window,ifi,screenYpixels,outlet,black);
%                 disp("visual Iso 3.2: "+num2str(toc))
             
                pause(timeTrials)
%                 typeStimuli=[typeStimuli markers];
%                 timeStimuli=[timeStimuli timeMarkers];
            end 
        DrawFormattedText(window, '', 'center', 'center', white);
        Screen('Flip', window);    
        %Auditory, 1Hz, polyrhythmic  
%         [time,type]=
        StimuliIn(1999,outlet);
%         timeStimuli=[timeStimuli time];
%         typeStimuli=[typeStimuli type];
        pause(timeTask)
        case 5
           
            freqSound=1;
            time=16;
            for j=1:numTrials
%                 tic
%                 [markers,timeMarkers]=
                auditoryStimuliPoly_v3(freqSound,nameSignal2,nameSignal1,time,outlet);
%                 disp("Audi poly 1: "+num2str(toc))
                pause(timeTrials)
%                 typeStimuli=[typeStimuli markers];
%                 timeStimuli=[timeStimuli timeMarkers];
            end 
                
            
        %Auditory, 3.2Hz, polyrhythmic
        DrawFormattedText(window, '', 'center', 'center', white);
        Screen('Flip', window);
%         [time,type]=
        StimuliIn(1999,outlet);
%         timeStimuli=[timeStimuli time];
%         typeStimuli=[typeStimuli type];
        pause(timeTask)
        case 6
           freqSound=3.2;
            time=16;
            for j=1:numTrials
%                 tic
%                 [markers,timeMarkers]=
                auditoryStimuliPoly_v3(freqSound,nameSignal2,nameSignal1,time,outlet);
%                 disp("Audi poly 3.2: "+num2str(toc))
                pause(timeTrials)
%                 typeStimuli=[typeStimuli markers];
%                 timeStimuli=[timeStimuli timeMarkers];
            end  
           
        %Visual, 1Hz, polyrhythmic
        DrawFormattedText(window, '', 'center', 'center', white);
        Screen('Flip', window);
%         [time,type]=
        StimuliIn(1999,outlet);
%         timeStimuli=[timeStimuli time];
%         typeStimuli=[typeStimuli type];
        pause(timeTask)
        case 7
            freqVisual=1;
            time=16;
            for j=1:numTrials
                tic
%                 [markers,timeMarkers]=
                visualStimuli_Poly_v3(freqVisual,time,window,ifi,xCenter,yCenter,screenXpixels, screenYpixels,outlet,black);
%                 disp("visual poly 1: "+num2str(toc))
                pause(timeTrials)
%                 typeStimuli=[typeStimuli markers];
%                 timeStimuli=[timeStimuli timeMarkers];
            end 
        DrawFormattedText(window, '', 'center', 'center', white);
        Screen('Flip', window);
        %Visual, 3.2Hz, polyrhythmic
%         [time,type]=
        StimuliIn(1999,outlet);
%         timeStimuli=[timeStimuli time];
%         typeStimuli=[typeStimuli type];
        pause(timeTask)
        case 8
            freqVisual=3.2;
            time=16;
            for j=1:numTrials
%                 tic
%                 [markers,timeMarkers]=
                visualStimuli_Poly_v3(freqVisual,time,window,ifi,xCenter,yCenter,screenXpixels, screenYpixels,outlet,black);
%                 disp("visual poly 3.2: "+num2str(toc))
                pause(timeTrials)
%                 typeStimuli=[typeStimuli markers];
%                 timeStimuli=[timeStimuli timeMarkers];
            end 
        DrawFormattedText(window, '', 'center', 'center', white);
        Screen('Flip', window);
%         [time,type]=
        StimuliIn(1999,outlet);
%         timeStimuli=[timeStimuli time];
%         typeStimuli=[typeStimuli type];
        pause(timeTask)
        
        
    end 
    selectedStimuli=0;
end 
StimuliIn(1019,outlet);
% timeStimuli=[timeStimuli time];
% typeStimuli=[typeStimuli type]; 
line1='The experiment is finish';
line2= '\n Thank you for participating';
DrawFormattedText(window, [line1 line2], 'center', 'center', white);
Screen('Flip', window);  
WaitSecs(4)
save('varibles.mat')
outlet.delete();

sca;