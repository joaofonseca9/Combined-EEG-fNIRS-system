function [markers, timeMarkers]=auditoryStimuliPoly_v3(freqSound,nameSignal1,nameSignal2,time,outlet)
%%
% LOAD BEEP & CREATE STRUCT
% > generate SS-struct 
%   > wavedata = sound data 
%   > fs = sample frequency [Hz]
%   > nrChan = number of channels [#]
%CreateWAVstruct is an external function
%There is a confussion in the name of the signal. 146 should be nameSignal2
%Nevertheless, leave like this because everything was build with the error.
%So it's running correctly. Only the names are not correct
[Beep146]    = CreateWAVstruct (nameSignal1); %This really is the secondary sound
[Beep880]   = CreateWAVstruct (nameSignal2); %This is really the primary
% [BeepSum]=CreateWAVstruct (nameSignal3);
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
% fs=0;  %This value was selected as 4 times the frequency the highest beat.
% %In this case: secondary beat 3/2*main and if main is 3.2 this will be 4.8
% %This multiply by 4 is 19.2. So we use 20.
% id='sdfwerr32432';
% info = lsl_streaminfo(lib,'EOEC','Markers',1,fs,'cf_int32',id);
% outlet = lsl_outlet(info); % thing you push your data through
% info2 = lsl_streaminfo(lib,'TimeAudyPoly','Markers',1,fs,'cf_int32',id);
% outlet2 = lsl_outlet(info2); % thing you push your data through
%% Definition markers
A1PIM=1104; %This is the auditory isorhythmic 1 Hz initial marker
A1PFM=1112;%This is the auditory isorhythmic 1 Hz final marker
A3PIM=1106;%This is the auditory isorhythmic 3.2 Hz initial marker
A3PFM=1114;%This is the auditory isorhythmic 3.2 Hz final marker
A1PIS=1120; %This is the auditory isorhythmic 1 Hz initial marker
A1PFS=1128;%This is the auditory isorhythmic 1 Hz final marker
A3PIS=1122;%This is the auditory isorhythmic 3.2 Hz initial marker
A3PFS=1130;%This is the auditory isorhythmic 3.2 Hz final marker
%% CREATE AND FILL AUDIO BUFFER
%Two vector to save the markers and the time of the markers in Matlab
markers=[];
timeMarkers=[];
% load PsychPortAudio sound driver
% > high precision, low latency, Multichannel sound playback and recording
if freqSound==3.2
    markerstart=A3PIM;
    startTask=1008;
    endTask=1009;

else
    markerstart=A1PIM;
    startTask=1004;
    endTask=1005;

    
end 
InitializePsychSound(1);
repetitions=1; % number of repetitions of the wav-file
% numberBeeps=freqSound*time-1;
% freq2=3/2*freqSound;
latency = 1;                       % 0 = don't care about latency or timing;Level1=Try to get the lowest latency
%that is possible under the constraint of reliable playback, freedom of choice
%for all parameters and interoperability with other applications. Level2=Take full control over the audio device, even if this causes other sound
%applications to fail or shutdown 
freqMain=freqSound; %At the end the ratio of the rhythms will be [2:3] 
%This means that for a freqSound of 1Hz, we will have one beat for the
%main beat at two for the secondary in one second. For two seconds, we will
%have two main beats and three secondary beats. 
beepBetween= (1/freqMain)/2; % Length of the pause between beeps (main beat) 
timeVector1=0:beepBetween:(time); %Time vector associated with the main beat 
if freqSound==3.2
    timeVector1=0:beepBetween:(time); %Time vector associated with the main beat 
else
    timeVector1=0:beepBetween:(time); %Time vector associated with the main beat 
end

% beepPauseTime2=1/freq2;% Length of the pause between beeps (secondary beat) 
% timeVector2=0:beepPauseTime2:time; %Time vector associated with secondary beat
% startCue=0; %Start inmediatly when running the code
%This for is to generate the timeVector for the secondary beat. The time
%will be half of the time associated with the main frequency, after the
%time associated with main frequency and so on
waitForDeviceStart = 1;
PsychPortAudio('Verbosity',1);      % verbosity = "wordiness" -> 1= print errors

% Get audio device
h_device = PsychPortAudio ('GetDevices'); %Array of structures. Each structure associated with one PortAudioDevice

% Open handle
h_Beep146   = PsychPortAudio('Open', [], [], latency, Beep146.fs, Beep146.nrChan);
h_Beep880=PsychPortAudio('Open', [], [], latency, Beep880.fs, Beep880.nrChan);


TimefreqInit=1/freqSound; %This is the time in between the main beat 
%this variable should be the frequency at which we want to create the 

% Fill buffer
PsychPortAudio('FillBuffer', h_Beep146, Beep146.wavedata);
PsychPortAudio('FillBuffer', h_Beep880, Beep880.wavedata);


timeVectorTot=timeVector1;
numberBeeps=freqSound*time;
secondsynch=0;
i=1;
m=0;
% disp(['freq:',num2str(freqMain),' time:',num2str(time)])
% disp(['numberBeeps:',num2str(numberBeeps),' timeV:',num2str(length(timeVector1))])
while m<=numberBeeps
    
    
    if timeVectorTot(i)==0
        startCue=0;
        m=m+1;
        outlet.push_sample(startTask)
        outlet.push_sample(markerstart)
        PsychPortAudio('Start',h_Beep880,repetitions,startCue,waitForDeviceStart);
        PsychPortAudio('Start',h_Beep146,repetitions,startCue,waitForDeviceStart);
        PsychPortAudio('Stop', h_Beep880, 1, 1);
        [actualStartTime2, ~, ~, estStopTime2] = PsychPortAudio('Stop', h_Beep146, 1, 1);
        i=i+1;
        
        
    else
        if rem(timeVectorTot(i),TimefreqInit)==0 && secondsynch==1
            startCue=actualStartTime2+timeVectorTot(i);
            
            outlet.push_sample(markerstart)


            PsychPortAudio('Start', h_Beep880, repetitions, startCue, waitForDeviceStart);
            PsychPortAudio('Start',h_Beep146,repetitions,startCue,waitForDeviceStart);

            
            PsychPortAudio('Stop', h_Beep880, 1, 1);
            PsychPortAudio('Stop', h_Beep146, 1, 1);
            secondsynch=0;
            m=m+1;
            
            i=i+1;
        elseif rem(timeVectorTot(i),TimefreqInit)==0 && secondsynch==0
            startCue=actualStartTime2+timeVectorTot(i);
            
            outlet.push_sample(markerstart)

            PsychPortAudio('Start', h_Beep880, repetitions, startCue, waitForDeviceStart);
            PsychPortAudio('Stop', h_Beep880, 1, 1);
            m=m+1;
            secondsynch=1;
            i=i+1;

        else
            startCue=actualStartTime2+timeVectorTot(i);

            PsychPortAudio('Start', h_Beep146, repetitions, startCue, waitForDeviceStart);
            PsychPortAudio('Stop', h_Beep146, 1, 1);
            i=i+1;

        end
        
    end
    if KbCheck
        sca;
        return
        
    end
    % disp(n) 
end
% outlet.delete()
% outlet2.delete()
outlet.push_sample(endTask)
end 
