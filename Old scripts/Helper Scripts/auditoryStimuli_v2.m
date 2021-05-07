% Creation auditory stimulation
% This code is based on provided code by Janne Heijs.
function [markers,timeMarkers]=auditoryStimuli_v2(freqSound,nameSignal,time,outlet)
%%
%This section is a modification of the code developed by Janne Heijs
% load lab streaming layer library
lib = lsl_loadlib();
% make a new stream outlet

% info = lsl_streaminfo([lib handle],[name],[type],[channelcount],[fs],[channelformat],[sourceid])

%   > name = name of stream; describes device/product

%   > type = content type of stream (EEG, Markers)

%   > channelcount = nr of channels per sample

%   > fs = samplking rate (Hz) as advertized by data source

%   > channelformat = cf_float32, cf__double64, cf_string, cf_int32, cf_int16

%   > sourceid = unique identifier for source or device, if available
fs=0;
id='sdfwerr32432';
info = lsl_streaminfo(lib,'EOEC','Markers',1,fs,'cf_int32',id);
info2=lsl_streaminfo(lib,'TimeAudiIso','Markers',1,fs,'cf_int32',id);
outlet = lsl_outlet(info); % thing you push your data through/ This is the
%marker associated with the task
% outlet2 = lsl_outlet(info2); % thing you push your data through. / This is 
%the marker associated with the time.
%% Definition markers
A1IIM=1100; %This is the auditory isorhythmic 1 Hz initial marker
A1IFM=1108;%This is the auditory isorhythmic 1 Hz final marker
A3IIM=1102;%This is the auditory isorhythmic 3.2 Hz initial marker
A3IFM=1110;%This is the auditory isorhythmic 3.2 Hz final marker

%%
% LOAD BEEP & CREATE STRUCT
% > generate SS-struct 
%   > wavedata = sound data 
%   > fs = sample frequency [Hz]
%   > nrChan = number of channels [#]
%CreateWAVstruct is an external function
[Beep]    = CreateWAVstruct (nameSignal);
% [Beep800]   = CreateWAVstruct ('Beep800Hz.wav');

%% CREATE AND FILL AUDIO BUFFER

% load PsychPortAudio sound driver
% > high precision, low latency, Multichannel sound playback and recording
%We create two vector to save the markers and the times at which they
%appeared.
markers=[];
timeMarkers=[];
if freqSound==3.2
    markerstart=A3IIM;
    startTask=1006;
    endTask=1007;
%     markerend=A3IFM;
else
    markerstart=A1IIM;
    startTask=1002;
    endTask=1003;
%     markerend=A1IFM;
    
end 

InitializePsychSound(1);
repetitions=1; % number of repetitions of the wav-file
numberBeeps=freqSound*time;
latency = 1;                       % 0 = don't care about latency or timing;Level1=Try to get the lowest latency
%that is possible under the constraint of reliable playback, freedom of choice
%for all parameters and interoperability with other applications. Level2=Take full control over the audio device, even if this causes other sound
%applications to fail or shutdown 
beepPauseTime = 1/freqSound; % Length of the pause between beeps
if freqSound==3.2
    timeVector1=0:beepPauseTime:(time+1); %Time vector associated with the main beat 
else
    timeVector1=0:beepPauseTime:(time); %Time vector associated with the main beat 
end

startCue=0; %Start inmediatly when running the code

waitForDeviceStart = 1;
PsychPortAudio('Verbosity',1);      % verbosity = "wordiness" -> 1= print errors

% Get audio device
h_device = PsychPortAudio ('GetDevices'); %Array of structures. Each structure associated with one PortAudioDevice

% Open handle
h_Beep   = PsychPortAudio('Open', [], [], latency, Beep.fs, Beep.nrChan);
% h_SSstart  = PsychPortAudio('Open', [], [], latency, SSstart.fs, SSstart.nrChan);
% h_SSstop   = PsychPortAudio('Open', [], [], latency, SSstop.fs, SSstop.nrChan);

% Fill buffer
PsychPortAudio('FillBuffer', h_Beep, Beep.wavedata);

% disp(['Frequency:',num2str(freqSound),' time:',num2str(time),' beeps:',num2str(numberBeeps)])
% disp(['startCue: ',num2str(startCue),' I:out'])
% disp(['numberBeeps:',num2str(numberBeeps),' timeVector:',num2str(length(timeVector1))])
for i=1:(numberBeeps+1)
    
    if i == 1
        
        startCue = GetSecs;
        outlet.push_sample(startTask)
        outlet.push_sample(markerstart)
        PsychPortAudio('Start',h_Beep,repetitions,startCue,waitForDeviceStart)
        [actualStartTime, ~, ~, estStopTime] = PsychPortAudio('Stop', h_Beep, 1, 1);

    else

        % Compute new start time for follow-up beep, beepPauseTime after end of
        % previous one
        startCue = actualStartTime + timeVector1(i);

        outlet.push_sample(markerstart)
        PsychPortAudio('Start', h_Beep, repetitions, startCue, waitForDeviceStart);
        % Wait for stop of playback
        PsychPortAudio('Stop', h_Beep, 1, 1);
    end
%     toc
    if KbCheck
        sca;
        return
        
    end


end 
%Send the end task marker
outlet.push_sample(endTask)
% Close the audio device
PsychPortAudio('Close', h_Beep);
% outlet.delete()
end



