%% EXP_AUDITORYSTEADYSTATEREPONSE

close all; clear all; clc;

%% SELECT SETTINGS

% indicate sessions
SEAT        = 'yes';        % seated session yes/no
WALK        = 'no';         % walking session yes/no

% indicate nr of trials/stimuli
nrRuns      = 3;            % nr of runs per session
nrStim      = 25;           % nr of stimuli per run (approximately)

% stimuli type
stimtype = 'ASSR_AM40Hz_500CF.wav'; 

% indicate path to save data
saveData = 'yes';               % save data yes/no
savePath = 'D:\dryEEG_data\ID11';  % subject data folder 

%% LSL OUTLET SENDING EVENTS

% load lsl library and make a new stream outlet
lib = lsl_loadlib();            % load labstreaming layer library
info = lsl_streaminfo(lib,'ASSR','Markers',1,0.0,'cf_int32','sdfwerr32432');
outlet = lsl_outlet(info);      % open an outlet

% create marker id's 
Marker_assr = 1266;             % ASSR stimuli
if strcmp(WALK,'no')==1 && strcmp(SEAT,'yes')==1
    Marker_startS = 1333;       % start signal SEATED
    Marker_stopS = 1300;        % stop signal SEATED
elseif strcmp(WALK,'yes')==1 && strcmp(SEAT,'no')==1
    Marker_startW = 1444;       % start signal WALKING
    Marker_stopW = 1400;        % stop signal WALKING
elseif strcmp(WALK,'yes')==1 && strcmp(SEAT,'yes')==1
    Marker_startS = 1333;       % start signal SEATED    
    Marker_stopS = 1300;        % stop signal SEATED
    Marker_startW = 1444;       % start signal WALKING
    Marker_stopW = 1400;        % stop signal WALKING
end

%% LOAD BEEP & CREATE STRUCT
% > generate SS-struct (sound stimuli)
%   > wavedata = sound data
%   > fs = sample frequency [Hz]
%   > nrChan = number of channels [#]

[SSassr]    = CreateWAVstruct (stimtype);
[SSstart]   = CreateWAVstruct ('WAVstart.wav');
[SSstop]    = CreateWAVstruct ('WAVstop.wav');

Lstim       = length(SSassr.wavedata)/SSassr.fs;    % stimuli length [s]

%% CREATE AND FILL AUDIO BUFFER

% load PsychPortAudio sound driver
% > high precision, low latency, Multichannel sound playback and recording
InitializePsychSound(1);
priority = 0;                       % 0 = better quality, increased latency; 1 = minimum latency
duration = 1;                       % number of repetitions of the wav-file
PsychPortAudio('Verbosity',1);      % verbosity = "wordiness" -> 1= print errors

% Get audio device
h_device = PsychPortAudio ('GetDevices');

% Open handle
h_assr   = PsychPortAudio('Open', [], [], priority, SSassr.fs, SSassr.nrChan);
h_start  = PsychPortAudio('Open', [], [], priority, SSstart.fs, SSstart.nrChan);
h_stop   = PsychPortAudio('Open', [], [], priority, SSstop.fs, SSstop.nrChan);

% Fill buffer
PsychPortAudio('FillBuffer', h_assr, SSassr.wavedata);
PsychPortAudio('FillBuffer', h_start, SSstart.wavedata);
PsychPortAudio('FillBuffer', h_stop, SSstop.wavedata);

%% SET NR OF STIMULI PER RUN
% > To keep attention to the task, subjects will be asked to count the
% > number of ASSR stimuli given during a run. Therefore the number of
% > stimuli should vary between runs.

nrStimSeat = []; nrStimWalk = [];
if strcmp(SEAT,'yes')==1
    for iRun = 1:nrRuns-1
        nrStimSeat(iRun,:) = randi([nrStim-5 nrStim+5]);
    end
    nrStimSeat(nrRuns,:) = (nrStim*nrRuns)-sum(nrStimSeat);
end

if strcmp(WALK,'yes')==1
    for iRun = 1:nrRuns-1
        nrStimWalk(iRun,:) = randi([nrStim-5 nrStim+5]);
    end
    nrStimWalk(nrRuns,:) = (nrStim*nrRuns)-sum(nrStimWalk);
end

%% GENERAL VARIABLES
% > 3 Runs seated (25x3=75 beeps) >> session = 0
% > 3 Runs gait (25x3=75 beeps) >> session = 1
% > randomized order

idxRunSeat = 1;                         % idx nr of runs in seat-session
idxRunWalk = 1;                         % idx nr of runs in walk-session
ses         = eye(2);                   % to randomize the order of sessions ...
idxWS       = ses(randi([1 2],1,1),:);  %   idx for order of sessions (flip[0,1])
clear ses

t_baseline  = 10;                       % time baseline [s]

% RunDets = details of the run (empty struct)
runDets = struct('nrID',[],'nrFile',[],'session',[],'nrStim',[],'timeStart',[],'timeStop',[],'timeBeeps',[],'waitbeeps',[]);

%% >> EXPERIMENT << %%
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% >>> WALKING & SEATED <<<   
if strcmp(WALK,'yes')==1 && strcmp(SEAT,'yes')==1
    
    nrID = num2str( input('Enter ID number: '), '%02d');                                % input ID number
    
    for iRun = 1:nrRuns*2
        
        if idxWS(1) == 1 % >>>>>>>>>> walking
            runDets.nrID = nrID;                                                        % fill details of run (Dets)
            runDets.nrFile = iRun;
            runDets.session = 'ASSR_walk';
            runDets.nrStim = nrStimWalk(idxRunWalk);
            fprintf('\nRun %g -   session: WALKING: ASSR paradigm \n\n', iRun);
            
            %%%>> START RUN <<%%%
            fprintf('press ENTER to start run\n'); KbStrokeWait;                        % press enter to start Run
            runDets.timeStart = GetSecs(PsychPortAudio('Start', h_start, 1, [], []));   % 'start'-sound 
            outlet.push_sample(Marker_startW);                                          % send LSL-marker (start)
            WaitSecs(t_baseline);                                                       % wait 10s for first stimuli
            
            %%%>> ASSR <<%%%
            for iStim = 1:nrStimWalk(idxRunWalk)
                runDets.timeBeeps(iStim,1) = GetSecs(PsychPortAudio('Start', h_assr, 1, [], []));   % play stimuli, GetSecs
                outlet.push_sample(Marker_assr);                                        % send LSL-marker (ASSR)
                runDets.waitBeeps(iStim,1) = WaitSecs(Lstim+2+rand(1));                 % interstim interval = 2.5 +/- 0.5
                fprintf('number %g out of %g stimuli\n', iStim, nrStimWalk(idxRunWalk));
            end
            
            WaitSecs(t_baseline/2);                                                     % wait 5s to end run
            runDets.timeStop = GetSecs(PsychPortAudio('Start', h_stop, 1));             % 'stop'-sound 
            outlet.push_sample(Marker_stopW);                                           % send LSL-marker (stop)
            
            %%%>> SAVE DATA <<%%%
            if strcmp(saveData,'yes')==1
                saveName = ['ID',runDets.nrID,'_',runDets.session,'_T',mat2str(idxRunWalk)];
                save([savePath, saveName], 'runDets');
            end
            idxRunWalk = idxRunWalk + 1;
            
            %%%>> END RUN <<%%%
            
        elseif idxWS(1) == 2 % >>>>>>>>>> seated
            runDets.nrID = nrID;                                                        % fill details of run (Dets)
            runDets.nrFile = iRun;
            runDets.session = 'ASSR_seat';
            runDets.nrStim = nrStimSeat(idxRunSeat);
            fprintf('\nRun %g -   session: SEATED: ASSR paradigm \n\n', iRun);
            
            %%%>> START RUN <<%%%
            fprintf('press ENTER to start run\n'); KbStrokeWait;                        % press enter to start Run
            runDets.timeStart = GetSecs(PsychPortAudio('Start', h_start, 1, [], []));   % 'start'-sound 
            outlet.push_sample(Marker_startS);                                          % send LSL-marker (start)
            WaitSecs(t_baseline);                                                       % wait 10s for first stimuli
                        
            %%%>> ASSR <<%%%
            for iStim = 1:nrStimSeat(idxRunSeat)
                runDets.timeBeeps(iStim,1) = GetSecs(PsychPortAudio('Start', h_assr, 1, [], [])); % play stimuli, GetSecs
                outlet.push_sample(Marker_assr);                                        % send LSL-marker (ASSR)
                runDets.waitBeeps(iStim,1) = WaitSecs(Lstim+2+rand(1));                 % interstim interval = 2.5 +/- 0.5
                fprintf('number %g out of %g stimuli\n', iStim, nrStimSeat(idxRunSeat));
            end
            
            WaitSecs(t_baseline/2);                                                     % wait 5s to end run
            runDets.timeStop = GetSecs(PsychPortAudio('Start', h_stop, 1));             % 'stop'-sound
            outlet.push_sample(Marker_stopS);                                           % send LSL marker (stop)
            
            %%%>> SAVE DATA <<%%%
            if strcmp(saveData,'yes')==1
                saveName = ['ID',runDets.nrID,'_',runDets.session,'_T',mat2str(idxRunSeat)];
                save([savePath, saveName], 'runDets');
            end
            idxRunSeat = idxRunSeat + 1;
            
            %%%>> END RUN <<%%%
            
        end
        idxWS = flip(idxWS);
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% >>> SEATED <<<       
elseif strcmp(WALK,'no')==1 && strcmp(SEAT,'yes')==1
    
    
    nrID = num2str( input('Enter ID number: '), '%02d');                                % input ID number
        
    for iRun = 1:nrRuns
        runDets.nrID = nrID;                                                            % fill details of run (Dets)
        runDets.nrFile = iRun;
        runDets.session = 'ASSR_seat';
        runDets.nrStim = nrStimSeat(idxRunSeat);
        fprintf('\nRun %g -   session: SEATED: ASSR paradigm \n\n', iRun);
        
        %%%>> START RUN <<%%%
        fprintf('press ENTER to start run\n'); KbStrokeWait;                            % press enter to start Run
        runDets.timeStart = GetSecs(PsychPortAudio('Start', h_start, 1, [], []));       % 'start'-sound
        outlet.push_sample(Marker_startS);                                              % send LSL-marker (start)
        WaitSecs(t_baseline);                                                           % wait 10s for first stimuli
                    
        %%%>> ASSR <<%%%
        for iStim = 1:nrStimSeat(idxRunSeat)
            runDets.timeBeeps(iStim,1) = GetSecs(PsychPortAudio('Start', h_assr, 1, [], [])); % play stimuli, GetSecs
            outlet.push_sample(Marker_assr);                                            % send LSL-trigger (ASSR)
            runDets.waitBeeps(iStim,1) = WaitSecs(Lstim+2+rand(1));                     % interstim interval = 2.5 +/- 0.5
            fprintf('number %g out of %g stimuli\n', iStim, nrStimSeat(idxRunSeat));
        end
        
        WaitSecs(t_baseline/2);                                                         % wait 5s to end run
        runDets.timeStop = GetSecs(PsychPortAudio('Start', h_stop, 1));                 % 'stop'-sound
        outlet.push_sample(Marker_stopS);                                               % send LSL-marker (stop)
        
        %%%>> SAVE DATA <<%%%
        if strcmp(saveData,'yes')==1
            saveName = ['ID',runDets.nrID,'_',runDets.session,'_T',mat2str(idxRunSeat)];
            save([savePath, saveName], 'runDets');
        end
        idxRunSeat = idxRunSeat + 1;
    end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% >>> WALKING <<<       

elseif strcmp(WALK,'yes')==1 && strcmp(SEAT,'no')==1
    
    % get subject ID number
    nrID = num2str( input('Enter ID number: '), '%02d');                                % input ID number

    for iRun = 1:nrRuns
        runDets.nrID = nrID;                                                            % fill details of run (Dets)
        runDets.nrFile = iRun;
        runDets.session = 'ASSR_walk';
        runDets.nrStim = nrStimWalk(idxRunWalk);
        fprintf('\nRun %g -   session: WALKING: ASSR paradigm \n\n', iRun);
        
        %%%>> START RUN <<%%%
        fprintf('press ENTER to start run\n'); KbStrokeWait;                            % press enter to start Run
        runDets.timeStart = GetSecs(PsychPortAudio('Start', h_start, 1, [], []));       % 'start'-sound
        outlet.push_sample(Marker_startW);                                              % send LSL-marker (start)
        WaitSecs(t_baseline);                                                           % wait 10s for first stimuli
        
        %%%>> ASSR <<%%%    
        for iStim = 1:nrStimWalk(idxRunWalk)
            runDets.timeBeeps(iStim,1) = GetSecs(PsychPortAudio('Start', h_assr, 1, [], [])); % play stimuli, GetSecs
            outlet.push_sample(Marker_assr);                                            % send LSL-marker (ASSR)
            runDets.waitBeeps(iStim,1) = WaitSecs(Lstim+2+rand(1));                     % interstim interval = 2.5 +/- 0.5 
            fprintf('number %g out of %g stimuli\n', iStim, nrStimWalk(idxRunWalk));
        end
        
        WaitSecs(t_baseline/2);                                                         % wait 5s to end run
        runDets.timeStop = GetSecs(PsychPortAudio('Start', h_stop, 1));                 % 'stop'-sound
        outlet.push_sample(Marker_stopW);                                               % send LSL-marker (stop)
        
        %%%>> SAVE DATA <<%%%
        if strcmp(saveData,'yes')==1
            saveName = ['ID',runDets.nrID,'_',runDets.session,'_T',mat2str(idxRunWalk)];
            save([savePath, saveName], 'runDets');
        end
        idxRunWalk = idxRunWalk + 1;
        
        %%%>> END RUN <<%%%
    end
end
PsychPortAudio('Close')

%% NOTES-SCRIPT

% info = lsl_streaminfo([lib handle],[name],[type],[channelcount],[fs],[channelformat],[sourceid])
%   > name = name of stream; describes device/product 
%   > type = content type of stream (EEG, Markers)
%   > channelcount = nr of channels per sample
%   > fs = samplking rate (Hz) as advertized by data source 
%   > channelformat = cf_float32, cf__double64, cf_string, cf_int32, cf_int16
%   > sourceid = unique identifier for source or device, if available


