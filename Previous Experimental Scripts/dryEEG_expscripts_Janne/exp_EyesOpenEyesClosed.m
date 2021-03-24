%% EXP_EyesOpenEyesClosed

clear all; close all; clc;

%% SELECT SETTINGS

% indicate sessions
SEAT        = 'yes';            % seated session yes/no
WALK        = 'no';             % walking session yes/no

% indicate path to save data
saveData = 'yes';               % save data yes/no
savePath = 'D:\internship I\Measurements\Subject_1';  % subject data folder 

% indicate nr of trials/stimuli
nrRuns      = 2;                % nr of runs per session (2 in exp)
nrStim      = 30;               % nr of stimuli per run (30 in exp)

%% LSL OUTLET SENDING EVENTS

% load lsl library and make a new stream outlet
lib = lsl_loadlib();            % load labstreaming layer library
info = lsl_streaminfo(lib,'EOEC','Markers',1,0.0,'cf_int32','sdfwerr32432');
outlet = lsl_outlet(info);      % open an outlet

% create marker id's 
Marker_EO = 1233;               % eyes open
Marker_EC = 1244;               % eyes closed
if strcmp(WALK,'no')==1 && strcmp(SEAT,'yes')==1
    Marker_startS = 1111;       % start signal SEATED
    Marker_stopS = 1100;        % stop signal SEATED
elseif strcmp(WALK,'yes')==1 && strcmp(SEAT,'no')==1
    Marker_startW = 1222;       % start signal WALKING
    Marker_stopW = 1200;        % stop signal WALKING
elseif strcmp(WALK,'yes')==1 && strcmp(SEAT,'yes')==1
    Marker_startS = 1111;       % start signal SEATED
    Marker_stopS = 1100;        % stop signal SEATED
    Marker_startW = 1222;       % start signal WALKING
    Marker_stopW = 1200;        % stop signal WALKING
end

%% LOAD VOICE COMMANDS

[WAVstart.wave,WAVstart.fs]     = audioread('WAVstart.wav');
[WAVstop.wave,WAVstop.fs]       = audioread('WAVstop.wav');
[WAVopen.wave,WAVopen.fs]       = audioread('WAVopen.wav');
[WAVclose.wave,WAVclose.fs]     = audioread('WAVdicht.wav');

% change rows<>columns & double channels
WAVstart.wave = WAVstart.wave';     WAVstart.wave(2,:) = WAVstart.wave(1,:);     WAVstart.nrChan = 2;
WAVstop.wave = WAVstop.wave';       WAVstop.wave(2,:) = WAVstop.wave(1,:);       WAVstop.nrChan = 2;
WAVopen.wave = WAVopen.wave';       WAVopen.wave(2,:) = WAVopen.wave(1,:);       WAVopen.nrChan = 2;
WAVclose.wave = WAVclose.wave';     WAVclose.wave(2,:) = WAVclose.wave(1,:);     WAVclose.nrChan = 2;

%% CREATE AND FILL AUDIO BUFFER

% Initialize Sounddriver
% This routine loads the PsychPortAudio sound driver for high precision, low latency, 
% multichannel sound playback and recording
% Call it at the beginning of your experiment script, optionally providing the 
% 'reallyneedlowlatency'-flag set to 1 to push really hard for low latency
InitializePsychSound(1);

priority = 0;                       % 0 = better quality, increased latency; 1 = minimum latency
duration = 1;                       % number of repetitions of the wav-file
PsychPortAudio('Verbosity',1);      % verbosity = "wordiness" -> 1= print errors

% Get audio device
h_device = PsychPortAudio ('GetDevices');

% Open handle
h_start  = PsychPortAudio('Open', [], [], priority, WAVstart.fs, WAVstart.nrChan);
h_stop   = PsychPortAudio('Open', [], [], priority, WAVstop.fs, WAVstop.nrChan);
h_open   = PsychPortAudio('Open', [], [], priority, WAVopen.fs, WAVopen.nrChan);
h_close  = PsychPortAudio('Open', [], [], priority, WAVclose.fs, WAVclose.nrChan);

% Fill buffer
PsychPortAudio('FillBuffer', h_start, WAVstart.wave);
PsychPortAudio('FillBuffer', h_stop, WAVstop.wave);
PsychPortAudio('FillBuffer', h_open, WAVopen.wave);
PsychPortAudio('FillBuffer', h_close, WAVclose.wave);

%% GENERAL VARIABLES
% > Eyes open (EO) / Eyes closed (EC)

idxRunSeat  = 1;                         % idx nr of runs in seat-session
idxRunWalk  = 1;                         % idx nr of runs in walk-session
idxEOEC     = [1 2];                     % idx (1) eyes open; (2) eyes closed
ses         = eye(2);                    % to randomize the order of sessions ...
idxWS       = ses(randi([1 2],1,1),:);   % idx for order of sessions (flip[0,1])
clear ses

t_baseline  = 10;                        % time baseline [s]
t_EOEC      = [5 3];                     % time eyes open / eyes closed [s]

% Dets - details of the run (empty struct)
runDets     = struct('nrID',[],'nrFile',[],'session',[],'timeStart',[],'timeStop',[],'timeEO',[],'timeEC',[]);

%% >> EXPERIMENT << %%
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% >>> SEATED <<<
if strcmp(WALK,'no')==1 && strcmp(SEAT,'yes')==1

    nrID = num2str( input('Enter ID number: '), '%02d');                                    % input ID number

    for iRun = 1:nrRuns
        runDets.nrID = nrID;                                                                % fill details of run (Dets)
        runDets.nrFile = iRun;
        runDets.session = 'EOEC_seat';
        fprintf('\nRun %g -   session: SEATED: EYES OPEN / EYES CLOSED\n\n', iRun);

        %%%>> START RUN <<%%%
        fprintf('press ENTER to start run\n'); KbStrokeWait;                                % press enter to start Run
        runDets.timeStart = GetSecs(PsychPortAudio('Start', h_start, 1, [], []));           % 'start'-sound
        outlet.push_sample(Marker_startS);                                                  % send LSL-marker (start)
        WaitSecs(t_baseline);                                                               % wait 10s for first beep

        for iStim = 1:nrStim
            % EYES CLOSED 
            runDets.timeEC(iStim,1) = GetSecs(PsychPortAudio('Start', h_close, 1));         % 'dicht'-sound
            outlet.push_sample(Marker_EC);                                                  % send LSL-marker (EC)
            WaitSecs(t_EOEC(1,2));                                                          % interstimuli interval 
            % EYES OPEN 
            runDets.timeEO(iStim,1) = GetSecs(PsychPortAudio('Start', h_open, 1));          % 'open'-sound
            outlet.push_sample(Marker_EO);                                                  % send LSL-marker (EO)
            WaitSecs(t_EOEC(1,1));                                                          % interstimuli interval

            fprintf('number %g out of %g stimuli\n', iStim, nrStim);
        end 

        WaitSecs(t_baseline/2);                                                             % wait 5s to end run
        runDets.timeStop = GetSecs(PsychPortAudio('Start', h_stop, 1));                     % 'stop'-sound
        outlet.push_sample(Marker_stopS);                                                   % send LSL-marker (stop)
        fprintf('\nEND TRIAL\n\n');
        
        % SAVE DATA
        if strcmp(saveData,'yes')==1
              saveName = ['ID',runDets.nrID,'_',runDets.session,'_T',mat2str(idxRunSeat)];
              save([savePath, saveName], 'runDets');
        end
        
        idxRunSeat = idxRunSeat + 1;
        
        %%%>> END RUN <<%%%

    end %iRun
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% >>> WALKING <<<   
elseif strcmp(WALK,'yes')==1 && strcmp(SEAT,'no')==1

    nrID = num2str( input('Enter ID number: '), '%02d');                                    % input ID number

    for iRun = 1:nrRuns
        runDets.nrID = nrID;                                                                % fill details of run (Dets)
        runDets.nrFile = iRun;
        runDets.session = 'EOEC_walk';
        fprintf('\nRun %g -   session: WALKING: EYES OPEN / EYES CLOSED\n\n', iRun);

        %%%>> START RUN <<%%%
        fprintf('press ENTER to start run\n'); KbStrokeWait;                                % press enter to start Run
        runDets.timeStart = GetSecs(PsychPortAudio('Start', h_start, 1, [], []));           % 'start'-sound
        outlet.push_sample(Marker_startW);                                                  % send LSL-marker (start)
        WaitSecs(t_baseline);                                                               % wait 10s for first beep

        for iStim = 1:nrStim
            % EYES CLOSED
            runDets.timeEC(nrStim,1) = GetSecs(PsychPortAudio('Start', h_close, 1));        % 'dicht'-sound
            outlet.push_sample(Marker_EC);                                                  % send LSL-marker (EC)
            WaitSecs(t_EOEC(1,2));                                                          % interstimuli interval
            % EYES OPEN
            runDets.timeEC(nrStim,1) = GetSecs(PsychPortAudio('Start', h_open, 1));         % 'open'-sound
            outlet.push_sample(Marker_EO);                                                  % send LSL-marker (EO)
            WaitSecs(t_EOEC(1,1));                                                          % interstimuli interval

            fprintf('number %g out of %g stimuli\n', iStim, nrStim);
        end 

        WaitSecs(t_baseline/2);                                                             % wait 5s to end run
        runDets.timeStop = GetSecs(PsychPortAudio('Start', h_stop, 1));                     % 'stop'-sound
        outlet.push_sample(Marker_stopW);                                                   % send LSL-marker (stop)
                
        % SAVE DATA
        if strcmp(saveData,'yes')==1
            saveName = ['ID',runDets.nrID,'_',runDets.session,'_T',mat2str(idxRunWalk)];
            save([savePath, saveName], 'runDets');
        end
        idxRunWalk = idxRunWalk + 1;
        
        %%%>> END RUN <<%%%

    end %iRun
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% >>> WALKING & SEATED <<<   
elseif strcmp(WALK,'yes')==1 && strcmp(SEAT,'yes')==1
   
    nrID = num2str( input('Enter ID number: '), '%02d');                                    % input ID number

    for iRun = 1:nrRuns*2
        
        if idxWS(1) == 1 % walking
            runDets.nrID = nrID;                                                            % fill details of run (Dets)
            runDets.nrFile = iRun;
            runDets.session = 'EOEC_walk';
            fprintf('\nRun %g -   session: WALKING: EYES OPEN / EYES CLOSED\n\n', iRun);

            %%%>> START RUN <<%%%
            fprintf('press ENTER to start run\n'); KbStrokeWait;                            % press enter to start Run
            runDets.timeStart = GetSecs(PsychPortAudio('Start', h_start, 1, [], []));       % 'start'-sound 
            outlet.push_sample(Marker_startW);                                              % send LSL-marker (start)
            WaitSecs(t_baseline);                                                           % wait 10s for first stimuli

            for iStim = 1:nrStim
                % EYES CLOSED
                runDets.timeEC(nrStim,1) = GetSecs(PsychPortAudio('Start', h_close, 1));    % 'dicht'-sound
                outlet.push_sample(Marker_EC);                                              % send LSL-marker (EC)
                WaitSecs(t_EOEC(1,2));                                                      % interstimuli interval
                % EYES OPEN
                runDets.timeEC(nrStim,1) = GetSecs(PsychPortAudio('Start', h_open, 1));     % 'open'-sound
                outlet.push_sample(Marker_EO);                                              % send LSL-marker (EO)
                WaitSecs(t_EOEC(1,1));                                                      % interstimuli interval

                fprintf('number %g out of %g stimuli\n', iStim, nrStim);
            end 

            WaitSecs(t_baseline/2);                                                         % wait 5s to end run
            runDets.timeStop = GetSecs(PsychPortAudio('Start', h_stop, 1));                 % 'stop'-sound
            outlet.push_sample(Marker_stopW);                                               % send LSL-marker (stop)

            % SAVE DATA
            if strcmp(saveData,'yes')==1
                saveName = ['ID',runDets.nrID,'_',runDets.session,'_T',mat2str(idxRunWalk)];
                save([savePath, saveName], 'runDets');
            end
            idxRunWalk = idxRunWalk + 1;
            
            %%%>> END RUN <<%%%

        elseif idxWS(1) == 2
            runDets.nrID = nrID;                                                            % fill details of run (Dets)
            runDets.nrFile = iRun;
            runDets.session = 'EOEC_seat';
            fprintf('\nRun %g -   session: SEATED: EYES OPEN / EYES CLOSED\n\n', iRun);
            
            %%%>> START RUN <<%%%
            fprintf('press ENTER to start run\n'); KbStrokeWait;                            % press enter to start Run
            runDets.timeStart = GetSecs(PsychPortAudio('Start', h_start, 1, [], []));       % 'start'-sound 
            outlet.push_sample(Marker_startS);                                              % send LSL-marker (start)
            WaitSecs(t_baseline);                                                           % wait 10s for first beep
            
            for iStim = 1:nrStim
                % EYES CLOSED 
                runDets.timeEC(nrStim,1) = GetSecs(PsychPortAudio('Start', h_close, 1));    % 'dicht'-sound
                outlet.push_sample(Marker_EC);                                              % send LSL-marker (EC)
                WaitSecs(t_EOEC(1,2));                                                      % interstimuli interval
                % EYES OPEN
                runDets.timeEC(nrStim,1) = GetSecs(PsychPortAudio('Start', h_open, 1));     % 'open'-sound
                outlet.push_sample(Marker_EO);                                              % send LSL-marker (EO)
                WaitSecs(t_EOEC(1,1));                                                      % interstimuli interval
                
                fprintf('number %g out of %g stimuli\n', iStim, nrStim);
            end 
            
            WaitSecs(t_baseline/2);                                                         % wait 5s to end run
            runDets.timeStop = GetSecs(PsychPortAudio('Start', h_stop, 1));                 % 'stop'-sound
            outlet.push_sample(Marker_stopS);                                               % send LSL-marker (stop)
            
            % SAVE DATA
            if strcmp(saveData,'yes')==1
                saveName = ['ID',runDets.nrID,'_',runDets.session,'_T',mat2str(idxRunSeat)];
                save([savePath, saveName], 'runDets');
            end
            idxRunSeat = idxRunSeat + 1;
            
            %%%>> END RUN <<%%%
            
        end %idxWS
        
        idxWS = flip(idxWS);
        
    end %iRun
    
end

PsychPortAudio('Close');

%% NOTES SCRIPT

% info = lsl_streaminfo([lib handle],[name],[type],[channelcount],[fs],[channelformat],[sourceid])
%   > name = name of stream; describes device/product 
%   > type = content type of stream (EEG, Markers)
%   > channelcount = nr of channels per sample
%   > fs = samplking rate (Hz) as advertized by data source 
%   > channelformat = cf_float32, cf__double64, cf_string, cf_int32, cf_int16
%   > sourceid = unique identifier for source or device, if available


