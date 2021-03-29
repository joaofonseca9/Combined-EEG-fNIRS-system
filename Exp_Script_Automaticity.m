%% EXPERIMENTAL SCRIPT FOR COMBINED EEG/FNIRS ANALYSYS %%

%TASKS IMPLEMENTED IN THIS CODE:

%Automaticity Task / Double Task paradigm:

%Checkerboard Task:

%Automatic & Non-automatic task / Single Task paradigm:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Settings:

% Clear the workspace and the screen
%sca; close all; clear all; clc

% MATLAB SCRIPT FOR THE AUTOMATICITY TEST
% A dual-task paradigm in which we test whether a 12 digit prelearned
% sequence (sequenceauto) has become an automatic movement for a finger
% tapping and a foot stomping task.

%Before starting the automaticity test, clear the workspace.
clear all

%Synch test skip => comment when actually testing patient
Screen('Preference', 'SkipSyncTests', 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SET UP PARAMETERS

%The sequences used for this study (automatic and non-automatic sequences
%randomized between participants)

%Sequences used in order to be able to print in the command window if
%to generate a new sequence use randi([1 4], 1, 12)
sequencesprint = {('4 3 4 1 4 1 2 4 3 2 1 2'),('2 1 2 3 2 1 3 2 4 2 4 1')};

sequences = {split(sequenceprint(1))',split(sequenceprint(2))'} ;

% % Set the right sequence that was studied at home (= automatic)
% sequenceauto = sequenceA;
% sequenceautoprint = sequenceprintA;

%Randomize the first sequence to be tested
%1=auto sequence, 2 non-auto sequence
sequence_idx=randi([1,2]);

%Parameters for the resting period in between the trials
t1 = 20; %Resting period in seconds
t2 = 5;  %Random interval around the resting period time

%Amount of letters presented during test for automaticity for one trial.
%Should be adjusted when letter presenting speed is changed!
N_letters=8; % 8 letters presented during a trial
N_trials=4; % number of trials 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LSL SETUP
% LSL outlet sending events

% Create and load the lab streaming layer library

addpath(genpath('C:\Users\joaop\Downloads\liblsl-Matlab'));
lib = lsl_loadlib(); version = lsl_library_version(lib);
lib = lsl_loadlib();

% Make a new stream outlet.
% info = lsl_streaminfo([lib handle],[name],[type],[channelcount],[fs],[channelformat],[sourceid])
% > name = name of stream; describes device/product
% > type = content type of stream (EEG, Markers)
% > channelcount = nr of channels per sample
% > fs = samplking rate (Hz) as advertized by data source
% > channelformat = cf_float32, cf__double64, cf_string, cf_int32, cf_int16
% > sourceid = unique identifier for source or device, if available
info    = lsl_streaminfo(lib, 'Dual Task', 'Markers', 1, 0.0, 'cf_int32', 'ReactionTime'); 

% Open an outlet for the data to run through.
outlet = lsl_outlet(info);

% MARKER SETUP
% Block related
instructions = 'instructions'; %NEVER USED (?)
finger_test='finger_test';
Marker_StartBlockCue1_5HzAddition       = 7000;         
Marker_EndBlockCue1_5HzAddition         = 7001;

Marker_StartBlock_AutomaticSequence     = 7002;
Marker_StartBlock_NonAutomaticSequence  = 7003;

Marker_EndBlock_AutomaticSequence       = 7012;
Marker_EndBlock_NonAutomaticSequence       = 7013;

%Open Pshychtoolbox.
PsychDefaultSetup(2);
KbName('UnifyKeyNames'); %Links the key presses to the key board names
KbQueueCreate;
KbQueueStart; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUDIO PREPARATION

%LOAD METRONOME SOUNDS (PsychToolbox)
audio_dir='.\Previous Experimental Scripts\Experiment_ME\metronomesounds';
cd(audio_dir)
[WAVMetronome8.wave,WAVMetronome8.fs]       = audioread('Metronome8.wav');

% change rows<>columns
WAVMetronome8.wave = WAVMetronome8.wave';         WAVMetronome8.nrChan=2;

%LOAD CUES

% Get Cueing Files
Cue1_5Hz       = 'Metronome120.wav';
[Cue1_5Hz]     = CreateWAVstruct(Cue1_5Hz);
Cue1_5HzLength = length(Cue1_5Hz.wavedata)/Cue1_5Hz.fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
PPA_device = PsychPortAudio ('GetDevices');

% Open handle
PPA_cue1Hz   = PsychPortAudio('Open', [], [], priority, WAVMetronome8.fs, WAVMetronome8.nrChan);
PPA_cue1_5Hz = PsychPortAudio('Open', [], [], priority, Cue1_5Hz.fs, Cue1_5Hz.nrChan);

% Fill buffer
PsychPortAudio('FillBuffer', PPA_cue1Hz, WAVMetronome8.wave);
PsychPortAudio('FillBuffer', PPA_cue1_5Hz, Cue1_5Hz.wavedata);

%AudioFile
file = [PPA_cue1Hz; PPA_cue1_5Hz];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE FILES IN FOLDER

fprintf('Select the project directory \n')
root_dir=uigetdir('C:\Users\mtabo\Documents\TryOutScript\', 'Select the project directory');

complete=0;
while complete==0
    sub_ID=input('What is the subject ID (2 digit number) \n', 's');
    sub=sprintf('sub-%s', sub_ID);
        rec_n=input('What is the number of the recording? \n');
        rec=sprintf('rec-%.2d', rec_n);

    inf=fprintf('\n root_dir = %s \n sub = %s \n rec = %s \n', root_dir, sub, rec);
    correct=input('Is the above information correct? (y/n) \n', 's');
    if strcmp(correct, 'y')
        complete=1;
    else
        continue
    end
end

% go to subject folder
sub_dir=fullfile(root_dir, sub);
if ~exist(sub_dir)
    mkdir(sub_dir)
end
cd(sub_dir)

logname=sprintf('%s_%s_triggers.log', sub, rec); diary(logname);
% save current script in subject directory
script=mfilename('fullpath');
script_name=mfilename;
copyfile(sprintf('%s.m', script), fullfile(sub_dir, sprintf('%s_%s.m', sub, script_name)))

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

% Preparations for the fixation cross so we only need to do this once
fixCrossDimPix = 40; % Here we set the size of the arms of our fixation cross
xCoords = [-fixCrossDimPix fixCrossDimPix 0 0]; % Set the coordinates (these are all relative to zero we will let the drawing routine center the cross in the center of our monitor for us)
yCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
allCoords = [xCoords; yCoords];
lineWidthPix = 4;% Set the line width for the fixation cross

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% START TEST FOR AUTOMATICITY

%Empty structure for key presses, -> use later again so it saves the key
%presses within this structure -> save at the end
events_handautodual=struct([]); % same for the presented letters of the hand + the answer
events_handautodual(1).sequence_label='Automatic';
events_handautodual(2).sequence_label='Non-automatic';
%Instruction automaticity test
Screen('TextSize',window,25);
DrawFormattedText(window,sprintf('You will now start with the automaticity test in a dual task situation. \n \n You will be tested on the sequence you learned at home, as well as the sequence you learned today. \n For each, you will be tested %d times. \n\n Detailed instructions will be given at the start of each task. \n Press any key to continue.',ceil(N_trials/2)),'center', 'center', white);
vbl = Screen('Flip', window);
KbStrokeWait; %wait for response to terminate instructions

Screen('TextSize',window,25);
DrawFormattedText(window, sprintf('While you perform the task, letters will be shown on the screen (A,G,O,L). \n The goal is to perform the sequence tapping while counting how many times G is presented. \n After each time you tapped the full sequence, you should tell us how many times G was presented. \n For answering this question, \n keep in mind that when the answer is 4 you press 4 and not Return (Enter) on the keyboard. \n\n Note that during the tapping task you cannot talk. \n Try to keep your body movements as still as possible exept for the right hand. \n Keep your eyes open (also during the rest periods). \n\n In between the trials you will see a fixation cross for 20 seconds. \n During the first few seconds you will hear a metronome sound. \n Tap the sequence on this rhythm, which is the same as you studied at home. \n\n After the fixation cross, the first trial will start automatically. \n So start tapping the sequence as soon as a letter on the screen appears. \n When ready: press any key.'),'center','center', white);
vbl = Screen('Flip', window);
KbStrokeWait; %wait for response to terminate instructions

%Demonstrate how cueing works
Screen('TextSize',window,25);
DrawFormattedText(window,'While you perform the finger tapping sequence and count the letters, \n you may hear a rhytmic sound like this one. \n Or you may just hear nothing.','center', 'center', white);
vbl = Screen('Flip', window);
PsychPortAudio('Start', file(2), 1, [], []);
WaitSecs(5);
PsychPortAudio('Stop', file(2));
% Finger tapping test -> 20 trials, dual task, in which a participant
% taps a prelearned sequence, while also letters are presented on the
% screen in a randomized speed

for ii=1:2
    %Instruction automaticity task finger tapping
    % trig.beep(440, 0.2, 'instructions');
    Screen('TextSize',window,25);
    
    if sequence_idx==1
        DrawFormattedText(window, sprintf('You will now perform the pre-learned sequence for the FINGER tapping task: \n %s \n\n  When ready to start: press any key.', char(sequencesprint(sequence_idx))),'center','center', white);
        outlet.push_sample(Marker_StartBlock_AutomaticSequence);
    else
        DrawFormattedText(window, sprintf('You will now perform the sequence you learned today for the FINGER tapping task: \n %s \n\n  When ready: press any key.', char(sequencesprint(sequence_idx))),'center','center', white);
        outlet.push_sample(Marker_StartBlock_NonAutomaticSequence);
    end
    vbl = Screen('Flip', window);
    KbStrokeWait; %wait for response to terminate instructions

    %Start loop for the trials
    for j=1:N_trials/2
        %Presentation of the letters on the screen (dual task). -> is random.
        %Participant has to count the times that G was presented.
        Letterlist='AGOL';
        letter_order=randi(length(Letterlist), 1, N_letters);
        value={Letterlist(letter_order)};

        endOfTrial=0; %helper variable to terminate cueing
        %Always start with a 20-25 seconds fixation cross with 8 seconds of metronome
        %sound
    %     trig.beep(440, 0.2, 'rest');
        Screen('TextSize', window, 36);
        Screen('DrawLines', window, allCoords,...
            lineWidthPix, white, [xCenter yCenter], 2);
        Screen('Flip', window);
        PsychPortAudio('Start', file(1), 1, [], []); % Play metronome sound file (8 seconds)
        WaitSecs(t1+randi(t2))

        %Presentation of random letters on the screen during the finger
        %tapping test + recording of the key presses
    %     trig.beep(440, 0.2, 'finger_auto_dual');
        onset=GetSecs;

        %preallocate table with key presses
        keypresses=table('Size', [12, 3], 'VariableNames', {'onset', 'duration', 'value'}, 'VariableTypes', {'double', 'double', 'cell'});
        m=1; %first key press
        KbQueueFlush; % clear all previous key presses from the list


        %% CUEING
        cued=round(rand); %Cue=1 (true) Uncued=0 (false)
        if(cued==1) 
            %Start the Cue
            PsychPortAudio('Start', file(2), 1, [], []);
            outlet.push_sample(Marker_StartBlockCue1_5HzAddition);
            events_handautodual(ii).trial(j).cue='cued';
        else
            events_handautodual(ii).trial(j).cue='uncued';
        end

        %% LETTER PRESENTATION
        for n=1:N_letters
            %Present random letter on the screen
            Screen('TextSize', window, 100);
            DrawFormattedText(window, value{1}(n),'center','center', white);
            vbl = Screen('Flip', window);
            time_letter=rand(1)+0.5; %Speed with which the letters are presented

            %Meanwhile record key presses
            start_timer=GetSecs;
            while GetSecs-start_timer<time_letter
                [ pressed, firstPress, ~, lastPress, ~]=KbQueueCheck;
                if m<13 && pressed %not more than 12 keys can be saved
                    if isempty(find(firstPress~=lastPress)) % no key was pressed twice
                        keys=KbName(find(firstPress)); % find the pressed keys
                        [timing, idx]=sort(firstPress(find(firstPress))); % get timing of key presses in ascending order
                        if length(idx)>1
                            keys=keys(idx); % sort the pressed keys in ascending order
                        else
                            keys={keys};
                            key_n=length(keys); % number of pressed keys
                        end
                        for q=1:key_n
                            keypresses.onset(m)=timing(q); %store and record the timing
                            keyValue=keys(q);
                            % Get the numeric value of the response (clicking '2' leads to '2@')
                            keyValue=regexp(keyValue,'\d*','Match');
                            keypresses.value(m)=keyValue{:};%store and record the presses
                            m=m+1;
                            if m>12
                                break
                            end
                        end
                    else
                        error('key was pressed twice')
                    end
                end
            end

            %Between each letter show a red fixation cross
            Screen('DrawLines', window, allCoords,...
                lineWidthPix, [1 0 0], [xCenter yCenter], 2);
            Screen('Flip', window);
            WaitSecs (0.2);
        end

        %Present white fixation cross for some seconds to show that
        %trial is over
        duration=GetSecs-onset;
    %     trig.beep(440, 0.2, 'rest');
        Screen('TextSize', window, 36);
        Screen('DrawLines', window, allCoords,...
            lineWidthPix, white, [xCenter yCenter], 2);
        Screen('Flip', window);
        WaitSecs(5); % 5 seconds, so the nirs signal has time to go back to baseline

        %Stop cueing
        if (cued==1)
            PsychPortAudio('Stop', file(2));
        end

        %Ask how many G's were presented
        Screen('TextSize',window,30);
        DrawFormattedText(window, 'How many times was G presented? ','center','center', white);
        vbl = Screen('Flip', window);
        [secs, keyCode, deltaSecs]=KbWait;
        % Save the response and the key presses
        response={KbName(find(keyCode))}; 
        % Get the numeric value of the response (clicking '2' leads to '2@')
        response=regexp(response,'\d*','Match');
        response=response{:};
        events_handautodual(sequence_idx).trial(j).stimuli=table(onset,duration, value, response);
        events_handautodual(sequence_idx).trial(j).responses=keypresses;
        DrawFormattedText(window, ['Your answer: ' response{1} '\n Press any key to continue.'],'center','center', white);
        vbl = Screen('Flip', window);
        KbStrokeWait;

        if j<N_trials/2
            DrawFormattedText(window, 'Press any key to continue with the next trial. \n Note that you will first start with a fixation cross again. \n Start tapping the sequence as soon as a letter on the screen appears.' ,'center','center', white);
            vbl = Screen('Flip', window);
            KbStrokeWait;
        end

    end
    
    % Send end of block markers
    
    if sequence_idx==1
        outlet.push_sample(Marker_EndBlock_AutomaticSequence);
        sequence_idx=2; %set the next sequence
    else
        outlet.push_sample(Marker_EndBlock_NonAutomaticSequence);
        sequence_idx=1; %set the next sequence
    end
end



%After all trials completed, the end of the finger tapping task is
%reached
% Screen('TextSize',window,30);
% DrawFormattedText(window, 'This is the end of the automaticity test for the finger tapping task. \n You can take a rest if needed. \n When ready: press any key to end this session.' ,'center','center', white);
% vbl = Screen('Flip', window);
% save('events_handautodual.mat', 'events_handautodual'); % save the events
% KbStrokeWait; %wait for response to terminate instructions

%Show dual task performance in command window (finger tapping)
fprintf('%%%%%%%%%%%%%% Finger AutoDual %%%%%%%%%%%%%% \n')
for s=1:2
    fprintf('--- %s Sequence --- \n', events_handautodual(s).sequence_label)
    for h = 1:N_trials/2
        fprintf('Trial %d: \n', h)
        %Show if the answers for the number of G's presented were correct
        if str2num(events_handautodual(s).trial(h).stimuli.response{1})==length(strfind(events_handautodual(s).trial(h).stimuli.value{1}, 'G'))
            fprintf('G correct \n')
        else
            fprintf('G incorrect \n')
        end
        %Show if the tapping tempo was correct.
        margin=0.25; % margin of error: think about what is most convenient
        delay=mean(diff(events_handautodual(s).trial(h).responses.onset)-1/1.50);
        fprintf('the tempo was off with on average %f seconds \n', delay);
        %Show if the tapped sequence was correct
        correctSequence=sequences(s);
        if all(strcmp(events_handautodual(s).trial(h).responses.value,correctSequence{:}'))
            fprintf('Seq correct \n')
        else
            fprintf('Seq incorrect \n')
        end

    end
end

% End of automaticity test is reached 
Screen('TextSize',window,25);
DrawFormattedText(window,'You have completed the automaticity test. \n We will continue with the rest of the experiment. \n Press any key to end this session.','center', 'center', white);
vbl = Screen('Flip', window);
%Press key to end the session and return to the 'normal' screen.
KbStrokeWait;
sca

%% end the lsl session
% trig.pulseIR(3, 0.2); % stop trigger for the nirs recording
% delete(trig);
% ses.stop();
% diary off;

%% HELPER FUNCTIONS
function triglistener(src, event)
for ii=1:numel(event.Data)
  info=src.info;
  fprintf('   lsl event (%s) received @ %s with (uncorrected) timestamp %.3f \n',  event.Data{ii}, info.type, event.Timestamps(ii));
end
end

function cleanupFun()
delete(ses);
delete(trigstr{1});
delete(trigstr{2});
delete(info);
end

% To Play Back Sound
function [WAVstruct] = CreateWAVstruct(WAVfilename)
% This function creates a struct with the information from the wav-files.

    wav = WAVfilename;                                          
    WAVstruct = struct('wavedata',[],'fs',[],'nrChan',[]);      
    [WAVstruct.wavedata, WAVstruct.fs] = psychwavread(wav);     
    WAVstruct.wavedata = WAVstruct.wavedata';                   
    WAVstruct.nrChan = size(WAVstruct.wavedata,1);  
end