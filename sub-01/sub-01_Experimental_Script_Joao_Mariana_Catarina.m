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

sequences = {split(sequencesprint(1))',split(sequencesprint(2))'} ;

%Randomize the first sequence to be tested
%1=auto sequence, 2 non-auto sequence
order_sequence=randperm(2,2);


%Parameters for the resting period in between the trials
t1 = 20; %Resting period in seconds
t2 = 5;  %Random interval around the resting period time
t3 = 9.5; %Duration of a trial (tapping the sequence 1 time)

%Amount of letters presented during test for automaticity for one trial.
%Should be adjusted when letter presenting speed is changed!
N_letters=8; % 8 letters presented during a trial
N_trials=1; % number of trials per block

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LSL SETUP
% LSL outlet sending events

% Create and load the lab streaming layer library


addpath(genpath('C:\Users\joaop\Downloads\liblsl-Matlab'));
% addpath(genpath('C:\Users\catar\Downloads\liblsl-Matlab-master'));
% addpath(genpath('C:\Users\maria\OneDrive\Documentos\GitHub\liblsl-Matlab'));
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
info    = lsl_streaminfo(lib, 'Dual Task', 'Markers', 1, 0.0, 'cf_int32', 'Automaticity_DualTask'); 

% Open an outlet for the data to run through.
outlet = lsl_outlet(info);

% MARKER SETUP
% Block related
instructions = 'instructions'; %NEVER USED (?)
finger_test='finger_test';
Marker_StartBlock_Cue       = 1700;         
Marker_EndBlock_Cue         = 1701;

Marker_StartBlock_AutomaticSequence     = 1702;
Marker_StartBlock_NonAutomaticSequence  = 1703;

Marker_StartBlock_AutomaticSequence_Dual     = 1704;
Marker_StartBlock_NonAutomaticSequence_Dual  = 1705;

Marker_EndBlock_AutomaticSequence       = 1712;
Marker_EndBlock_NonAutomaticSequence       = 1713;

Marker_EndBlock_AutomaticSequence_Dual      = 1714;
Marker_EndBlock_NonAutomaticSequence_Dual      = 1715;

Marker_CHECK = 1255;        % checkerboard flip
Marker_start = 1555;        % start signal 
Marker_stop = 1500;         % stop signal 

%Open Pshychtoolbox.
PsychDefaultSetup(2);
KbName('UnifyKeyNames'); %Links the key presses to the key board names
KbQueueCreate;
KbQueueStart; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD METRONOME SOUNDS (PsychToolbox)
audio_dir='.\metronomesounds';
cd(audio_dir)
[WAVMetronome8.wave,WAVMetronome8.fs]       = audioread('Metronome8.wav');
[WAVMetronome600.wave,WAVMetronome600.fs]       = audioread('Metronome600.wav');
[WAVMetronome300.wave,WAVMetronome300.fs]       = audioread('Metronome300.wav');

% change rows<>columns
WAVMetronome8.wave = WAVMetronome8.wave';         WAVMetronome8.nrChan=2;
WAVMetronome600.wave = WAVMetronome600.wave';         WAVMetronome600.nrChan=2;
WAVMetronome300.wave = WAVMetronome300.wave';         WAVMetronome300.nrChan=2;

% Get Cueing Files
Cue1_25Hz       = 'Metronome120.wav';
[Cue1_25Hz]     = CreateWAVstruct(Cue1_25Hz);
Cue1_25HzLength = length(Cue1_25Hz.wavedata)/Cue1_25Hz.fs;

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
h_device = PsychPortAudio ('GetDevices');

% Open handle
h_Metronome8   = PsychPortAudio('Open', [], [], priority, WAVMetronome8.fs, WAVMetronome8.nrChan);
h_Metronome600   = PsychPortAudio('Open', [], [], priority, WAVMetronome600.fs, WAVMetronome600.nrChan);
h_Metronome300   = PsychPortAudio('Open', [], [], priority, WAVMetronome300.fs, WAVMetronome300.nrChan);
PPA_cue1_25Hz = PsychPortAudio('Open', [], [], priority, Cue1_25Hz.fs, Cue1_25Hz.nrChan);

% Fill buffer
PsychPortAudio('FillBuffer', h_Metronome8, WAVMetronome8.wave);
PsychPortAudio('FillBuffer', h_Metronome600, WAVMetronome600.wave);
PsychPortAudio('FillBuffer', h_Metronome300, WAVMetronome300.wave);
PsychPortAudio('FillBuffer', PPA_cue1_25Hz, Cue1_25Hz.wavedata);

%AudioFile
file = [h_Metronome8; PPA_cue1_25Hz; h_Metronome300; h_Metronome600];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE FILES IN FOLDER

fprintf('Select the project directory \n')
root_dir=uigetdir('C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\Combined-EEG-fNIRS-system', 'Select the project directory');

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

% preparations for the fixation cross so we only need to do this once
fixCrossDimPix = 40; % Here we set the size of the arms of our fixation cross
xCoords = [-fixCrossDimPix fixCrossDimPix 0 0]; % Set the coordinates (these are all relative to zero we will let the drawing routine center the cross in the center of our monitor for us)
yCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
allCoords = [xCoords; yCoords];
lineWidthPix = 4;% Set the line width for the fixation cross
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WELCOME SCREEN
%Instruction automaticity test
Screen('TextSize',window,45);
DrawFormattedText(window, 'Welcome to the experiment! \n \n If you have any questions please ask now. \n \n Thank you for participating!','center', 'center', white);
vbl = Screen('Flip', window);
WaitSecs(10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Demonstrate how cueing works
Screen('TextSize',window,25);
DrawFormattedText(window,'While you perform the tasks, \n you may hear a rhytmic sound like this one. \n Or you may just hear nothing.','center', 'center', white);
vbl = Screen('Flip', window);
PsychPortAudio('Start', file(2), 1, [], []);
WaitSecs(5);
PsychPortAudio('Stop', file(2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% START TEST FOR AUTOMATICITY

%Empty structure for key presses, -> use later again so it saves the key
%presses within this structure -> save at the end
events_handautodual=struct([]); % same for the presented letters of the hand + the answer
events_handautodual(1).sequence_label='Automatic';
events_handautodual(2).sequence_label='Non-automatic';

%Instruction automaticity test
Screen('TextSize',window,25);
DrawFormattedText(window,sprintf('You will now start with the automaticity test in a dual task situation. \n \n You will be tested on the sequence you learned at home. \n\n After the actual experiment, you will also be tested on the sequence you learned today. \n\n For each, you will be tested %d times. \n\n Detailed instructions will be given at the start of each task. \n Press any key to continue.',N_trials),'center', 'center', white);
vbl = Screen('Flip', window);
KbStrokeWait; %wait for response to terminate instructions

Screen('TextSize',window,25);
DrawFormattedText(window, sprintf('While you perform the task, letters will be shown on the screen (A,G,O,L). \n The goal is to perform the sequence tapping while counting how many times G is presented. \n\n After each time you tapped the full sequence, you should tell us how many times G was presented. \n\n For answering this question, \n keep in mind that when the answer is 4 you press 4 and not Return (Enter) on the keyboard. \n\n Note that during the tapping task you cannot talk. \n Try to keep your body movements as still as possible except for the right hand. \n Keep your eyes open (also during the rest periods). \n\n In between the trials you will see a fixation cross for 20 seconds. \n During the first few seconds you will hear a metronome sound. \n Tap the sequence on this rhythm, which is the same as you studied at home. \n\n After the fixation cross, the first trial will start automatically. \n So start tapping the sequence as soon as a letter on the screen appears. \n When ready: press any key.'),'center','center', white);
vbl = Screen('Flip', window);
KbStrokeWait; %wait for response to terminate instructions


%Show sequence
Screen('TextSize',window,25);
DrawFormattedText(window, sprintf('Sequence: \n %s \n\n  When ready to start: press any key.', char(sequencesprint(1))),'center','center', white);
vbl = Screen('Flip', window);
KbStrokeWait; %wait for response to terminate instructions

%Start loop for the trials
for j=1:N_trials
    %Presentation of the letters on the screen (dual task). -> is random.
    %Participant has to count the times that G was presented.
    Letterlist='AGOL';
    letter_order=randi(length(Letterlist), 1, N_letters);
    value={Letterlist(letter_order)};

    endOfTrial=0; %helper variable to terminate cueing
    
    %Always start with a 20-25 seconds fixation cross with 8 seconds of metronome
    %sound
    Screen('TextSize', window, 36);
    Screen('DrawLines', window, allCoords,...
        lineWidthPix, white, [xCenter yCenter], 2);
    Screen('Flip', window);
    
    %% CUEING
    %If it is an uncued trial, play a metronome sound for 8 seconds and use
    %the rest of the rest time for baseline. If it is cued, just start the
    %cueing and stop it at the end of the trial
    cued=round(rand); %Cue=1 (true) Uncued=0 (false)
    if cued==0
        PsychPortAudio('Start', file(1), 1, [], []); % Play metronome sound file (8 seconds)
        events_handautodual(1).trial(j).cue='uncued';
    else
        %Start the Cue
        PsychPortAudio('Start', file(2), 1, [], []);
        outlet.push_sample(Marker_StartBlock_Cue);
        events_handautodual(1).trial(j).cue='cued';   
    end
    
    WaitSecs(t1+randi(t2)) %time that the white cross is shown
    
    %Presentation of random letters on the screen during the finger
    %tapping test + recording of the key presses
    outlet.push_sample(Marker_StartBlock_AutomaticSequence_Dual);
    onset=GetSecs;

    %preallocate table with key presses
    keypresses=table('Size', [12, 3], 'VariableNames', {'onset', 'duration', 'value'}, 'VariableTypes', {'double', 'double', 'cell'});
    m=1; %first key press
    KbQueueFlush; % clear all previous key presses from the list



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
    
    
    %Push end of trial Marker
    outlet.push_sample(Marker_EndBlock_AutomaticSequence_Dual);
    
    %Stop cueing
    if (cued==1)
        PsychPortAudio('Stop', file(2));
        outlet.push_sample(Marker_EndBlock_Cue);
    end
    
    %Present white fixation cross for some seconds to show that
    %trial is over
    duration=GetSecs-onset;
    Screen('TextSize', window, 36);
    Screen('DrawLines', window, allCoords,...
        lineWidthPix, white, [xCenter yCenter], 2);
    Screen('Flip', window);
    WaitSecs(5); % 5 seconds, so the nirs signal has time to go back to baseline

    
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
    events_handautodual(1).trial(j).stimuli=table(onset,duration, value, response);
    events_handautodual(1).trial(j).responses=keypresses;
    DrawFormattedText(window, ['Your answer: ' response{1} '\n Press any key to continue.'],'center','center', white);
    vbl = Screen('Flip', window);
    KbStrokeWait;

    %If it's in the last trial of the block (where we change the
    %sequence), prompt user to continue to next trial
    if j<N_trials/2
        DrawFormattedText(window, 'Press any key to continue with the next trial. \n Note that you will first start with a fixation cross again. \n Start tapping the sequence as soon as a letter on the screen appears.' ,'center','center', white);
        vbl = Screen('Flip', window);
        KbStrokeWait;
    end

end


%Show dual task performance in command window for automatic sequence (finger tapping)
fprintf('%%%%%%%%%%%%%% Automatic Sequence Dual Task %%%%%%%%%%%%%% \n')
for h = 1:N_trials
    fprintf('Trial %d: \n', h)
    %Show if the answers for the number of G's presented were correct
    if str2num(events_handautodual(1).trial(h).stimuli.response{1})==length(strfind(events_handautodual(1).trial(h).stimuli.value{1}, 'G'))
        fprintf('G correct \n')
    else
        fprintf('G incorrect \n')
    end
    %Show if the tapping tempo was correct.
    margin=0.25; % margin of error: think about what is most convenient
    delay=mean(diff(events_handautodual(1).trial(h).responses.onset)-1/1.50);
    fprintf('the tempo was off with on average %f seconds \n', delay);
    %Show if the tapped sequence was correct
    correctSequence=sequences(1);
    if all(strcmp(events_handautodual(1).trial(h).responses.value,correctSequence{:}'))
        fprintf('Seq correct \n')
    else
        fprintf('Seq incorrect \n')
    end
end


% End of automaticity test is reached 
Screen('TextSize',window,25);
DrawFormattedText(window,'You have completed the automaticity test. \n We will continue with the rest of the experiment. \n Press any key to end this session.','center', 'center', white);
vbl = Screen('Flip', window);
%Press key to end the session and return to the 'normal' screen.
KbStrokeWait;
% sca


%Save Results
str1=['events_handautodual_sub',sub_ID,'_',rec];
save(str1,'events_handautodual');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% START OF THE EXPERIMENT

%Empty structure for key presses -> use later again so it saves the key
%presses within this structure -> save at the end
events_handauto=struct([]); 
events_handauto(1).sequence_label='Automatic';
events_handauto(2).sequence_label='Non-automatic';

%Instruction experiment
Screen('TextSize',window,25);
DrawFormattedText(window,'You will now start with the finger tapping task (single task). \n \n You will either start with the automatic tasks, performing the at home studied sequence, \n or with the non-automatic tasks, the new sequence you learned today. \n Note that for the non-automatic tasks you will also perform an automaticity test (dual task) after the single task experiment. \n This is the same test as you just did for the automatic (at home studied) sequence. \n\n Detailed instructions will appear at the start of each new task. \n You can take a break in between tasks. \n This will be indicated in the on-screen instructions. \n\n Press any key to continue and see with which test you start.','center', 'center', white);
vbl = Screen('Flip', window);
KbStrokeWait; %wait for response to terminate instructions

%Start the randomization loop between non-automatic (=1) and automatic (=2)
for sequence_idx=order_sequence %Either [1,2] or [2,1] -> determines the order of the tasks
    if sequence_idx==2
        %% Practice new sequence
        Screen('TextSize', window, 25);
        DrawFormattedText(window, 'You will now perform the non-automatic finger tapping task. \n \n For the next 5 minutes you can practice a new sequence for the finger tapping task, \n the same way you practiced at home. \n After that we will start with the finger tapping task. \n Press any key to see the new sequence and start practicing.', 'center', 'center', white);
        vbl= Screen('Flip', window);
        KbStrokeWait; %wait for response to terminate instructions

        %Presenting the new (non-automatic) sequence on the screen
        Screen('TextSize', window, 50);
        DrawFormattedText(window, sprintf('%s', char(sequencesprint(sequence_idx))), 'center', 'center', white); % is this the same for each participant?
        vbl= Screen('Flip', window);
        %PsychPortAudio('Start', h_Metronome300, 1, [], []); %Play metronome sound file (5 minutes)
        WaitSecs(3);

        Screen('TextSize', window, 25);
        DrawFormattedText(window, 'The time to practice the new sequence is over. \n Press any key to continue to the finger tapping experiment.', 'center', 'center', white);
        vbl= Screen('Flip', window);
        WaitSecs(5)
        KbStrokeWait; %wait for response to terminate instructions

        %% NON-AUTOMATICITY TASKS
        DrawFormattedText(window, sprintf('As a reminder, the sequence you will perform is: \n %s \n\n  In between each trial there is a rest period of 20 seconds. \n During this rest you will hear a metronome sound, tap the sequence according to this interval sound. \n Trials and rest periods are indicated with red(= trial) and white(= rest) fixation crosses showing on the screen. \n \n When the red cross appears, please start tapping the sequence. \n When ready to start: press any key.', char(sequencesprint(sequence_idx))),'center','center', white);
        outlet.push_sample(Marker_StartBlock_NonAutomaticSequence);
        vbl = Screen('Flip', window);
        KbStrokeWait; %wait for response to terminate instructions
        
        %% Stimulus for finger tapping non-automatic sequence
        for j = 1:N_trials
            keypresses=table('Size', [12, 3], 'VariableNames', {'onset', 'duration', 'value'}, 'VariableTypes', {'double', 'double', 'cell'});
            %Rest period between each sequence 20-25 seconds

            % Fixation cross during rest
            Screen('TextSize', window, 36);        
            Screen('DrawLines', window, allCoords,...
            lineWidthPix, white, [xCenter yCenter], 2);
            Screen('Flip', window);
            %% CUEING
            %If it is an uncued trial, play a metronome sound for 8 seconds and use
            %the rest of the rest time for baseline. If it is cued, just start the
            %cueing and stop it at the end of the trial
            cued=round(rand); %Cue=1 (true) Uncued=0 (false)
            if cued==0
                PsychPortAudio('Start', file(1), 1, [], []); % Play metronome sound file (8 seconds)
                events_handautodual(1).trial(j).cue='uncued';
            else
                %Start the Cue
                PsychPortAudio('Start', file(2), 1, [], []);
                outlet.push_sample(Marker_StartBlock_Cue);
                events_handautodual(1).trial(j).cue='cued';   
            end

            WaitSecs(t1+randi(t2)) %time that the white cross is shown

            %% Red fixation cross during finger tapping trial
            outlet.push_sample(Marker_StartBlock_NonAutomaticSequence);
            onset=GetSecs;
            Screen('TextSize', window, 36);
            Screen('DrawLines', window, allCoords,...
            lineWidthPix, [1 0 0], [xCenter yCenter], 2);
            Screen('Flip', window)
            m=1;

            %% Record key presses
            start_timer=GetSecs;
            while GetSecs-start_timer < t3
            if m<13
                [secs, keyCode, deltaSecs] = KbWait([],2);
                if any(keyCode)
                    key={KbName(find(keyCode))};
                    % Get the numeric value of the response (clicking '2' leads to '2@')
                    keyValue=regexp(key,'\d*','Match');
                    keypresses.onset(m)=secs;
                    keypresses.value(m)=keyValue{:};%store and record the presses;
                    m=m+1;
                elseif KbName('ESCAPE')
                elseif GetSecs-start_timer >= t3
                    break
                end
            end
            end
            
            %% Stop cueing
            if (cued==1)
                PsychPortAudio('Stop', file(2));
                outlet.push_sample(Marker_EndBlock_Cue);
            end

            %% Short white fix cross after trial
            outlet.push_sample(Marker_EndBlock_NonAutomaticSequence);
            duration=GetSecs-onset;
            Screen('TextSize', window, 36);
            Screen('DrawLines', window, allCoords,...
            lineWidthPix, white, [xCenter yCenter], 2);
            Screen('Flip', window)
            WaitSecs (5)

            
            %% Show sequence before new trial starts
            %If it's in the last trial of the block (where we change the
            %sequence), prompt user to continue to next trial
            if j<N_trials
                DrawFormattedText(window, sprintf('Sequence:\n %s \n\n Note that you will first start with a fixation cross again. \n Start tapping the sequence as soon as the red cross appears.\n \n Press any key to continue with the next trial.' , char(sequencesprint(sequence_idx))),'center','center', white);
                vbl = Screen('Flip', window);
                KbStrokeWait;
            end

            % save the response and the key presses
            value={'red X'};
            events_handauto(sequence_idx).trial(j).stimuli=table(onset,duration, value);
            events_handauto(sequence_idx).trial(j).responses=keypresses;
        end
        
        %% AUTOMATICITY for the Non-automatic sequence (dual-tasking)
        %Instruction automaticity test
        Screen('TextSize',window,25);
        DrawFormattedText(window,sprintf('You will now start an automaticity test this new sequence, in a dual task situation,\n as you were for the prelearned sequence. \n \n You will be tested %d times. \n\n Detailed instructions will be given at the start of each task. \n Press any key to continue.',N_trials),'center', 'center', white);
        vbl = Screen('Flip', window);
        KbStrokeWait; %wait for response to terminate instructions
        
        Screen('TextSize',window,25);
        DrawFormattedText(window, sprintf('While you perform the task, letters will be shown on the screen (A,G,O,L). \n The goal is to perform the sequence tapping while counting how many times G is presented. \n After each time you tapped the full sequence, you should tell us how many times G was presented. \n For answering this question, \n keep in mind that when the answer is 4 you press 4 and not Return (Enter) on the keyboard. \n\n Note that during the tapping task you cannot talk. \n Try to keep your body movements as still as possible except for the right hand. \n Keep your eyes open (also during the rest periods). \n\n In between the trials you will see a fixation cross for 20 seconds. \n During the first few seconds you will hear a metronome sound. \n Tap the sequence on this rhythm, which is the same as you studied at home. \n\n After the fixation cross, the first trial will start automatically. \n So start tapping the sequence as soon as a letter on the screen appears. \n When ready: press any key.'),'center','center', white);
        vbl = Screen('Flip', window);
        KbStrokeWait; %wait for response to terminate instructions
        
        for t=1:N_trials
            %Presentation of the letters on the screen (dual task). -> is random.
            %Participant has to count the times that G was presented.
            Letterlist='AGOL';
            letter_order=randi(length(Letterlist), 1, N_letters);
            value={Letterlist(letter_order)};

            %Always start with a 20-25 seconds fixation cross with 8 seconds of metronome
            %sound
        %     trig.beep(440, 0.2, 'rest');
            Screen('TextSize', window, 36);
            Screen('DrawLines', window, allCoords,...
                lineWidthPix, white, [xCenter yCenter], 2);
            Screen('Flip', window);
            %% CUEING
            %If it is an uncued trial, play a metronome sound for 8 seconds and use
            %the rest of the rest time for baseline. If it is cued, just start the
            %cueing and stop it at the end of the trial
            cued=round(rand); %Cue=1 (true) Uncued=0 (false)
            if cued==0
                PsychPortAudio('Start', file(1), 1, [], []); % Play metronome sound file (8 seconds)
                events_handautodual(1).trial(j).cue='uncued';
            else
                %Start the Cue
                PsychPortAudio('Start', file(2), 1, [], []);
                outlet.push_sample(Marker_StartBlock_Cue);
                events_handautodual(1).trial(j).cue='cued';   
            end

            WaitSecs(t1+randi(t2)) %time that the white cross is shown

            %Presentation of random letters on the screen during the finger
            %tapping test + recording of the key presses
            outlet.push_sample(Marker_StartBlock_NonAutomaticSequence_Dual);
            onset=GetSecs;

            %preallocate table with key presses
            keypresses=table('Size', [12, 3], 'VariableNames', {'onset', 'duration', 'value'}, 'VariableTypes', {'double', 'double', 'cell'});
            m=1; %first key press
            KbQueueFlush; % clear all previous key presses from the list

            %% LETTER PRESENTATION
            for n=1:N_letters
                %Present random letter on the screen
                Screen('TextSize', window, 100);
                DrawFormattedText(window, value{1}(n),'center','center', white);
                vbl = Screen('Flip', window);
                time_letter=rand(1)+0.5; %Speed with which the letters are presented

                %% Meanwhile record key presses
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

                %% Between each letter show a red fixation cross
                Screen('DrawLines', window, allCoords,...
                    lineWidthPix, [1 0 0], [xCenter yCenter], 2);
                Screen('Flip', window);
                WaitSecs (0.2);
            end
            
            %% Set marker of end of dual task with non auto sequence
            
            outlet.push_sample(Marker_EndBlock_NonAutomaticSequence_Dual);
            
            %% Stop cueing
            if (cued==1)
                PsychPortAudio('Stop', file(2));
                outlet.push_sample(Marker_EndBlock_Cue);
            end

            %% Present white fixation cross for some seconds to show that
            %trial is over
            duration=GetSecs-onset;
        %     trig.beep(440, 0.2, 'rest');
            Screen('TextSize', window, 36);
            Screen('DrawLines', window, allCoords,...
                lineWidthPix, white, [xCenter yCenter], 2);
            Screen('Flip', window);
            WaitSecs(5); % 5 seconds, so the nirs signal has time to go back to baseline
            %% Ask how many G's were presented
            Screen('TextSize',window,30);
            DrawFormattedText(window, 'How many times was G presented? ','center','center', white);
            vbl = Screen('Flip', window);
            [secs, keyCode, deltaSecs]=KbWait;
            % Save the response and the key presses
            response={KbName(find(keyCode))}; 
            % Get the numeric value of the response (clicking '2' leads to '2@')
            response=regexp(response,'\d*','Match');
            response=response{:};
            events_handautodual(sequence_idx).trial(t).stimuli=table(onset,duration, value, response);
            events_handautodual(sequence_idx).trial(t).responses=keypresses;
            DrawFormattedText(window, ['Your answer: ' response{1} '\n Press any key to continue.'],'center','center', white);
            vbl = Screen('Flip', window);
            KbStrokeWait;

            %If it's in the last trial of the block (where we change the
            %sequence), prompt user to continue to next trial
            if t<N_trials
                DrawFormattedText(window, 'Press any key to continue with the next trial. \n Note that you will first start with a fixation cross again. \n Start tapping the sequence as soon as a letter on the screen appears.' ,'center','center', white);
                vbl = Screen('Flip', window);
                KbStrokeWait;
            end
        end
    end
    
    %% AUTOMATIC TASKS
    if sequence_idx==1
        DrawFormattedText(window, sprintf('You will now perform the sequence you learned at home for the finger tapping task: \n %s \n\n In between each trial there is a rest period of 20 seconds. \n During this rest you will hear a metronome sound, tap the sequence according to this interval sound. \n Trials and rest periods are indicated with red(= trial) and white(= rest) fixation crosses showing on the screen.\n \n When the red cross appears, please start tapping the sequence. \n When ready: press any key.', char(sequencesprint(sequence_idx))),'center','center', white);
        vbl = Screen('Flip', window);
        KbStrokeWait; %wait for response to terminate instructions
        
        %% Stimulus for finger tapping automatic sequence
        for j = 1:N_trials
            keypresses=table('Size', [12, 3], 'VariableNames', {'onset', 'duration', 'value'}, 'VariableTypes', {'double', 'double', 'cell'});
            %Rest period between each sequence 20-25 seconds

            % Fixation cross during rest
            Screen('TextSize', window, 36);        
            Screen('DrawLines', window, allCoords,...
            lineWidthPix, white, [xCenter yCenter], 2);
            Screen('Flip', window);
            %% CUEING
            %If it is an uncued trial, play a metronome sound for 8 seconds and use
            %the rest of the rest time for baseline. If it is cued, just start the
            %cueing and stop it at the end of the trial
            cued=round(rand); %Cue=1 (true) Uncued=0 (false)
            if cued==0
                PsychPortAudio('Start', file(1), 1, [], []); % Play metronome sound file (8 seconds)
                events_handautodual(1).trial(j).cue='uncued';
            else
                %Start the Cue
                PsychPortAudio('Start', file(2), 1, [], []);
                outlet.push_sample(Marker_StartBlock_Cue);
                events_handautodual(1).trial(j).cue='cued';   
            end

            WaitSecs(t1+randi(t2)) %time that the white cross is shown

            %Red fixation cross during finger tapping trial
            outlet.push_sample(Marker_StartBlock_AutomaticSequence);
            onset=GetSecs;
            Screen('TextSize', window, 36);
            Screen('DrawLines', window, allCoords,...
            lineWidthPix, [1 0 0], [xCenter yCenter], 2);
            Screen('Flip', window)
            m=1;

            %% Record key presses
            start_timer=GetSecs;
            while GetSecs-start_timer < t3
            if m<13
                [secs, keyCode, deltaSecs] = KbWait([],2);
                if any(keyCode)
                    key={KbName(find(keyCode))};
                    keypresses.onset(m)=secs;
                    keypresses.value(m)=key;
                    m=m+1;
                elseif KbName('ESCAPE')
                elseif GetSecs-start_timer >= t3
                    break
                end
            end
            end
            
            outlet.push_sample(Marker_EndBlock_AutomaticSequence);

            %% Short white fix cross after trial
            duration=GetSecs-onset;
            Screen('TextSize', window, 36);
            Screen('DrawLines', window, allCoords,...
            lineWidthPix, white, [xCenter yCenter], 2);
            Screen('Flip', window)
            WaitSecs (5)

            

            % Show sequence before new trial starts
            %If it's in the last trial of the block (where we change the
            %sequence), prompt user to continue to next trial
            if j<N_trials
                DrawFormattedText(window, sprintf('Sequence:\n %s \n Press any key to continue with the next trial. \n Note that you will first start with a fixation cross again. \n Start tapping the sequence as soon as a letter on the screen appears.' , char(sequencesprint(sequence_idx))),'center','center', white);
                vbl = Screen('Flip', window);
                KbStrokeWait;
            end

            % save the response and the key presses
            value={'red X'};
            events_handauto(sequence_idx).trial(j).stimuli=table(onset,duration, value);
            events_handauto(sequence_idx).trial(j).responses=keypresses;
        end
    end         
end
               
%% RESULTS    
%% Show performance in dual tasking
fprintf('%%%%%%%%%%%%%% finger Automatic/Non-Automatic Sequence Dual Task %%%%%%%%%%%%%% \n')
sequence_idx=2;
fprintf('--- %s Sequence --- \n', events_handautodual(sequence_idx).sequence_label)
for h = 1:N_trials
    fprintf('Trial %d: \n', h)
    %Show if the answers for the number of G's presented were correct
    if str2num(events_handautodual(sequence_idx).trial(h).stimuli.response{1})==length(strfind(events_handautodual(sequence_idx).trial(h).stimuli.value{1}, 'G'))
        fprintf('G correct \n')
    else
        fprintf('G incorrect \n')
    end
    %Show if the tapping tempo was correct.
    margin=0.25; % margin of error: think about what is most convenient
    delay=mean(diff(events_handautodual(sequence_idx).trial(h).responses.onset)-1/1.50);
    fprintf('the tempo was off with on average %f seconds \n', delay);
    %Show if the tapped sequence was correct
    correctSequence=sequences(sequence_idx);
    if all(strcmp(events_handautodual(sequence_idx).trial(h).responses.value,correctSequence{:}'))
        fprintf('Seq correct \n')
    else
        fprintf('Seq incorrect \n')
    end
end


%% Show single task performance in command window (finger tapping)
fprintf('%%%%%%%%%%%%%% finger Automatic/Non-Automatic Single Task %%%%%%%%%%%%%% \n')
for sequence_idx=order_sequence
    fprintf('--- %s Sequence --- \n', events_handauto(sequence_idx).sequence_label)
    for h = 1:N_trials
      %% Show tempo  
      fprintf('Trial %d: \n', h)
      margin=0.25; % margin of error: think about what is most convenient
      delay=mean(diff((events_handauto(sequence_idx).trial(h).responses.onset))-1/1.50);
      fprintf('the tempo was off with on average %f seconds \n', delay);
      
      %% Show if the tapped sequence was correct
      correctSequence=sequences(sequence_idx);
      if all(strcmp(events_handauto(sequence_idx).trial(h).responses.value,correctSequence{:}'))
        fprintf('Seq correct \n')
      else
          fprintf('Seq incorrect \n')
      end   
    end
end


%% End of the automatic/non-automatic single task
Screen('TextSize',window,30);
DrawFormattedText(window,'This is the end of the fingertapping tasks. \n \n We will now move on to the checkerboard task. \n\n Press any key to continue', 'center', 'center', white);
vbl = Screen('Flip', window);
KbStrokeWait;

%Save the order of the experiment
str2=['order_sub',sub_ID,'_',rec];
save(str2,'order_sequence');

%% Save results

str1=['events_handautodual_sub',sub_ID,'_',rec];
save(str1,'events_handautodual');
str2=['events_handauto_sub',sub_ID,'_',rec];
save(str2,'events_handauto');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CHECKERBOARD TASK

%% GENERATE CHECKERBOARD 

XYChecks = [24 18];                                             % nr of checks 4:3
nrChecks = XYChecks(1)*XYChecks(2);                             % total number of checks
dim = screenYpixels/XYChecks(2);                                      % dimension of a check [nr of pixels]
if dim*XYChecks(1)>screenXpixels
    error('ERROR: dimensions of checks are not correct');
end
baseRect = [0 0 dim*XYChecks(1) dim*XYChecks(2)];               % rect of dim*nrXpixels by dim pixels*nrYpixels
[Xpos,Ypos] = meshgrid(0:1:XYChecks(1)-1, 0:1:XYChecks(2)-1);   % Make coordinates for our grid of squares

% Reshape the matrices of coordinates into a vector (reshape) & 
% Scale grid spacing to the size of our squares and centre (.*dim)
Xgrid = reshape(Xpos, 1, nrChecks) .*dim;    
Ygrid = reshape(Ypos, 1, nrChecks) .*dim;  

% Visual angle
% > visual angle (V) = 2*arctan(size on screen (S) / 2* distance to screen(D))
Vangle = 1;                                                     % visual angle [degree]
dist = 65;                                                      % distance to screen [cm]
Lcheck = (2*dist)*tand(Vangle/2);                               % length of check on screen [cm]

% Set the colors of each of our squares
% > 1= white; 0=black
for iNrChecks = 1:XYChecks(2)
    if mod(iNrChecks,2)==0
        checkColor(iNrChecks,:) = repmat([1 0],1,XYChecks(1)/2); % [1 0 1 0 ...]
    else
        checkColor(iNrChecks,:) = repmat([0 1],1,XYChecks(1)/2); % [0 1 0 1 ...]
    end
end
checkColor = reshape(checkColor, 1, nrChecks);

% pixelnumbers for the checks 
% > [(xpos;ypos;length;heigth) x nrChecks] 
checks = [Xgrid ; Ygrid ; Xgrid+dim ; Ygrid+dim];

% Create Grid for the color with size of 'checks' (GridColor) &
% Overlap the GridColor on the Grid of the Checkerboard (=GridChecks)
for i=1:nrChecks
    GridColor = meshgrid(checks(1,i):1:checks(3,i)-1, checks(2,i):1:checks(4,i)-1) *0 +checkColor(1,i);
    GridChecks(checks(2,i)+1:1:checks(4,i), checks(1,i)+1:1:checks(3,i)) = GridColor;
    % imshow(GridChecks)
end

% Create checkerboard texture (1) and inverse checkerboard (2)
Texture(1)  = Screen('MakeTexture', window, GridChecks);
Texture(2)  = Screen('MakeTexture', window, 1-GridChecks);

% create red focus dot
dotColor = [1 0 0];                 % Set the color of our dot to full red (RBG)
dotSizePix = 20;                    % Dot size in pixels

% Draw the dot to the screen. 
% > Screen('DrawDots', windowPtr, xy [,size] [,color] [,center] [,dot_type]);
% > dot_type: round dots (1,2,3) > 2 tries to use high-quality anti-aliasing
Screen('DrawDots', window, [xCenter yCenter], dotSizePix, dotColor, [], 2);

%% GENERAL VARIABLES

% experimental variables
% > 'start' = time between keypress to start and first checkerboard presentation; 
% > grey screen with focus dot will be presented
startTime = 10;                             % time between pressing start and first checkerboard presentation (=grey screen) [sec]
startFrame = round(startTime / ifi);        % time between pressing start and first checkerboard presentation [nrSamples]
startFrameCount = 0;                        % counter for nr of start frames

% > 'flip' = flipping the contrast of the checkerboard presentations
flipTime = 0.45;                            % duration of checkerboard presentation, before flipping [s]
flipFrame = round(flipTime / ifi);          % duration of checkerboard presentation, before flipping [nrSamples]
flipNr = 300;                          % total nr of flips made (300)
flipNrCount = 0;                            % counter for nr of flips made
flipFrameCount = 0;                         % counting the number of frames                         
flipFrameWait = 0;                          % time to wait in frames for a flip

textureCue = [1 2];                         % cue that determines which texture (checkerboard) will be shown
SaveFrameLog = [];                          % save Screen presentation > (start=0; checkerboard = 1/2)

%% >> EXPERIMENT << %%
%%%%%%%%%%%%%%%%%%%%%%
%% Instructions
Screen('TextSize',window,70);
DrawFormattedText(window, 'The only thing you need to do for this task is look at the screen \n \n A red dot will indicate where you should look. \n \n Press any key to continue.','center','center', white);
vbl = Screen('Flip', window); 
KbStrokeWai



%% %>> START SCREEN <<%%%
% > grey screen with focus dot with text: 'Focus on red dot; press any key to START'
% > after pressing any key, the text removes
Screen('TextSize',window,70);
DrawFormattedText(window, 'Focus on the red dot \n\n\n\n Press any key to START','center','center', white);
vbl = Screen('Flip', window); 
Screen('DrawDots', window, [xCenter yCenter], dotSizePix, dotColor, [], 2);
KbStrokeWait; 
outlet.push_sample(Marker_start); % wait for key press and send LSL-marker (start)

% > after pressing any key, the text removes
% > 10sec grey screen with focus dot before checkerboards are presented 
while startFrameCount ~= startFrame
    vb1 = Screen('Flip', window);                   % Sync us to the vertical retrace
    startFrameCount = startFrameCount + 1;          % increment of counter
    SaveFrameLog = [SaveFrameLog; 0];               % save screenpresentation (startstimuli=0) for each frame
end

%%%>> CHECKERBOARD <<%%%
if startFrameCount == startFrameCount
    while flipNrCount ~= flipNr
        flipFrameCount = flipFrameCount + 1;            % increment of counter
        SaveFrameLog = [SaveFrameLog; textureCue(1)];   % save screenpresentation (checkerboard 1/2) for each frame

        % Draw our texture to the screen
        Screen('DrawTexture', window, Texture(textureCue(1)));
        Screen('DrawDots', window, [xCenter yCenter], dotSizePix, dotColor, [], 2);
        vbl = Screen('Flip', window, vbl + flipFrameWait); % ?? (flipFrameWait - 0.5) * ifi) / 0= default; flip on the next possible video retrace
              
        % Reverse texture cue to show inverse checkerboard if the time is up
        if flipFrameCount == flipFrame
            outlet.push_sample(Marker_CHECK);           % send LSL-marker (checkerboard flip)
            textureCue = fliplr(textureCue);            % flip index of texture
            flipNrCount = flipNrCount + 1;              % increment counter
            flipFrameCount = 0;                         % set frame counter on zero
        end
    end
end

%%%>> END SCREEN <<%%%
% > grey screen with focus dot with text: End of experiment; Press any key to EXIT'
if flipNrCount == flipNr
    Screen('DrawDots', window, [xCenter yCenter], dotSizePix, dotColor, [], 2);
    Screen('TextSize',window,70); % Stop screen
    DrawFormattedText(window, 'End of experiment. \n \n Thank you for participating! \n\n Press any key to EXIT','center','center', white);
    vbl = Screen('Flip', window); 
    outlet.push_sample(Marker_stop); 
end
KbStrokeWait; % wait for keypress
sca; % close screen

close all;
clear all;


%% HELPER FUNCTIONS
% To Play Back Sound
function [WAVstruct] = CreateWAVstruct(WAVfilename)
% This function creates a struct with the information from the wav-files.

    wav = WAVfilename;                                          
    WAVstruct = struct('wavedata',[],'fs',[],'nrChan',[]);      
    [WAVstruct.wavedata, WAVstruct.fs] = psychwavread(wav);     
    WAVstruct.wavedata = WAVstruct.wavedata';                   
    WAVstruct.nrChan = size(WAVstruct.wavedata,1);  
end
