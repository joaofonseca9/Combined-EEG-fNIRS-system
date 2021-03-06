%% AUTOMATICITY TEST / DUAL TASK PARADIGM:
% A dual-task paradigm in which we test whether a 12 digit prelearned
% sequence (sequenceauto) has become an automatic movement for a finger
% tapping and a foot stomping task.

%% Settings:

% Clear the workspace and the screen
clear;

% Synch test skip => comment when actually testing patient
Screen('Preference', 'SkipSyncTests', 1);

root_dir='C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\Experiment\Combined-EEG-fNIRS-system';
addpath(genpath('C:\Users\joaop\Downloads\liblsl-Matlab'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SET UP PARAMETERS

% The sequences used for this study (automatic and non-automatic sequences
% randomized between participants)

% Sequences used in order to be able to print in the command window if
% to generate a new sequence use:
% newseq1=randi([1 4], 1, 12)
% newseq2=randi([1 4], 1, 12)
% sequencesprint={ num2str(reshape(newseq1', 1, [])), num2str(reshape(newseq2', 1, [])) }
% OUR SEQUENCES - script generateSequences
% 2 4 3 4 1 3 4 1 2 1 3 2
% 4 1 3 2 4 1 4 2 3 2 1 3
sequencesprint = {('2 4 3 4 1 3 4 1 2 1 3 2'),('4 1 3 2 4 1 4 2 3 2 1 3')};
sequences = {split(sequencesprint(1))',split(sequencesprint(2))'} ;

% Order of the sequences to be tested 
% 1=auto sequence, 2 non-auto sequence
order_sequence=[1,2];


% Parameters for the resting period in between the trials
t1 = 20; % Resting period in seconds
t2 = 5;  % Random interval around the resting period time
t3 = 12; % Duration of a trial (tapping the sequence 1 time)

% Amount of letters presented during test for automaticity for one trial.
% Should be adjusted when letter presenting speed is changed!
N_letters=15; % 8 letters presented during a trial
N_trials=2; % Number of trials per block

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LSL SETUP
% LSL outlet sending events

% Create and load the lab streaming layer library
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
info = lsl_streaminfo(lib, 'MJC_Automaticity', 'Markers', 1, 0.0, 'cf_int32', 'Automaticity_DualTask'); 

% Open an outlet for the data to run through.
outlet = lsl_outlet(info);


%% MARKER SETUP
% Block related
Marker_Test= 1600;

Marker_StartBlock_Metronome = 1698;
Marker_EndBlock_Metronome   = 1699;

Marker_StartBlock_Cue       = 1700;         
Marker_EndBlock_Cue         = 1701;

Marker_Keypress = 1777;
    
Marker_StartBlock_AutomaticSequence_Cued            = 1702;
Marker_StartBlock_AutomaticSequence_Uncued          = 1703;
Marker_StartBlock_NonAutomaticSequence_Cued         = 1704;
Marker_StartBlock_NonAutomaticSequence_Uncued       = 1705;

Marker_StartBlock_AutomaticSequence_Dual_Cued       = 1706;
Marker_StartBlock_AutomaticSequence_Dual_Uncued     = 1707;
Marker_StartBlock_NonAutomaticSequence_Dual_Cued    = 1708;
Marker_StartBlock_NonAutomaticSequence_Dual_Uncued  = 1709;

Marker_EndBlock_AutomaticSequence_Cued              = 1710;
Marker_EndBlock_AutomaticSequence_Uncued            = 1711;
Marker_EndBlock_NonAutomaticSequence_Cued           = 1712;
Marker_EndBlock_NonAutomaticSequence_Uncued         = 1713;

Marker_EndBlock_AutomaticSequence_Dual_Cued         = 1714;
Marker_EndBlock_AutomaticSequence_Dual_Uncued       = 1715;
Marker_EndBlock_NonAutomaticSequence_Dual_Cued      = 1716;
Marker_EndBlock_NonAutomaticSequence_Dual_Uncued    = 1717;

Marker_CHECK = 1255;        % Checkerboard flip
Marker_start = 1555;        % Start signal 
Marker_stop = 1500;         % Stop signal 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test Markers
done=0;
while (done==0)
    outlet.push_sample(Marker_Test);
    correct=input('Check if marker was sent \n (Y - continue; Other Key - Send another marker \n', 's');
    if strcmpi(correct, 'y')
        done=1;
    end
end

%% ASK TO START LABRECORDER
done=0;
while (done==0)
    correct=input('Please start up LabRecorder. After that, press Y to continue. \n', 's');
    if strcmpi(correct, 'y')
        done=1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test Markers again (After starting labrecorder)
done=0;
while (done==0)
    outlet.push_sample(Marker_Test);
    correct=input('Check if marker was sent \n (Y - continue; Other Key - Send another marker \n', 's');
    if strcmpi(correct, 'y')
        done=1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Open Pshychtoolbox.
PsychDefaultSetup(2);
KbName('UnifyKeyNames'); % Links the key presses to the key board names
KbQueueCreate;
KbQueueStart; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE FILES IN FOLDER

fprintf('Select the project directory, where the experimental scripts are \n')
% root_dir=uigetdir('C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\Combined-EEG-fNIRS-system', 'Select the project directory');
root_dir=uigetdir('C:\Users\maria\OneDrive\Documentos\GitHub\Combined-EEG-fNIRS-system', 'Select the project directory');
addpath(root_dir);
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

% Go to subject folder
%Folder Organization
%Combined-EEG -> Experimental (project directory)
%Data -> sub
sub_dir=fullfile(root_dir,'..','..','Data' ,sub);
if ~exist(sub_dir)
    mkdir(sub_dir)
end
cd(sub_dir)

logname=sprintf('%s_%s_logs.log', sub, rec); diary(logname);
% Save current script in subject directory
script=mfilename('fullpath');
script_name=mfilename;
copyfile(sprintf('%s.m', script), fullfile(sub_dir, sprintf('%s_%s.m', sub, script_name)))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD METRONOME SOUNDS (PsychToolbox)
audio_dir=fullfile(root_dir, 'metronomesounds');
addpath(audio_dir)
[WAVMetronome8.wave,WAVMetronome8.fs]       = audioread('Metronome8.wav');
[WAVMetronome300.wave,WAVMetronome300.fs]       = audioread('Metronome300.wav');

% Change rows<>columns
WAVMetronome8.wave = WAVMetronome8.wave';         WAVMetronome8.nrChan=2;
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
duration = 1;                       % Number of repetitions of the wav-file
PsychPortAudio('Verbosity',1);      % Verbosity = "wordiness" -> 1= print errors

% Get audio device
h_device = PsychPortAudio ('GetDevices');

% Open handle
h_Metronome8   = PsychPortAudio('Open', [], [], priority, WAVMetronome8.fs, WAVMetronome8.nrChan);
h_Metronome300   = PsychPortAudio('Open', [], [], priority, WAVMetronome300.fs, WAVMetronome300.nrChan);
PPA_cue1_25Hz = PsychPortAudio('Open', [], [], priority, Cue1_25Hz.fs, Cue1_25Hz.nrChan);

% Fill buffer
PsychPortAudio('FillBuffer', h_Metronome8, WAVMetronome8.wave);
PsychPortAudio('FillBuffer', h_Metronome300, WAVMetronome300.wave);
PsychPortAudio('FillBuffer', PPA_cue1_25Hz, Cue1_25Hz.wavedata);


% % Set Handle For Audio Capture (delay check)
% CAP_cue1_25Hz = PsychPortAudio('Open', [], 2, priority, Cue1_25Hz.fs, Cue1_25Hz.nrChan);

%AudioFile
file = [h_Metronome8; PPA_cue1_25Hz; h_Metronome300];
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
fixCrossDimPix = 50; % Here we set the size of the arms of our fixation cross
xCoords = [-fixCrossDimPix fixCrossDimPix 0 0]; % Set the coordinates (these are all relative to zero we will let the drawing routine center the cross in the center of our monitor for us)
yCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
allCoords = [xCoords; yCoords];
lineWidthPix = 10; % Set the line width for the fixation cross
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VIDEO PREPARATION
% Movie ID's to load
movie_id_auto=sort(randperm(40,20));
a=1:40;
movie_id_nonauto=a(~ismember(a,movie_id_auto));

videodir=fullfile(root_dir,'LetterPresentation');
iter=1;
for ii=movie_id_auto
    videofilename=['LetterPresentation_',num2str(ii,'%d'),'.mov'];
    moviename=fullfile(videodir, videofilename);
%     moviename = [ PsychtoolboxRoot 'PsychDemos/MovieDemos/DualDiscs.mov' ];
    [id,duration]=Screen('OpenMovie', window, moviename);
    moviePtr.auto.id(iter)=id;
    moviePtr.auto.duration(iter)=duration;
    moviePtr.auto.moviename(iter)=convertCharsToStrings(videofilename);
    iter=iter+1;
end
iter=1;
for ii=movie_id_nonauto
    videofilename=['LetterPresentation_',num2str(ii,'%d'),'.mov'];
    moviename=fullfile(videodir, videofilename);
    [id,duration]=Screen('OpenMovie', window, moviename);
    moviePtr.nonauto.id(iter)=id;
    moviePtr.nonauto.duration(iter)=duration;
    moviePtr.nonauto.moviename(iter)=convertCharsToStrings(videofilename);
    iter=iter+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PSEUDORANDOMIZATION

% Pseudorandomize which trials are cued and uncued (must be 50/50 split)
events_autodual=randCuedTrials(N_trials);
events_autosingle=randCuedTrials(N_trials);
events_nonautosingle=randCuedTrials(N_trials);
events_nonautodual=randCuedTrials(N_trials);


% Save randomizations 
str=['results_',sub,'_',rec];
save(str,'events_autodual','events_autosingle','events_nonautosingle','events_nonautodual');

% In case of a restart, if we want to reload the previous randomization:
% load(['results_',sub,'_',rec,'.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUTOMATICITY TEST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WELCOME SCREEN
HideCursor;
% Instruction automaticity test
Screen('TextSize',window,30);
DrawFormattedText(window, 'Welcome to the experiment! \n \n If you have any questions please ask now. \n \n Thank you for participating! \n \n Press any key to continue','center', 'center', white);
Screen('Flip', window);
KbStrokeWait;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Instruction automaticity test

% Page 1
Screen('TextSize',window,30);
DrawFormattedText(window,sprintf('You will start with the automaticity test in a dual task situation. \n \n You will be tested on the sequence you learned at home. \n\n After the actual experiment, you will also be tested on the sequence you will learn today. \n\n For each, you will be tested %d times. \n\n Detailed instructions will be given at the start of each task. \n Press any key to continue.',N_trials),'center', 'center', white);
 Screen('Flip', window);
KbStrokeWait; % Wait for response to terminate instructions


% Page 2
Screen('TextSize',window,30);
DrawFormattedText(window, sprintf('While you perform the task, letters will be shown on the screen (D, G, Q, O). \n The goal is to perform the sequence tapping while counting how many times G is presented. \n\n After each time you tapped the full sequence, you should tell us how many times G was presented. \n\n For answering this question, \n keep in mind that when the answer is 4 you press 4 and not Return (Enter) on the keyboard. \n \n Press any key for the next instructions.'),'center','center', white);
 Screen('Flip', window);
KbStrokeWait; % Wait for response to terminate instructions

%Page 3 - Demonstrate the letter sequence
Screen('TextSize',window,30);
DrawFormattedText(window,'You can see a demonstration of the mentioned letter presentations \n \n Press any key to see the letter presentation.','center', 'center', white);
Screen('Flip', window);
KbStrokeWait; % Wait for response to terminate instructions

videofilename='LetterPresentation_test.mov';
moviename=fullfile(videodir, videofilename);
[id]=Screen('OpenMovie', window, moviename);
playMovieTest(id,window);

% Page 4
Screen('TextSize',window,30);
DrawFormattedText(window, sprintf('In between the trials you will see a fixation cross for 20 seconds. \n \n  You will hear a metronome sound during the first few seconds \n or during the entire 20 seconds. \n \n Tap the sequence on this rhythm, which is the same as you studied at home. \n \n After the fixation cross, the first trial will start automatically. \n So start tapping the sequence as soon as a letter on the screen appears. \n \n Note that during the tapping task you cannot talk. \n Try to keep your body movements as still as possible exept for the right hand. \n Keep your eyes open (also during the rest periods). \n Press any key to continue.'),'center','center', white);
 Screen('Flip', window);
KbStrokeWait; % Wait for response to terminate instructions

% Page 5 - Demonstrate how cueing works
Screen('TextSize',window,30);
DrawFormattedText(window,'While you perform the tasks and during the white crosses, \n you may hear a rhytmic sound like this one. \n Or you may just hear nothing. \n \n In any situation, you should start tapping when the letters appear','center', 'center', white);
 Screen('Flip', window);
PsychPortAudio('Start', file(2), 1, [], []);
WaitSecs(8);
PsychPortAudio('Stop', file(2));

% Page 6 - Show sequence
Screen('TextSize',window,30);
DrawFormattedText(window, sprintf('Sequence: \n %s \n\n  When ready to start: press any key.', char(sequencesprint(1))),'center','center', white);
 Screen('Flip', window);
KbStrokeWait; % Wait for response to terminate instructions

%% START TEST FOR AUTOMATICITY
% Start loop for the trials
for j=1:N_trials
    % Get letter list and letter frame stamp
    load(fullfile(videodir,['LetterPresentation_isLetterFrame_',num2str(movie_id_auto(j),'%d'),'.mat']));
    load(fullfile(videodir,['LetterPresentation_',num2str(movie_id_auto(j),'%d'),'.mat']));
    value=letters;
    
    % Preallocate table with key presses for speed 
    keypresses=table('Size', [12, 3], 'VariableNames', {'onset', 'duration', 'value'}, 'VariableTypes', {'double', 'double', 'cell'});
    
    % Always start with a 20-25 seconds fixation cross with 8 seconds of metronome
    % sound
    Screen('TextSize', window, 30);
    Screen('DrawLines', window, allCoords,...
        lineWidthPix, white, [xCenter yCenter], 2);
    Screen('Flip', window);
    
    %% CUEING
    % If it is an uncued trial, play a metronome sound for 8 seconds and use
    % the rest of the rest time for baseline. If it is cued, just start the
    % cueing and stop it at the end of the trial
    cued=events_autodual.trial(j).cue; % Cue=1 (true) Uncued=0 (false)
    if cued==0 % Uncued trial
        outlet.push_sample(Marker_StartBlock_Metronome);
        PsychPortAudio('Start', file(1), 1, [], []); % Play metronome sound file (8 seconds)
        WaitSecs(8);
        outlet.push_sample(Marker_EndBlock_Metronome);
        WaitSecs(t1-8+randi(t2)); % Time that the white cross is shown
        outlet.push_sample(Marker_StartBlock_AutomaticSequence_Dual_Uncued);
    else % Cued trial
        % Start the Cue
%         PsychPortAudio('GetAudioData',CAP_cue1_25Hz,120,[],[],[]);
%         startrecord=GetSecs;
%         PsychPortAudio('Start',CAP_cue1_25Hz,1, [], []);
%         WaitSecs(1)
        
        outlet.push_sample(Marker_StartBlock_Cue);
%         startmoment = GetSecs;
        startTime=PsychPortAudio('Start', file(2), 1, [], []);
        WaitSecs(t1+randi(t2)) % Time that the white cross is shown
        outlet.push_sample(Marker_StartBlock_AutomaticSequence_Dual_Cued);
    end
    
    

    
    %% START VIDEO PRESENTATION
    % Presentation of random letters on the screen during the finger
    % tapping test + recording of the key presses
    onset=GetSecs;
    
    keypresses=playMovie(moviePtr.auto.id(j),window, outlet, isLetterFrame,keypresses,cued,'auto', file);
    
    duration=GetSecs-onset;

    %% Ask how many G's were presented
    Screen('TextSize',window,30);
    DrawFormattedText(window, 'How many times was G presented? ','center','center', white);
     Screen('Flip', window);
    [secs, keyCode, deltaSecs]=KbWait;
    
    % Save the response and the key presses
    response={KbName(find(keyCode))}; 

    
    % Get the numeric value of the response (clicking '2' leads to '2@')
    try
        response=regexp(response,'\d*','Match');
        response=response{:};
    catch 
        response={[]};
    end
    
    if isempty(response)
        response={[]};
    end
    %% Convert fingertapping responses to numerical values
    keypresses = convertKeypresses(keypresses);
%     keypresses=convertKeypresses_DEV(keypresses);
    moviename=moviePtr.auto.moviename(j);
    events_autodual.trial(j).stimuli=table(onset,duration, value, response, moviename);
    events_autodual.trial(j).responses=keypresses;
    DrawFormattedText(window, ['Your answer: ' response{1} '\n Press any key to continue.'],'center','center', white);
     Screen('Flip', window);
    KbStrokeWait;

    % If it's in the last trial of the block (where we change the
    % sequence), prompt user to continue to next trial
    if j<N_trials
        DrawFormattedText(window, 'Press any key to continue with the next trial. \n Note that you will first start with a fixation cross again. \n Start tapping the sequence as soon as a letter on the screen appears.' ,'center','center', white);
         Screen('Flip', window);
        KbStrokeWait;
    end

end

%% End of automaticity test is reached 
Screen('TextSize',window,30);
DrawFormattedText(window,'You have completed the automaticity test. \n The researcher will now take a look at your results. \n You can rest for as long as you need. \n Press any key to continue.','center', 'center', white);
Screen('Flip', window);
KbStrokeWait;


%% Show dual task performance in command window for automatic sequence (finger tapping)
Screen('TextSize',window,30);
for h = 1:N_trials
    trial=sprintf('Trial %d: \n \n',h);
    % Show if the answers for the number of G's presented were correct
    numberOfG=strfind(events_autodual.trial(h).stimuli.value, 'G');
    
    if ~isempty(events_autodual.trial(h).stimuli.response{1})
        if str2num(events_autodual.trial(h).stimuli.response{1})==length(numberOfG{1,1})
            g_count='G correct \n';
        else
            g_count='G incorrect \n';
        end
    else
        g_count='G incorrect \n';
    end
    % Show if the tapping tempo was correct.
    margin=0.25; % margin of error: think about what is most convenient
    delay=mean(diff(events_autodual.trial(h).responses.onset)-1/1.50);
    tempo=sprintf('The tempo was off with on average %f seconds \n', delay);
    % Show if the tapped sequence was correct
    correctSequence=sequences(1);
    if all(strcmp(events_autodual.trial(h).responses.value,correctSequence{:}'))
        seq='Sequence is correct \n \n';
    else
        seq='Sequence is incorrect \n \n';
    end
    
    DrawFormattedText(window,sprintf('%s %s %s %s Press any key to continue.',trial,g_count,tempo,seq), 'center', 'center', white);
    
    Screen('Flip', window);
    KbStrokeWait;
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SINGLE TASKING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HideCursor;
%% Instructions
% Page 1
Screen('TextSize',window,30);
DrawFormattedText(window,'You will now start with the finger tapping task. \n \n You will either start with the sequence you studied at home, \n or with a new sequence you will learn today for 5 minutes. \n \n Note that for this new sequence you will also perform an automaticity test. \n This is the same test as the one you just did for the \n automatic (at home studied) sequence. \n\n Detailed instructions will appear at the start of each new task. \n \n Press any key to continue and see with which test you start.','center', 'center', white);
vbl = Screen('Flip', window);
KbStrokeWait; % Wait for response to terminate instructions

% Start the randomization loop between non-automatic (=1) and automatic (=2)
for sequence_idx=order_sequence % Either [1,2] or [2,1] -> determines the order of the tasks
    if sequence_idx==2
        %% Practice new sequence
        Screen('TextSize', window, 30);
        DrawFormattedText(window, 'You will now perform the non-automatic finger tapping task. \n \n For the next 5 minutes you can practice a new sequence for the finger tapping task, \n the same way you practiced at home. \n In the first 3 minutes there is no cue \n and for the last 2 minutes of practice you will hear a cue. \n After that we will start with the finger tapping task. \n \n Press any key to see the new sequence and start practicing.', 'center', 'center', white);
        vbl= Screen('Flip', window);
        KbStrokeWait; % Wait for response to terminate instructions

        % Presenting the new (non-automatic) sequence on the screen
        Screen('TextSize', window, 50);
        DrawFormattedText(window, sprintf('%s', char(sequencesprint(sequence_idx))), 'center', 'center', white); % is this the same for each participant?
        vbl= Screen('Flip', window);
        WaitSecs(180); %180 
        PsychPortAudio('Start', file(2)); %Play metronome sound file (2 minutes)
        WaitSecs(120) %120
        PsychPortAudio('Stop', file(2));
        
        Screen('TextSize', window, 30);
        DrawFormattedText(window, 'The time to practice the new sequence is over. \n Press any key to continue to the finger tapping experiment.', 'center', 'center', white);
        vbl= Screen('Flip', window);
        KbStrokeWait; % Wait for response to terminate instructions

        %% NON-AUTOMATICITY TASKS
        DrawFormattedText(window, sprintf('As a reminder, the sequence you will perform is: \n %s \n\n  In between each trial there is a rest period of 20 seconds. \n During this rest you will hear a metronome sound for a few seconds \n or for the entire rest period. \n Tap the sequence according to this interval sound. \n \n Trials and rest periods are indicated with red(= trial) \n and white(= rest) fixation crosses showing on the screen. \n \n When the red cross appears, please start tapping the sequence. \n When ready to start: press any key.', char(sequencesprint(sequence_idx))),'center','center', white);
        vbl = Screen('Flip', window);
        KbStrokeWait; % Wait for response to terminate instructions
        
        %% Stimulus for finger tapping non-automatic sequence
        for j = 1:N_trials
            keypresses=table('Size', [12, 3], 'VariableNames', {'onset', 'duration', 'value'}, 'VariableTypes', {'double', 'double', 'cell'});
            % Rest period between each sequence 20-25 seconds

            % Fixation cross during rest
            Screen('TextSize',window,30);       
            Screen('DrawLines', window, allCoords,...
            lineWidthPix, white, [xCenter yCenter], 2);
            Screen('Flip', window);
            %% CUEING
            % If it is an uncued trial, play a metronome sound for 8 seconds and use
            % the rest of the rest time for baseline. If it is cued, just start the
            % cueing and stop it at the end of the trial
            cued=events_nonautosingle.trial(j).cue; % Cue=1 (true) Uncued=0 (false)
            if cued==0
                outlet.push_sample(Marker_StartBlock_Metronome);
                PsychPortAudio('Start', file(1), 1, [], []); % Play metronome sound file (8 seconds)
                WaitSecs(8);
                outlet.push_sample(Marker_EndBlock_Metronome);
                WaitSecs(t1-8+randi(t2)) % Time that the white cross is shown
                outlet.push_sample(Marker_StartBlock_NonAutomaticSequence_Uncued);
            else
                % Start the Cue
                PsychPortAudio('Start', file(2), 1, [], []);
                outlet.push_sample(Marker_StartBlock_Cue);
                WaitSecs(t1+randi(t2)) % Time that the white cross is shown
                outlet.push_sample(Marker_StartBlock_NonAutomaticSequence_Cued);
            end

            %% Red fixation cross during finger tapping trial
            onset=GetSecs;
            Screen('TextSize', window, 30);
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
                    outlet.push_sample(Marker_Keypress);
                    keypresses.value(m)={KbName(find(keyCode))};
                    keypresses.onset(m)=secs;
                    m=m+1;
                elseif KbName('ESCAPE')
                elseif GetSecs-start_timer >= t3
                    break
                end
            end
            end
            
            %% Stop cueing & push end trial markers
            if (cued==1)
                PsychPortAudio('Stop', file(2));
                outlet.push_sample(Marker_EndBlock_Cue);
                outlet.push_sample(Marker_EndBlock_NonAutomaticSequence_Cued);
            else
                outlet.push_sample(Marker_EndBlock_NonAutomaticSequence_Uncued);
            end

            %% Short white fix cross after trial
            duration=GetSecs-onset;
            Screen('TextSize', window, 30);
            Screen('DrawLines', window, allCoords,...
            lineWidthPix, white, [xCenter yCenter], 2);
            Screen('Flip', window)
            WaitSecs((8-5).*rand(1) + 5)

            
            %% Show sequence before new trial starts
            % If it's in the last trial of the block (where we change the
            % sequence), prompt user to continue to next trial
            if j<N_trials
                DrawFormattedText(window, sprintf('Sequence:\n %s \n\n Note that you will first start with a fixation cross again. \n Start tapping the sequence as soon as the red cross appears.\n \n Press any key to continue with the next trial.' , char(sequencesprint(sequence_idx))),'center','center', white);
                vbl = Screen('Flip', window);
                KbStrokeWait;
            end

            % Save the response and the key presses
            value={'red X'};
            %% Convert fingertapping responses to numerical values
            keypresses = convertKeypresses(keypresses);
%             keypresses=convertKeypresses_DEV(keypresses);
            events_nonautosingle.trial(j).stimuli=table(onset,duration, value);
            events_nonautosingle.trial(j).responses=keypresses;
        end
        
        %% AUTOMATICITY TEST for the Non-automatic sequence (dual-tasking)
        % Instruction automaticity test
        Screen('TextSize',window,30);
        DrawFormattedText(window,sprintf('You will now start an automaticity test for this new sequence, in a dual task situation,\n  like you did for the sequence you learned at home. \n \n You will be tested %d times. \n\n Detailed instructions will be given at the start of each task. \n Press any key to continue.',N_trials),'center', 'center', white);
        vbl = Screen('Flip', window);
        KbStrokeWait; % Wait for response to terminate instructions
        
        Screen('TextSize',window,30);
        DrawFormattedText(window, sprintf('While you perform the task, letters will be shown on the screen (D, G, Q, O). \n The goal is to perform the sequence tapping while counting how many times G is presented. \n After each time you tapped the full sequence, you should tell us how many times G was presented. \n For answering this question, \n keep in mind that when the answer is 4 you press 4 and not Return (Enter) on the keyboard. \n\n Note that during the tapping task you cannot talk. \n Try to keep your body movements as still as possible except for the right hand. \n Keep your eyes open (also during the rest periods). \n\n In between the trials you will see a fixation cross for 20 seconds. \n During the first few seconds you will hear a metronome sound. \n Tap the sequence on this rhythm, which is the same as you studied at home. \n\n After the fixation cross, the first trial will start automatically. \n So start tapping the sequence as soon as a letter on the screen appears. \n When ready: press any key.'),'center','center', white);
        vbl = Screen('Flip', window);
        KbStrokeWait; % Wait for response to terminate instructions
        
        for t=1:N_trials
            % Presentation of the letters on the screen (dual task). -> is random.
            % Participant has to count the times that G was presented.
            
            % Get letter list and letter frame stamps
            load(fullfile(videodir,['LetterPresentation_',num2str(movie_id_nonauto(t),'%d'),'.mat']));
            load(fullfile(videodir,['LetterPresentation_isLetterFrame_',num2str(movie_id_nonauto(t),'%d'),'.mat']));
            value=letters;
            
            % Preallocate table with key presses for speed 
            keypresses=table('Size', [12, 3], 'VariableNames', {'onset', 'duration', 'value'}, 'VariableTypes', {'double', 'double', 'cell'});
    
            % Always start with a 20-25 seconds fixation cross with 8 seconds of metronome
            % sound
        %     trig.beep(440, 0.2, 'rest');
            Screen('TextSize', window, 30);
            Screen('DrawLines', window, allCoords,...
                lineWidthPix, white, [xCenter yCenter], 2);
            Screen('Flip', window);
            %% CUEING
            % If it is an uncued trial, play a metronome sound for 8 seconds and use
            % the rest of the rest time for baseline. If it is cued, just start the
            % cueing and stop it at the end of the trial
            cued=events_nonautodual.trial(t).cue; % Cue=1 (true) Uncued=0 (false)
            if cued==0
                outlet.push_sample(Marker_StartBlock_Metronome);
                PsychPortAudio('Start', file(1), 1, [], []); % Play metronome sound file (8 seconds)
                WaitSecs(8);
                outlet.push_sample(Marker_EndBlock_Metronome);
                WaitSecs(t1-8+randi(t2)) % Time that the white cross is shown
                outlet.push_sample(Marker_StartBlock_NonAutomaticSequence_Dual_Uncued);
            else
                % Start the Cue
                PsychPortAudio('Start', file(2), 1, [], []);
                outlet.push_sample(Marker_StartBlock_Cue);  
                WaitSecs(t1+randi(t2)) % Time that the white cross is shown
                outlet.push_sample(Marker_StartBlock_NonAutomaticSequence_Dual_Cued);
            end

            

            %% START VIDEO PRESENTATION
            % Presentation of random letters on the screen during the finger
            % tapping test + recording of the key presses 
            onset=GetSecs;
            
            keypresses=playMovie(moviePtr.nonauto.id(t),window, outlet, isLetterFrame, keypresses, cued, 'nonauto', file);
            
            duration=GetSecs-onset;

            %% Ask how many G's were presented
            Screen('TextSize',window,30);
            DrawFormattedText(window, 'How many times was G presented? ','center','center', white);
            vbl = Screen('Flip', window);
            [secs, keyCode, deltaSecs]=KbWait;
            % Save the response and the key presses
            response={KbName(find(keyCode))}; 
            
            % Get the numeric value of the response (clicking '2' leads to '2@')
            try
                response=regexp(response,'\d*','Match');
                response=response{:};
            catch ME
                % If an error is spotted, like missclick, make that response
                % an empty cell
                response={[]};
            end
            
            if isempty(response)
                response={[]};
            end
            
            %% Convert fingertapping responses to numerical values
            keypresses = convertKeypresses(keypresses);
%             keypresses=convertKeypresses_DEV(keypresses);
            moviename=moviePtr.nonauto.moviename(t);
            events_nonautodual.trial(t).stimuli=table(onset,duration, value, response, moviename);
            events_nonautodual.trial(t).responses=keypresses;
            DrawFormattedText(window, ['Your answer: ' response{1} '\n Press any key to continue.'],'center','center', white);
            vbl = Screen('Flip', window);
            KbStrokeWait;

            % If it's in the last trial of the block (where we change the
            % sequence), prompt user to continue to next trial
            if t<N_trials
                DrawFormattedText(window, 'Press any key to continue with the next trial. \n Note that you will first start with a fixation cross again. \n Start tapping the sequence as soon as a letter on the screen appears.' ,'center','center', white);
                vbl = Screen('Flip', window);
                KbStrokeWait;
            end
        end
    end
    
    %% AUTOMATIC TASKS
    if sequence_idx==1
        DrawFormattedText(window, sprintf('You will now perform the sequence you learned at home for the finger tapping task: \n %s \n\n In between each trial there is a rest period of 20 seconds. \n During this rest you will hear a metronome sound during the first few seconds \n or during the entire rest period. \n Tap the sequence according to this interval sound. \n Trials and rest periods are indicated with red(= trial) \n and white(= rest) fixation crosses showing on the screen.\n \n When the red cross appears, please start tapping the sequence. \n When ready: press any key.', char(sequencesprint(sequence_idx))),'center','center', white);
        vbl = Screen('Flip', window);
        KbStrokeWait; % Wait for response to terminate instructions
        
        %% Stimulus for finger tapping automatic sequence
        for j = 1:N_trials
            keypresses=table('Size', [12, 3], 'VariableNames', {'onset', 'duration', 'value'}, 'VariableTypes', {'double', 'double', 'cell'});
            % Rest period between each sequence 20-25 seconds

            % Fixation cross during rest
            Screen('TextSize',window,30);       
            Screen('DrawLines', window, allCoords,...
            lineWidthPix, white, [xCenter yCenter], 2);
            Screen('Flip', window);
            %% CUEING
            % If it is an uncued trial, play a metronome sound for 8 seconds and use
            % the rest of the rest time for baseline. If it is cued, just start the
            % cueing and stop it at the end of the trial
            cued=events_autosingle.trial(j).cue; % Cue=1 (true) Uncued=0 (false)
            if cued==0
                outlet.push_sample(Marker_StartBlock_Metronome);
                PsychPortAudio('Start', file(1), 1, [], []); % Play metronome sound file (8 seconds)
                WaitSecs(8);
                outlet.push_sample(Marker_EndBlock_Metronome);
                WaitSecs(t1-8+randi(t2)) % Time that the white cross is shown
                outlet.push_sample(Marker_StartBlock_AutomaticSequence_Uncued);
            else
                % Start the Cue
                PsychPortAudio('Start', file(2), 1, [], []);
                outlet.push_sample(Marker_StartBlock_Cue);
                WaitSecs(t1+randi(t2)) % Time that the white cross is shown
                outlet.push_sample(Marker_StartBlock_AutomaticSequence_Cued);
            end

            % Red fixation cross during finger tapping trial
            onset=GetSecs;
            Screen('TextSize', window, 30);
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
                    outlet.push_sample(Marker_Keypress);
                    keypresses.value(m)={KbName(find(keyCode))};
                    keypresses.onset(m)=secs;
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
                outlet.push_sample(Marker_EndBlock_AutomaticSequence_Cued);
            else
                outlet.push_sample(Marker_EndBlock_AutomaticSequence_Uncued);
            end

            %% Short white fix cross after trial
            duration=GetSecs-onset;
            Screen('TextSize', window, 30);
            Screen('DrawLines', window, allCoords,...
            lineWidthPix, white, [xCenter yCenter], 2);
            Screen('Flip', window)
            WaitSecs((8-5).*rand(1) + 5)

            

            % Show sequence before new trial starts
            % If it's in the last trial of the block (where we change the
            % sequence), prompt user to continue to next trial
            if j<N_trials
                DrawFormattedText(window, sprintf('Sequence:\n %s \n Press any key to continue with the next trial. \n Note that you will first start with a fixation cross again. \n Start tapping the sequence as soon as a red cross appears on the screen.' , char(sequencesprint(sequence_idx))),'center','center', white);
                vbl = Screen('Flip', window);
                KbStrokeWait;
            end

         %% Convert fingertapping responses to numerical values and save results
            keypresses = convertKeypresses(keypresses);
%             keypresses=convertKeypresses_DEV(keypresses);
            value={'red X'};
            events_autosingle.trial(j).stimuli=table(onset,duration, value);
            events_autosingle.trial(j).responses=keypresses;
        end
    end         
end

%% End of the automatic/non-automatic single task
Screen('TextSize',window,30);
DrawFormattedText(window,'This is the end of the fingertapping tasks. \n \n We will now move on to the checkerboard task. \n\n Press any key to continue', 'center', 'center', white);
vbl = Screen('Flip', window);
KbStrokeWait;

%% RESULTS    
%% Show performance in dual tasking
fprintf('%%%%%%%%%%%%%% Finger Automatic/Non-Automatic Sequence Dual Task %%%%%%%%%%%%%% \n')
fprintf('--- Non-Automatic Sequence --- \n')
for h = 1:N_trials
    fprintf('Trial %d: \n', h)
    % Show if the answers for the number of G's presented were correct
    numberOfG=strfind(events_nonautodual.trial(h).stimuli.value, 'G');
     if ~isempty(events_nonautodual.trial(h).stimuli.response{1})
        if str2num(events_nonautodual.trial(h).stimuli.response{1})==length(numberOfG{1,1})
            fprintf('G correct \n')
        else
            fprintf('G incorrect \n')
        end
    else
        fprintf('G incorrect \n')
    end
    % Show if the tapping tempo was correct.
    margin=0.25; % Margin of error: think about what is most convenient
    delay=mean(diff(events_nonautodual.trial(h).responses.onset)-1/1.50);
    fprintf('the tempo was off with on average %f seconds \n', delay);
    % Show if the tapped sequence was correct
    correctSequence=sequences(sequence_idx);
    if all(strcmp(events_nonautodual.trial(h).responses.value,correctSequence{:}'))
        fprintf('Seq correct \n')
    else
        fprintf('Seq incorrect \n')
    end
end


%% Show single task performance in command window (finger tapping)
fprintf('%%%%%%%%%%%%%% Finger Automatic/Non-Automatic Single Task %%%%%%%%%%%%%% \n')
fprintf('--- Automatic Sequence --- \n')
for h = 1:N_trials
  %% Show tempo  
  fprintf('Trial %d: \n', h)
  margin=0.25; % Margin of error: think about what is most convenient
  delay=mean(diff((events_autosingle.trial(h).responses.onset))-1/1.50);
  fprintf('the tempo was off with on average %f seconds \n', delay);

  %% Show if the tapped sequence was correct
  correctSequence=sequences(1);
  if all(strcmp(events_autosingle.trial(h).responses.value,correctSequence{:}'))
    fprintf('Seq correct \n')
  else
      fprintf('Seq incorrect \n')
  end   
end

fprintf('--- Non Automatic Sequence --- \n')
for h = 1:N_trials
  %% Show tempo  
  fprintf('Trial %d: \n', h)
  margin=0.25; % Margin of error: think about what is most convenient
  delay=mean(diff((events_nonautosingle.trial(h).responses.onset))-1/1.50);
  fprintf('the tempo was off with on average %f seconds \n', delay);

  %% Show if the tapped sequence was correct
  correctSequence=sequences(2);
  if all(strcmp(events_nonautosingle.trial(h).responses.value,correctSequence{:}'))
    fprintf('Seq correct \n')
  else
      fprintf('Seq incorrect \n')
  end   
end



%% Save results - overwrite the structs, now with the results 
str=['results_',sub,'_',rec];
save(str,'events_autodual','events_autosingle','events_nonautosingle','events_nonautodual');

% save('cueRecording_Automaticity.mat','startmoment','startTime','startrecord','stoprecord','audio1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CHECKERBOARD TASK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Instructions
Screen('TextSize',window,30);
DrawFormattedText(window, 'Welcome to the checkerboard task. \n \n The only thing you need to do for this task is look at the screen \n \n A red dot will indicate where you should look. \n \n Press any key to continue.','center','center', white);
vbl = Screen('Flip', window); 
KbStrokeWait

%% GENERATE CHECKERBOARD 

XYChecks = [24 18];                                             % Nr of checks 4:3
nrChecks = XYChecks(1)*XYChecks(2);                             % Total number of checks
dim = screenYpixels/XYChecks(2);                                      % Dimension of a check [nr of pixels]
if dim*XYChecks(1)>screenXpixels
    error('ERROR: dimensions of checks are not correct');
end
baseRect = [0 0 dim*XYChecks(1) dim*XYChecks(2)];               % Rect of dim*nrXpixels by dim pixels*nrYpixels
[Xpos,Ypos] = meshgrid(0:1:XYChecks(1)-1, 0:1:XYChecks(2)-1);   % Make coordinates for our grid of squares

% Reshape the matrices of coordinates into a vector (reshape) & 
% Scale grid spacing to the size of our squares and centre (.*dim)
Xgrid = reshape(Xpos, 1, nrChecks) .*dim;    
Ygrid = reshape(Ypos, 1, nrChecks) .*dim;  

% Visual angle
% > visual angle (V) = 2*arctan(size on screen (S) / 2* distance to screen(D))
Vangle = 1;                                                     % Visual angle [degree]
dist = 65;                                                      % Distance to screen [cm]
Lcheck = (2*dist)*tand(Vangle/2);                               % Length of check on screen [cm]

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

% Pixelnumbers for the checks 
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

% Create red focus dot
dotColor = [1 0 0];                 % Set the color of our dot to full red (RBG)
dotSizePix = 20;                    % Dot size in pixels

% Draw the dot to the screen. 
% > Screen('DrawDots', windowPtr, xy [,size] [,color] [,center] [,dot_type]);
% > dot_type: round dots (1,2,3) > 2 tries to use high-quality anti-aliasing
Screen('DrawDots', window, [xCenter yCenter], dotSizePix, dotColor, [], 2);


%% GENERAL VARIABLES

% Experimental variables
% > 'start' = time between keypress to start and first checkerboard presentation; 
% > grey screen with focus dot will be presented
startTime = 10;                             % Time between pressing start and first checkerboard presentation (=grey screen) [sec]
startFrame = round(startTime / ifi);        % Time between pressing start and first checkerboard presentation [nrSamples]
startFrameCount = 0;                        % Counter for nr of start frames

% > 'flip' = flipping the contrast of the checkerboard presentations
flipTime = 0.45;                            % Duration of checkerboard presentation, before flipping [s]
flipFrame = round(flipTime / ifi);          % Duration of checkerboard presentation, before flipping [nrSamples]
flipNr = 300;                          % Total nr of flips made (300)
flipNrCount = 0;                            % Counter for nr of flips made
flipFrameCount = 0;                         % Counting the number of frames                         
flipFrameWait = 0;                          % Time to wait in frames for a flip

textureCue = [1 2];                         % Cue that determines which texture (checkerboard) will be shown
SaveFrameLog = [];                          % Save Screen presentation > (start=0; checkerboard = 1/2)

%% >> EXPERIMENT << %%
%%%%%%%%%%%%%%%%%%%%%%
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
    startFrameCount = startFrameCount + 1;          % Increment of counter
    SaveFrameLog = [SaveFrameLog; 0];               % Save screenpresentation (startstimuli=0) for each frame
end

%%%>> CHECKERBOARD <<%%%
if startFrameCount == startFrameCount
    while flipNrCount ~= flipNr
        flipFrameCount = flipFrameCount + 1;            % Increment of counter
        SaveFrameLog = [SaveFrameLog; textureCue(1)];   % Save screenpresentation (checkerboard 1/2) for each frame

        % Draw our texture to the screen
        Screen('DrawTexture', window, Texture(textureCue(1)));
        Screen('DrawDots', window, [xCenter yCenter], dotSizePix, dotColor, [], 2);
        vbl = Screen('Flip', window, vbl + flipFrameWait); % ?? (flipFrameWait - 0.5) * ifi) / 0= default; flip on the next possible video retrace
              
        % Reverse texture cue to show inverse checkerboard if the time is up
        if flipFrameCount == flipFrame
            outlet.push_sample(Marker_CHECK);           % Send LSL-marker (checkerboard flip)
            textureCue = fliplr(textureCue);            % Flip index of texture
            flipNrCount = flipNrCount + 1;              % Increment counter
            flipFrameCount = 0;                         % Set frame counter on zero
        end
    end
end

%%%>> END SCREEN <<%%%
% > grey screen with focus dot with text: End of experiment; Press any key to EXIT'
if flipNrCount == flipNr
%     Screen('DrawDots', window, [xCenter yCenter], dotSizePix, dotColor, [], 2);
    Screen('TextSize',window,30); % Stop screen
    DrawFormattedText(window, 'End of experiment. \n \n Thank you for participating! \n\n Press any key to EXIT','center','center', white);
    vbl = Screen('Flip', window); 
    outlet.push_sample(Marker_stop); 
end
KbStrokeWait; % Wait for keypress
sca; % Close screen

close all;
clear all;

%% HELPER FUNCTIONS
function events=randCuedTrials(n)
    % Generate a vector where half is 0's and half is 1's
    % 1=cued, 0=uncued
    isCued = zeros(n, 1);
    isCued(randperm(numel(isCued), round(n/2))) = 1;
    for j=1:n
        events.trial(j).cue=isCued(j);
    end
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


% Convert keypresses of the experimental laptop to the expected numerical
% values
function keypresses = convertKeypresses(keypresses)
% Get the numeric value of the responses to keypressing as well
    for h=1:length(keypresses.value)
        keyValue=cell2mat(keypresses.value(h));
        if isempty(keyValue)
            continue;
        else
           switch keyValue
            case '4'
                keypresses.value(h)={'1'};
            case '5'
                keypresses.value(h)={'2'};
            case '6'
                keypresses.value(h)={'3'};
            case 'BackSpace'  
                keypresses.value(h)={'4'};
            otherwise % The subject clicked another numerical button besides the ones he was supposed to
                keypresses.value(h)={[]};
           end
        end
    end
end

% FOR DEVELOPMENTAL USE ONLY
% Created because for the PC where the code was created, pressing the
% numerical keys 1,2,3,4 gave the following responses 1!, 2", 3#, 4$, so we
% had to extract the numerical values
% Get the numeric value of the responses to keypressing as well
function keypresses=convertKeypresses_DEV(keypresses)
    for h=1:length(keypresses.value)
        try
            keyValue=regexp(keypresses.value(h),'\d*','Match'); 
            keypresses.value(h)=keyValue{:};
        catch % If regexp returns an error, the keypress does not have a numeric value
            keypresses.value(h)={[]};
        end
    end
end