%% AUTOMATICITY TEST / DUAL TASK PARADIGM:
%A dual-task paradigm in which we test whether a 12 digit prelearned
% sequence (sequenceauto) has become an automatic movement for a finger
% tapping and a foot stomping task.

%% Settings:

% Clear the workspace and the screen
clear;

%Synch test skip => comment when actually testing patient
Screen('Preference', 'SkipSyncTests', 1);

root_dir='C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\Experiment\Combined-EEG-fNIRS-system';
addpath(genpath('C:\Users\joaop\Downloads\liblsl-Matlab'));
% addpath(genpath('C:\Users\catar\Downloads\liblsl-Matlab-master'));
% addpath(genpath('C:\Users\maria\OneDrive\Documentos\GitHub\liblsl-Matlab'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SET UP PARAMETERS

%The sequences used for this study (automatic and non-automatic sequences
%randomized between participants)

%Sequences used in order to be able to print in the command window if
%to generate a new sequence use:
% newseq1=randi([1 4], 1, 12)
% newseq2=randi([1 4], 1, 12)
% sequencesprint={ num2str(reshape(newseq1', 1, [])), num2str(reshape(newseq2', 1, [])) }
sequencesprint = {('4  1  4  2  1  2  3  2  2  4  3  3'),('4  2  4  4  2  3  1  1  3  4  4  1')};

sequences = {split(sequencesprint(1))',split(sequencesprint(2))'} ;

%Parameters for the resting period in between the trials
t1 = 20; %Resting period in seconds
t2 = 5;  %Random interval around the resting period time
t3 = 10; %Duration of a trial (tapping the sequence 1 time)

%Amount of letters presented during test for automaticity for one trial.
%Should be adjusted when letter presenting speed is changed!
N_letters=8; % 8 letters presented during a trial
N_trials=1; % number of trials per block

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
info    = lsl_streaminfo(lib, 'MJC_Automaticity', 'Markers', 1, 0.0, 'cf_int32', 'Automaticity_DualTask'); 

% Open an outlet for the data to run through.
outlet = lsl_outlet(info);

%% ASK TO START LABRECORDER
done=0;
while (done==0)
    correct=input('Please start up LabRecorder. After that, press Y to continue. \n', 's');
    if strcmpi(correct, 'y')
        done=1;
    end
end
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

Marker_CHECK = 1255;        % checkerboard flip
Marker_start = 1555;        % start signal 
Marker_stop = 1500;         % stop signal 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Open Pshychtoolbox.
PsychDefaultSetup(2);
KbName('UnifyKeyNames'); %Links the key presses to the key board names
KbQueueCreate;
KbQueueStart; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE FILES IN FOLDER

fprintf('Select the project directory \n')
root_dir=uigetdir('C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\Combined-EEG-fNIRS-system', 'Select the project directory');
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

% go to subject folder
sub_dir=fullfile(root_dir,'Data' ,sub);
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
%% LOAD METRONOME SOUNDS (PsychToolbox)
audio_dir=fullfile(root_dir, 'metronomesounds');
addpath(audio_dir)
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


% % Set Handle For Audio Capture (delay check)
% CAP_cue1_25Hz = PsychPortAudio('Open', [], 2, priority, Cue1_25Hz.fs, Cue1_25Hz.nrChan);

%AudioFile
file = [h_Metronome8; PPA_cue1_25Hz; h_Metronome300; h_Metronome600];
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
fixCrossDimPix = 50; % Here we set the size of the arms of our fixation cross
xCoords = [-fixCrossDimPix fixCrossDimPix 0 0]; % Set the coordinates (these are all relative to zero we will let the drawing routine center the cross in the center of our monitor for us)
yCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
allCoords = [xCoords; yCoords];
lineWidthPix = 10;% Set the line width for the fixation cross
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VIDEO PREPARATION
%Movie ID's to load
movie_id_automaticity=sort(randperm(40,20));
a=1:40;
movie_id_single=a(~ismember(a,movie_id_automaticity));

save('movie_id_automaticity.mat','movie_id_automaticity');
save('movie_id_single.mat','movie_id_single');

videodir=fullfile(root_dir,'LetterPresentation');
iter=1;
for ii=movie_id_automaticity
    videofilename=['LetterPresentation_',num2str(ii,'%d'),'.mov'];
    moviename=fullfile(videodir, videofilename);
%     moviename = [ PsychtoolboxRoot 'PsychDemos/MovieDemos/DualDiscs.mov' ];
    [id,duration]=Screen('OpenMovie', window, moviename);
    moviePtr.id(iter)=id;
    moviePtr.duration(iter)=duration;
    moviePtr.moviename(iter)={moviename};
    iter=iter+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PSEUDORANDOMIZATION
%Pseudorandomize which trials are cued and uncued (must be 50/50 split)
events_autodual=randCuedTrials(N_trials);

% Save the randomizations

str=['events_autodual_',sub,'_',rec];
save(str,'events_autodual');

%In case of a restart, if we want to reload the previous randomization:
% load(['events_autodual_',sub,'_',rec,'.mat']);
%% WELCOME SCREEN
HideCursor;
%Instruction automaticity test
outlet.push_sample(Marker_Test);
Screen('TextSize',window,30);
DrawFormattedText(window, 'Welcome to the experiment! \n \n If you have any questions please ask now. \n \n Thank you for participating! \n \n Press any key to continue','center', 'center', white);
Screen('Flip', window);
KbStrokeWait;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Instruction automaticity test

%Page 1
outlet.push_sample(Marker_Test);
Screen('TextSize',window,30);
DrawFormattedText(window,sprintf('You will start with the automaticity test in a dual task situation. \n \n You will be tested on the sequence you learned at home. \n\n After the actual experiment, you will also be tested on the sequence you learned today. \n\n For each, you will be tested %d times. \n\n Detailed instructions will be given at the start of each task. \n Press any key to continue.',N_trials),'center', 'center', white);
 Screen('Flip', window);
KbStrokeWait; %wait for response to terminate instructions


%Page 2
outlet.push_sample(Marker_Test);
Screen('TextSize',window,30);
DrawFormattedText(window, sprintf('While you perform the task, letters will be shown on the screen (D, G, Q, O). \n The goal is to perform the sequence tapping while counting how many times G is presented. \n\n After each time you tapped the full sequence, you should tell us how many times G was presented. \n\n For answering this question, \n keep in mind that when the answer is 4 you press 4 and not Return (Enter) on the keyboard. \n \n Press any key for the next instructions.'),'center','center', white);
 Screen('Flip', window);
KbStrokeWait; %wait for response to terminate instructions

%Page 3
outlet.push_sample(Marker_Test);
Screen('TextSize',window,30);
DrawFormattedText(window, sprintf('In between the trials you will see a fixation cross for 20 seconds. \n \n  You will hear a metronome sound during the first few seconds \n or during the entire 20 seconds. \n \n Tap the sequence on this rhythm, which is the same as you studied at home. \n \n After the fixation cross, the first trial will start automatically. \n So start tapping the sequence as soon as a letter on the screen appears. \n \n Note that during the tapping task you cannot talk. \n Try to keep your body movements as still as possible exept for the right hand. \n Keep your eyes open (also during the rest periods). \n Press any key to continue.'),'center','center', white);
 Screen('Flip', window);
KbStrokeWait; %wait for response to terminate instructions

%Page 4 - Demonstrate how cueing works
outlet.push_sample(Marker_Test);
Screen('TextSize',window,30);
DrawFormattedText(window,'While you perform the tasks and during the white crosses, \n you may hear a rhytmic sound like this one. \n Or you may just hear nothing. \n \n In any situation, you should start tapping when the letters appear','center', 'center', white);
 Screen('Flip', window);
PsychPortAudio('Start', file(2), 1, [], []);
WaitSecs(8);
PsychPortAudio('Stop', file(2));

%Page 5 - Show sequence
outlet.push_sample(Marker_Test);
Screen('TextSize',window,30);
DrawFormattedText(window, sprintf('Sequence: \n %s \n\n  When ready to start: press any key.', char(sequencesprint(1))),'center','center', white);
 Screen('Flip', window);
KbStrokeWait; %wait for response to terminate instructions

%% START TEST FOR AUTOMATICITY
% Start loop for the trials
for j=1:N_trials
    %Get letter list and letter frame stamp
    load(fullfile(videodir,['LetterPresentation_isLetterFrame_',num2str(movie_id_automaticity(j),'%d'),'.mat']));
    load(fullfile(videodir,['LetterPresentation_',num2str(movie_id_automaticity(j),'%d'),'.mat']));
    value=letters;
    
    %preallocate table with key presses for speed 
    keypresses=table('Size', [12, 3], 'VariableNames', {'onset', 'duration', 'value'}, 'VariableTypes', {'double', 'double', 'cell'});
    
    %Always start with a 20-25 seconds fixation cross with 8 seconds of metronome
    %sound
    Screen('TextSize', window, 30);
    Screen('DrawLines', window, allCoords,...
        lineWidthPix, white, [xCenter yCenter], 2);
    Screen('Flip', window);
    
    %% CUEING
    %If it is an uncued trial, play a metronome sound for 8 seconds and use
    %the rest of the rest time for baseline. If it is cued, just start the
    %cueing and stop it at the end of the trial
    cued=events_autodual.trial(j).cue; %Cue=1 (true) Uncued=0 (false)
    if cued==0 %uncued trial
        outlet.push_sample(Marker_StartBlock_Metronome);
        PsychPortAudio('Start', file(1), 1, [], []); % Play metronome sound file (8 seconds)
        WaitSecs(8);
        outlet.push_sample(Marker_EndBlock_Metronome);
        WaitSecs(t1-8+randi(t2)); %time that the white cross is shown
        outlet.push_sample(Marker_StartBlock_AutomaticSequence_Dual_Uncued);
    else %
        %Start the Cue
%         PsychPortAudio('GetAudioData',CAP_cue1_25Hz,120,[],[],[]);
%         startrecord=GetSecs;
%         PsychPortAudio('Start',CAP_cue1_25Hz,1, [], []);
%         WaitSecs(1)
        
        outlet.push_sample(Marker_StartBlock_Cue);
%         startmoment = GetSecs;
        startTime=PsychPortAudio('Start', file(2), 1, [], []);
        WaitSecs(t1+randi(t2)) %time that the white cross is shown
        outlet.push_sample(Marker_StartBlock_AutomaticSequence_Dual_Cued);
    end
    
    

    
    %% START VIDEO PRESENTATION
    %Presentation of random letters on the screen during the finger
    %tapping test + recording of the key presses
    onset=GetSecs;
    
    keypresses=playMovie(moviePtr.id(j),window, outlet, Marker_Keypress, isLetterFrame,keypresses);

    %% Stop cueing & push end trial markers
    if (cued==1)
%         audio1 = PsychPortAudio('GetAudioData',CAP_cue1_25Hz);
%         stoprecord = GetSecs;
        PsychPortAudio('Stop', file(2));
        outlet.push_sample(Marker_EndBlock_Cue);
        outlet.push_sample(Marker_EndBlock_AutomaticSequence_Dual_Cued);
    else
        outlet.push_sample(Marker_EndBlock_AutomaticSequence_Dual_Uncued);
    end
    
    %% Present white fixation cross for some seconds to show that
    %trial is over
    duration=GetSecs-onset;
    Screen('TextSize', window, 30);
    Screen('DrawLines', window, allCoords,...
        lineWidthPix, white, [xCenter yCenter], 2);
    Screen('Flip', window);
    WaitSecs((8-5).*rand(1) + 5); % 5 seconds, so the nirs signal has time to go back to baseline

    
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
%     keypresses = convertKeypresses(keypresses);
    keypresses=convertKeypresses_DEV(keypresses);
    moviename=moviePtr.moviename(j);
    events_autodual.trial(j).stimuli=table(onset,duration, value, response, moviename);
    events_autodual.trial(j).responses=keypresses;
    DrawFormattedText(window, ['Your answer: ' response{1} '\n Press any key to continue.'],'center','center', white);
     Screen('Flip', window);
    KbStrokeWait;

    %If it's in the last trial of the block (where we change the
    %sequence), prompt user to continue to next trial
    if j<N_trials/2
        DrawFormattedText(window, 'Press any key to continue with the next trial. \n Note that you will first start with a fixation cross again. \n Start tapping the sequence as soon as a letter on the screen appears.' ,'center','center', white);
         Screen('Flip', window);
        KbStrokeWait;
    end

end


%% Show dual task performance in command window for automatic sequence (finger tapping)
fprintf('%%%%%%%%%%%%%% Automatic Sequence Dual Task %%%%%%%%%%%%%% \n')
for h = 1:N_trials
    fprintf('Trial %d: \n', h)
    %Show if the answers for the number of G's presented were correct
    numberOfG=strfind(events_autodual.trial(h).stimuli.value, 'G');
    
    if ~isempty(events_autodual.trial(h).stimuli.response{1})
        if str2num(events_autodual.trial(h).stimuli.response{1})==length(numberOfG{1,1})
            fprintf('G correct \n')
        else
            fprintf('G incorrect \n')
        end
    else
        fprintf('G incorrect \n')
    end
    %Show if the tapping tempo was correct.
    margin=0.25; % margin of error: think about what is most convenient
    delay=mean(diff(events_autodual.trial(h).responses.onset)-1/1.50);
    fprintf('the tempo was off with on average %f seconds \n', delay);
    %Show if the tapped sequence was correct
    correctSequence=sequences(1);
    if all(strcmp(events_autodual.trial(h).responses.value,correctSequence{:}'))
        fprintf('Seq correct \n')
    else
        fprintf('Seq incorrect \n')
    end
end


%% End of automaticity test is reached 
Screen('TextSize',window,30);
DrawFormattedText(window,'You have completed the automaticity test. \n We will continue with the rest of the experiment. \n Press any key to end this session.','center', 'center', white);
 Screen('Flip', window);
%Press key to end the session and return to the 'normal' screen.
KbStrokeWait;
sca

%% Save Results
str=['events_autodual_',sub,'_',rec];
save(str,'events_autodual');


% save('cueRecording_Automaticity.mat','startmoment','startTime','startrecord','stoprecord','audio1');
%% HELPER FUNCTIONS
function events=randCuedTrials(n)
    %Generate a vector where half is 0's and half is 1's
    %1=cued, 0=uncued
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


%Convert keypresses of the experimental laptop to the expected numerical
%values
function keypresses = convertKeypresses(keypresses)
%Get the numeric value of the responses to keypressing as well
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
            otherwise %the subject clicked another numerical button besides the ones he was supposed to
                keypresses.value(h)={[]};
           end
        end
    end
end

%FOR DEVELOPMENTAL USE ONLY
%Created because for the PC where the code was created, pressing the
%numerical keys 1,2,3,4 gave the following responses 1!, 2", 3#, 4$, so we
%had to extract the numerical values
%Get the numeric value of the responses to keypressing as well
function keypresses=convertKeypresses_DEV(keypresses)
    for h=1:length(keypresses.value)
        try
            keyValue=regexp(keypresses.value(h),'\d*','Match'); 
            keypresses.value(h)=keyValue{:};
        catch %if regexp returns an error, the keypress does not have a numeric value
            keypresses.value(h)={[]};
        end
    end
end