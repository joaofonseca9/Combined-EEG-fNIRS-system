% MATLAB SCRIPT FOR THE AUTOMATICITY TEST
% A dual-task paradigm in which we test whether a 12 digit prelearned
% sequence has become an automatic movement for a finger tapping and a foot
% stomping task. The script is randomized in such a way that it either
% starts with the finger tapping or the foot stomping task. After that the
% experiment will automatically proceed for the other limb.

%Before starting the automaticity test, clear the workspace.
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LSL SETUP
% LSL outlet sending events
%
% Create and load the lab streaming layer library
lib = lsl_loadlib();
%
% Make a new stream outlet.
% info = lsl_streaminfo([lib handle],[name],[type],[channelcount],[fs],[channelformat],[sourceid])
% > name = name of stream; describes device/product
% > type = content type of stream (EEG, Markers)
% > channelcount = nr of channels per sample
% > fs = samplking rate (Hz) as advertized by data source
% > channelformat = cf_float32, cf__double64, cf_string, cf_int32, cf_int16
% > sourceid = unique identifier for source or device, if available
info = lsl_streaminfo(lib,'AutovsNAuto','Markers',1,0.0,'cf_string','sdfwerr32432');
%
% Open an outlet for the data to run through.
outlet = lsl_outlet(info);
%
% Create marker id's
instructions = 'instructions';
finger_test='finger_test';
foot_test='foot_test';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INITIALISATION

%Open Phsychtoolbox.
PsychDefaultSetup(2);
KbName('UnifyKeyNames'); %Links the key presses to the key board names
% ListenChar; % option A
KbQueueCreate; % option B
KbQueueStart; % option B

%Skip screen synchronization to prevent Pyshtoolbox for freezing
Screen('Preference', 'SkipSyncTests', 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LOAD METRONOME SOUNDS (PsychToolbox)
audio_dir='C:\Users\mtabo\Documents\TryOutScript\metronomesounds';
cd(audio_dir)
[WAVMetronome5.wave,WAVMetronome5.fs]       = audioread('Metronome5.wav');

% change rows<>columns
WAVMetronome5.wave = WAVMetronome5.wave';         WAVMetronome5.nrChan=2;

% CREATE AND FILL AUDIO BUFFER
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
h_Metronome5   = PsychPortAudio('Open', [], [], priority, WAVMetronome5.fs, WAVMetronome5.nrChan);

% Fill buffer
PsychPortAudio('FillBuffer', h_Metronome5, WAVMetronome5.wave);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SAVE FILES IN FOLDER

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RANDOMIZATION

%Amount of letters presented during test for automaticity for one trial.
%Should be adjusted when letter presenting speed is changed!
N_letters=8; % 8 letters presented during a trial, because a trial takes 7 seconds to perform
N_trials=8; % number of trials performed for each limb (8?)

%Create a vector to represent the two different options (1=finger tapping
%test, 2= foot stomping test).
order_test=[1,2];
%Randomize and determine the order of the two different tests. Either
%[1,2] or [2,1].
%order_test=stim_test(randperm(length(stim_test)));
%Save the order of the automaticity test experiment
save('order_test.mat', 'order_test');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCREEN PREPARATION

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
%START TEST FOR AUTOMATICITY

%Empty structure for key presses -> use later again so it saves the key
%presses within this structure -> save at the end
presses_handtest=struct([]);
letters_handtest=struct([]); % same for the presented letters of the hand + the answer
letters_foottest=struct([]); % same for the presented letters of the foot + the answer

%Instruction automaticity test
Screen('TextSize',window,25);
DrawFormattedText(window,'You will now start with the automaticity test. \n You will either start with the finger tapping or foot stomping task. \n Instructions will be given at the start of each task. \n Press any key to continue.','center', 'center', white);
vbl = Screen('Flip', window);
KbStrokeWait; %wait for response to terminate instructions

%Start the randomization loop
for i=order_test %Either [1,2] or [2,1] -> determines the order of the tasks
  
  % Finger tapping test -> 8 trials, presents letters upon randomized speed
  if i==1
    %Instruction automaticity task finger tapping
    outlet.push_sample({instructions})
    Screen('TextSize',window,25);
    DrawFormattedText(window, 'You will now perform the pre-learned sequence for the FINGER tapping task. \n  Letters will be shown on the screen (A,G,O,L) while you perform the task. \n The goal is to perform the sequence tapping while counting how many times G is presented. \n After each time you tapped the full sequence, you should tell us how many times G was presented. \n We will perform 8 trails. \n\n Note that during the tapping task you cannot talk. \n Try to keep your body movements as still as possible exept for the right hand. \n\n In between the trials you will see a fixation cross for 10 seconds. \n During the first 5 seconds you will hear a metronome sound. \n Tap the sequence on this rhythm, which is the same as you studied at home. \n\n We will start with a fixation cross on the screen for 10 seconds. \n After that the first trial will start automatically. \n So start tapping the sequence as soon as a letter on the screen appears. \n Press any key to continue and start the test.','center','center', white);
    vbl = Screen('Flip', window);
    KbStrokeWait; %wait for response to terminate instructions
    
    for j=1:N_trials
      %Presentation of the letters on the screen (dual task). -> is random.
      %Participant has to count the amount that G was presented.
      Letterlist= {'A', 'G', 'O', 'L'};
      letter_order=randi(length(Letterlist), 1, N_letters);
      letters_handtest(j).presented =Letterlist(letter_order);
      
      % always start with a 10 seconds fixation cross with 5 seconds of metronome
      % sound
      Screen('TextSize', window, 36);
      Screen('DrawLines', window, allCoords,...
        lineWidthPix, white, [xCenter yCenter], 2);
      Screen('Flip', window);
      PsychPortAudio('Start', h_Metronome5, 1, [], []); % Play metronome sound file
      WaitSecs(10);
      
      %Presentation of random letters on the screen during the finger
      %tapping test + recording of the key presses
      outlet.push_sample({finger_test})
      m=1; % first key press
      %           FlushEvents('keyDown'); % option A: clear all previous key presses from the list
      KbQueueFlush; % option B: clear all previous key presses from the list
      for n=1:N_letters
        % present random letter
        Screen('TextSize', window, 100);
        DrawFormattedText(window, [cell2mat(letters_handtest(j).presented(n))],'center','center', white);
        vbl = Screen('Flip', window);
        time_letter=rand(1)+0.5; %Speed with which the letters are presented = A randomized value between 0 and 1, + 0.5 sec
        
        % record key presses
        start_timer=GetSecs;
        while GetSecs-start_timer<time_letter
          %                   if CharAvail % option A
          %                     [ch, when]=GetChar;
          %                     fprintf('the key value was: %s \n', ch);
          %                     presses_handtest(j).key{m}=ch;
          %                     presses_handtest(j).secs(m)=when.secs;% better than using string arrays
          %                     m=m+1;
          %                   end
          [ pressed, firstPress, ~, lastPress, ~]=KbQueueCheck; % option B
          if pressed
            if isempty(find(firstPress~=lastPress)) % no key was pressed twice
              keys=KbName(find(firstPress)); % find the pressed keys
              [timing, idx]=sort(firstPress(find(firstPress))); % get timing of key presses in ascending order
              keys=keys(idx); % sort the pressed keys in ascending order
              key_n=length(keys); % number of pressed keys
              presses_handtest(j).key(m:m+key_n-1)=keys;
              presses_handtest(j).secs(m:m+key_n-1)=timing;
              m=m+key_n;
              KbQueueFlush;
            else
              error('key was pressed twice') % if this error occurs we need to find a way to handle this
            end
          end
        end
        
        % between each letter red fixation cross
        Screen('DrawLines', window, allCoords,...
          lineWidthPix, [1 0 0], [xCenter yCenter], 2);
        Screen('Flip', window);
        WaitSecs (0.2);
      end
      
      % present white fixation cross for some seconds to show that
      % trial is over
      Screen('TextSize', window, 36);
      Screen('DrawLines', window, allCoords,...
        lineWidthPix, white, [xCenter yCenter], 2);
      Screen('Flip', window);
      WaitSecs(2);
      
      % show feedback
      % ask how many G's were presented
      Screen('TextSize',window,30);
      DrawFormattedText(window, 'How many times was G presented? ','center','center', white);
      vbl = Screen('Flip', window);
      [secs, keyCode, deltaSecs]=KbWait;
      letters_handtest(j).reported_G=KbName(find(keyCode));
      DrawFormattedText(window, ['Your answer: ' letters_handtest(j).reported_G '\n Press any key to continue.'],'center','center', white);
      vbl = Screen('Flip', window);
      KbStrokeWait; %wait for response to terminate instruction
      % Show if the reported number of Gs is correct
      if str2num(letters_handtest(j).reported_G)==sum(strcmp(letters_handtest(j).presented, 'G'))
        DrawFormattedText(window, ['This is correct. \n Press any key to continue.' ],'center','center', white);
        vbl = Screen('Flip', window);
      else
        DrawFormattedText(window, ['This is incorrect. \n Press any key to continue.' ],'center','center', white);
        vbl = Screen('Flip', window);
      end
      KbStrokeWait;
      % Show if the tempo was correct. ! reflect if you want to show
      % this, because we cannot show this result for the foot stomping
      margin=0.1; % margin of error: think about what is most convenient
      if all(abs(diff(presses_handtest(j).secs)-1/1.5)<margin)
        DrawFormattedText(window, ['Your tempo was ok. \n Press any key to continue.' ],'center','center', white);
        vbl = Screen('Flip', window);
      else
        DrawFormattedText(window, ['Your tempo was not ok. \n Press any key to continue.' ],'center','center', white);
        vbl = Screen('Flip', window);
      end
      KbStrokeWait;
      DrawFormattedText(window, 'Press any key to continue with the next trail. \n Note that you will first start with a fixation cross again. \n Start tapping the sequence as soon as a letter on the screen appears.' ,'center','center', white);
      vbl = Screen('Flip', window);
      KbStrokeWait;
    end
    
    % After all trials completed, the end of the finger tapping task is
    % reached.
    Screen('TextSize',window,30);
    DrawFormattedText(window, 'This is the end of the automaticity test for the finger tapping task. \n  Press any key to end this session.' ,'center','center', white);
    vbl = Screen('Flip', window);
    save('letters_handtest.mat', 'letters_handtest'); % save the letters that were presented and the reported number of g's
    save('presses_handtest.mat', 'presses_handtest'); %Save which keys where pressed during the experiment
    KbStrokeWait; %wait for response to terminate instructions
    
  elseif i==2 % foot test, presents letters upon randomized speed
    % Instruction automaticity task foot stomping
    outlet.push_sample({instructions})
    Screen('TextSize',window,25);
    DrawFormattedText(window, 'You will now perform the pre-learned sequence for the FOOT stomping task. \n  Letters will be shown on the screen (A,G,O,L) while you perform the task. \n The goal is to perform the sequence stomping while counting how many times G is presented. \n After each time you stomped the full sequence, you should tell us how many times G was presented. \n We will perform 8 trials. \n\n Note that during the stomping task you cannot talk. \n Try to keep your body movements as still as possible exept for your right leg. \n\n In between the trials you will see a fixation cross for 10 seconds. \n During the first 5 seconds you will hear a metronome sound. \n Stomp the sequence on this rhythm, which is the same as you studied at home. \n\n We will start with a fixation cross on the screen for 10 seconds. \n After that the first trial will start automatically. \n So start stomping the sequence as soon as a letter on the screen appears. \n Press any key to continue and start the test.','center','center', white);
    vbl = Screen('Flip', window);
    KbStrokeWait; %wait for response to terminate instructions
    
    for j=1:N_trials
      %Presentation of the letters on the screen (dual task). -> is random.
      %Participant has to count the amount that G was presented.
      Letterlist= {'A', 'G', 'O', 'L'};
      letter_order=randi(length(Letterlist), 1, N_letters);
      letters_foottest(j).presented =Letterlist(letter_order);
      
      % always start with a fixation cross and 5 seconds of metronome
      % sound
      Screen('TextSize', window, 36);
      Screen('DrawLines', window, allCoords,...
        lineWidthPix, white, [xCenter yCenter], 2);
      Screen('Flip', window);
      PsychPortAudio('Start', h_Metronome5, 1, [], []); % Play metronome sound file
      WaitSecs(10);
      
      %Presentation of random letters on the screen during the foot
      %stomping test
      outlet.push_sample({foot_test})
      for n=1:N_letters
        % present random letter
        Screen('TextSize', window, 100);
        DrawFormattedText(window, [cell2mat(letters_foottest(j).presented(n))],'center','center', white);
        vbl = Screen('Flip', window);
        time_letter=rand(1)+0.5; %Speed with which the letters are presented = A randomized value between 0 and 1, + 0.5 sec
        WaitSecs(time_letter);
        
        % between each letter red fixation cross
        Screen('DrawLines', window, allCoords,...
          lineWidthPix, [1 0 0], [xCenter yCenter], 2);
        Screen('Flip', window);
        WaitSecs (0.2);
      end
      
      % present white fixation cross for some seconds to show that
      % trial is over
      Screen('TextSize', window, 36);
      Screen('DrawLines', window, allCoords,...
        lineWidthPix, white, [xCenter yCenter], 2);
      Screen('Flip', window);
      WaitSecs(2);
      
      % show feedback
      % ask how many G's were presented
      Screen('TextSize',window,30);
      DrawFormattedText(window, 'How many times was G presented? ','center','center', white);
      vbl = Screen('Flip', window);
      [secs, keyCode, deltaSecs]=KbWait;
      letters_foottest(j).reported_G=KbName(find(keyCode));
      DrawFormattedText(window, ['Your answer: ' letters_foottest(j).reported_G '\n Press any key to continue.'],'center','center', white);
      vbl = Screen('Flip', window);
      KbStrokeWait; %wait for response to terminate instruction
      % Show if the reported number of Gs is correct
      if str2num(letters_foottest(j).reported_G)==sum(strcmp(letters_foottest(j).presented, 'G'))
        DrawFormattedText(window, ['This is correct. \n Press any key to continue.' ],'center','center', white);
        vbl = Screen('Flip', window);
      else
        DrawFormattedText(window, ['This is incorrect. \n Press any key to continue.' ],'center','center', white);
        vbl = Screen('Flip', window);
      end
      KbStrokeWait;
      DrawFormattedText(window, 'Press any key to continue with the next trail. \n Note that you will first start with a fixation cross again. \n Start tapping the sequence as soon as a letter on the screen appears.' ,'center','center', white);
      vbl = Screen('Flip', window);
      KbStrokeWait;
    end
    
    % After all trials completed, the end of the foot stomping task is reached.
    Screen('TextSize',window,25);
    DrawFormattedText(window, 'End of the automaticity test for the foot stomping task. \n Press any key to end this session.','center','center', white);
    vbl = Screen('Flip', window);
    save('letters_foottest.mat', 'letters_foottest'); % save the letters that were presented and the reported number of g's
    KbStrokeWait; %wait for response to terminate instructions
  end
end

% End of automaticity test is reached (both limbs are tested)
Screen('TextSize',window,25);
DrawFormattedText(window,'You have completed the automaticity test. \n We will start with the preparation of the experiment now. \n Press any key to continue.','center', 'center', white);
vbl = Screen('Flip', window);
%Press key to end the session and return to the 'normal' screen.
KbStrokeWait; %wait for response to terminate instructions
sca