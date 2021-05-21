%% AUTOMATIC SINGLE TASK

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

%% START OF THE TASK
Screen('TextSize',window,25);
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
        events_autosingle.trial(j).cue='uncued';
    else
        %Start the Cue
        PsychPortAudio('Start', file(2), 1, [], []);
        outlet.push_sample(Marker_StartBlock_Cue);
        events_autosingle.trial(j).cue='cued';   
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
    events_autosingle.trial(j).stimuli=table(onset,duration, value);
    events_autosingle.trial(j).responses=keypresses;
end

%% End of the automatic single task
Screen('TextSize',window,30);
DrawFormattedText(window,'This is the end of the automatic fingertapping task. \n \n Please rest as long as you want. \n\n Press any key to continue', 'center', 'center', white);
vbl = Screen('Flip', window);
KbStrokeWait;


%% Save results

str1=['events_autosingle',sub,'_',rec];
save(str1,'events_autosingle');

sca;