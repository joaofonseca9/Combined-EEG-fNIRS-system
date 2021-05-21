%% NON-AUTOMATIC FINGERTAPPING TASKS
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

%% NON-AUTOMATIC FINGERTAPPING TASKS
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
        events_nonautosingle.trial(j).cue='uncued';
    else
        %Start the Cue
        PsychPortAudio('Start', file(2), 1, [], []);
        outlet.push_sample(Marker_StartBlock_Cue);
        events_nonautosingle.trial(j).cue='cued';   
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
            keypresses.value(m)=keyValue{:}; %store and record the presses;
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
    events_nonautosingle.trial(j).stimuli=table(onset,duration, value);
    events_nonautosingle.trial(j).responses=keypresses;
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
        events_nonautodual.trial(j).cue='uncued';
    else
        %Start the Cue
        PsychPortAudio('Start', file(2), 1, [], []);
        outlet.push_sample(Marker_StartBlock_Cue);
        events_nonautodual.trial(j).cue='cued';   
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
    events_nonautodual.trial(t).stimuli=table(onset,duration, value, response);
    events_nonautodual.trial(t).responses=keypresses;
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

%% Show dual task performance in command window for automatic sequence (finger tapping)
fprintf('%%%%%%%%%%%%%% Automatic Sequence Dual Task %%%%%%%%%%%%%% \n')
for h = 1:N_trials
    fprintf('Trial %d: \n', h)
    %Show if the answers for the number of G's presented were correct
    if str2num(events_nonautodual.trial(h).stimuli.response{1})==length(strfind(events_nonautodual.trial(h).stimuli.value{1}, 'G'))
        fprintf('G correct \n')
    else
        fprintf('G incorrect \n')
    end
    %Show if the tapping tempo was correct.
    margin=0.25; % margin of error: think about what is most convenient
    delay=mean(diff(events_nonautodual.trial(h).responses.onset)-1/1.50);
    fprintf('the tempo was off with on average %f seconds \n', delay);
    %Show if the tapped sequence was correct
    correctSequence=sequences(1);
    if all(strcmp(events_nonautodual.trial(h).responses.value,correctSequence{:}'))
        fprintf('Seq correct \n')
    else
        fprintf('Seq incorrect \n')
    end
end

%% Save results

str1=['events_nonautosingle',sub,'_',rec];
save(str1,'events_nonautosingle');

str1=['events_nonautodual',sub,'_',rec];
save(str1,'events_nonautodual');

sca;