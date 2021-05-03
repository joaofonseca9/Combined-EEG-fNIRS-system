function keypresses=playMovie(movie, win, outlet, marker, isLetterFrame)
% Most simplistic demo on how to play a movie.
%
% SimpleMovieDemo(moviename [, windowrect=[]]);
%
% This bare-bones demo plays a single movie whose name has to be provided -
% including the full filesystem path to the movie - exactly once, then
% exits. This is the most minimalistic way of doing it. For a more complex
% demo see PlayMoviesDemo. The remaining demos show more advanced concepts
% like proper timing etc.
%
% The demo will play our standard DualDiscs.mov movie if the 'moviename' is
% omitted.
%

% History:
% 02/05/2009  Created. (MK)
% 06/17/2013  Cleaned up. (MK)
% 25/04/2021 Changed. Joao Fonseca
Marker_StartMovie = 1798;
Marker_StopMovie = 1799;
Marker_Letter = 1797;

%preallocate table with key presses
keypresses=table('Size', [12, 3], 'VariableNames', {'onset', 'duration', 'value'}, 'VariableTypes', {'double', 'double', 'cell'});
m=1; %first key press
KbQueueFlush; % clear all previous key presses from the list
keyReleased=1;
n=1;
try
    % Start playback engine:
    Screen('PlayMovie', movie, 1);
    outlet.push_sample(Marker_StartMovie);
    % Playback loop: Runs until end of movie:
    while 1==1
        % Wait for next movie frame, retrieve texture handle to it
        [tex] = Screen('GetMovieImage', win, movie);
        % Valid texture returned? A negative value means end of movie reached:
        if tex<=0
            % We're done, break out of loop:
            break;
        end
        if(isLetterFrame(n)==1)
            outlet.push_sample(Marker_Letter);
        end
        n=n+1;
        
        % Draw the new texture immediately to screen:
        Screen('DrawTexture', win, tex);
%         frame_time=frame_idx-time_idx;
        
        [keyIsDown, secs, keyCode] = KbCheck;
        if keyIsDown==1 && keyReleased==1 && m<13  %only register the keypress once (wait for release of key to register again)
            outlet.push_sample(marker);
            keypresses.onset(m)=secs;
            keypresses.value(m)={KbName(find(keyCode))};
            keyReleased=0;
            m=m+1;
        elseif keyIsDown==0
            keyReleased=1;
        end
        % Update display:
        Screen('Flip', win);
        
        % Release texture:
        Screen('Close', tex);
        
        
    end
    
    % Stop playback:
    Screen('PlayMovie', movie, 0);
    outlet.push_sample(Marker_StopMovie);
    
    % Close movie:
    Screen('CloseMovie', movie);
    
    
catch %#ok<CTCH>
    sca;
    psychrethrow(psychlasterror);
end

return
