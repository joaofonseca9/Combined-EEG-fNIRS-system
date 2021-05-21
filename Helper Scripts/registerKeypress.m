function keypresses=registerKeypress(outlet, marker, isLetterFrame,keypresses)

Marker_Letter = 1797;


m=1; %first key press
KbQueueFlush; % clear all previous key presses from the list
keyReleased=1;
n=1;
if(isLetterFrame(n)==1)
    outlet.push_sample(Marker_Letter);
end
n=n+1;

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
end