


%% EXP_AUDITORYSTEADYSTATEREPONSE

close all; clear all; clc;

stimtype = 'ASSR_AM40Hz_500CF.wav'; 
[SSassr]    = CreateWAVstruct (stimtype);

fprintf('\n\npress ENTER to start example stimuli\n'); KbStrokeWait;
sound(SSassr.wavedata,SSassr.fs)
fprintf('\n\npress ENTER to stop example stimuli\n'); KbStrokeWait;
clear sound




