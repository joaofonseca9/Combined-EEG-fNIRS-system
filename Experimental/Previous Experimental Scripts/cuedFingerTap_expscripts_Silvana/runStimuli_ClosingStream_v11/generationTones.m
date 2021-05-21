%% Generate steady-state beeps
% > left and right ear 600 Hz


% clear all; close all; clc;
% % cd 'C:\1. PROMPT\A. Dry-EEG\Methods'
% 
% playSound = 2;
% saveSound = 1;


%% stimulus variables
function [beep_main,beep_sec,beep_sum]=generationTones(fs)
playSound = 0;
saveSound = 1;
% Amplitude of sinus in V
A_left  = 6.5;
A_right = 6.5;
% Frequency of sinus in Hz
f_left_main  = 146.83;
f_right_main = 146.83;
f_left_sec = 880;
f_right_sec = 880;


% Create WAV file
% fs = 44100; % sample frequency [Hz].
length_wav = 0.15; % length of the wav [s] The lenght of the beep/sound
t = 0:1/fs:length_wav-1/fs;

% prevent start/stop artifacts (fade in/out)
fade_length = 0.00005; % [sec]
fade_samples = round(fs*fade_length);   % nr of fade samples
fade_in = [0: (1/fade_samples) :fade_length].^(1/2);       % fade in from 0-1 in nr of fade_samples > flipped is fade_out
fade = ones(numel(t),numel(A_left));    % straight line (1) with length of signal
fade(1:numel(fade_in),:) = repmat(fade_in',1,numel(A_left)); % fade in
fade(size(fade,1)-numel(fade_in)+1:end,:) = repmat((fliplr(fade_in))',1,numel(A_left)); % fade out

% create signal
beep_left_main = fade .* (A_left.*sin(f_left_main*2*pi.*t'));        % fade * A*sin(f*2pi*t)
beep_right_main = fade .* (A_right.*sin(f_right_main*2*pi.*t'));     % fade * A*sin(f*2pi*t)
beep_main = [beep_left_main, beep_right_main];                % combine left/right frequencies in a matrix


beep_left_main_nonFade =(A_left.*sin(f_left_main*2*pi.*t'));
beep_right_main_nonFade =(A_right.*sin(f_right_main*2*pi.*t'));


beep_left_sec = fade .* (A_left.*sin(f_left_sec*2*pi.*t'));        % fade * A*sin(f*2pi*t)
beep_right_sec = fade .* (A_right.*sin(f_right_sec*2*pi.*t'));     % fade * A*sin(f*2pi*t)
beep_sec = [beep_left_sec, beep_right_sec]; % combine left/right frequencies in a matrix

beep_main_NonFade=[beep_left_main_nonFade,beep_right_main_nonFade];
% figure(1)
% plot(t,beep_main_NonFade)
% figure(2)
% plot(t,beep_main)
beep_sum=[beep_left_main+beep_left_sec,beep_right_sec];
%% double for stereo

if playSound == 1
    sound(beep_main);
%     sound(beep_sec);
    
end
% playSound=2;
if playSound ==2
    sound(beep_sec);
end
%% save as wav-file
if saveSound == 1
    audiowrite('Beep146Hz.wav',beep_main,fs)
    audiowrite('Beep880Hz.wav',beep_sec,fs)
%     audiowrite('BeepSum.wav',beep_sum,fs)
end
end 
