clc; close all; clear all;

load("cueRecording_Automaticity.mat");
fs = 48000;
% fs = 4*2^13;

onset = (startmoment-startrecord)*fs;
audio = audio1(1,onset:end);
audio = -audio;

t = (0:length(audio)-1)./fs;

cueonset = 0.022*fs;
% cueonset = 1;
cues = zeros(size(t));
cues(:,cueonset:(fs)/1.25:end) = 1;

figure
plot(t,audio)
% xlim([0 2])

% audio = audio(5*fs:25*fs);
% t = (0:length(audio)-1)./fs;
% cues = cues(5*fs:25*fs);

[pks,locs] = findpeaks(audio,'MinPeakHeight',.05,'MinPeakDistance',0.7*fs);
[pks2,locs2] = findpeaks(cues);

% pks(end) = []; locs(end) = [];
% pks2(1) = []; locs2(1) = [];

difference = (locs - locs2)/fs;
distance = diff(locs/fs);

mudif = mean(difference)
mudis = mean(distance)

sigdif = std(difference)
sigdis = std(distance)

figure
hold on
plot(t,audio,'b')
plot(locs/fs,pks,'*')
plot(t,cues)
plot(locs2/fs,pks2,'k*')
hold off
xlabel('seconden')
% xlim([0 3])
ylim([0 max(audio)])

[Pxx,Fxx] = periodogram(audio,[],[],48000);
figure
plot(Fxx,Pxx)

[b,a] = butter(4,[300 375]./(fs/2),'bandpass');

audiofilter = filtfilt(b,a,audio);

figure
hold on
plot(t,audio);
plot(t,cues);
hold off