%% PLOTTING RECORDED CUES AND CHECKING LATENCY
clear
%Get sound data
load('cueRecording.mat');

fs=48000;
onset = (startmoment-startrecord)*fs;
audio=audio1;
audio = audio(1,onset:end);
audio = -audio;

time=(0:length(audio)-1)./fs;

%Generate cue signal
freq=1.25;
cues=zeros(size(audio));
cues(round((startmoment-startrecord)*fs):round(fs*1/freq):end) = max(max(audio))/2;

theoreticalcuestart=zeros(size(audio));
theoreticalcuestart(round((startmoment-startrecord)*fs))=max(max(audio))/2;

recordstart=zeros(size(audio));
recordstart(1)=max(max(audio))/2;

recordstop=zeros(size(audio));
recordstop(round((stoprecord-startrecord)*fs))=max(max(audio))/2;

%Plotting

a=plot(time, audio,'b');
hold on
b=plot(time, cues,'r');
hold on

c=plot(time(1:stoprecord),theoreticalcuestart(1:stoprecord),'g');
hold on
d=plot(time(1:stoprecord),recordstart(1:stoprecord),'c');
hold on
e=plot(time(1:stoprecord),recordstop(1:stoprecord),'y');
legend([a(1), b(1),c(1),d(1),e(1)], 'Sound Recording', 'Theoretical cues','Theoretical cues Start', 'Start Record', 'Stop Record')
hold off