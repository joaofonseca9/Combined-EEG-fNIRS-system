%% PLOTTING RECORDED CUES AND CHECKING LATENCY

%Get sound data
load('cueRecJoao.mat');
load('startmoment.mat');

load('stoprecord.mat');

load('startrecord.mat');

file=file.Cue1_25Hz;
fs=file.fs;

sound=file.wavedata;
sound = sound(1,onset:end);
sound = -sound;

time=(0:length(sound)-1)./fs;

%Generate cue signal
freq=1.25;
cues=zeros(size(sound));
cues(round((startmoment-startrecord)*fs):2*round(fs*1/freq):end) = max(max(sound))/2;

theoreticalcuestart=zeros(size(sound));
theoreticalcuestart(round((startmoment-startrecord)*fs))=max(max(sound))/2;

recordstart=zeros(size(sound));
recordstart(1)=max(max(sound))/2;

recordstop=zeros(size(sound));
recordstop(round((stoprecord-startrecord)*fs))=max(max(sound))/2;

stop=round((stoprecord-startrecord)*fs);
%Plotting

a=plot(time(1:stop), sound(1:stop),'b');
hold on
b=plot(time(1:stop), cues(1:stop),'r');
hold on

c=plot(time(1:stop),theoreticalcuestart(1:stop),'g');
hold on
d=plot(time(1:stop),recordstart(1:stop),'c');
hold on
e=plot(time(1:stop),recordstop(1:stop),'y');
legend([a(1), b(1),c(1),d(1),e(1)], 'Sound Recording', 'Theoretical cues','Theoretical cues Start', 'Start Record', 'Stop Record')
hold off