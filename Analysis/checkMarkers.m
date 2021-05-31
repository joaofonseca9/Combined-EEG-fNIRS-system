function [marker_table]=checkMarkers(EEG,nirs_data_raw, nirs_events)
%% Verify markers
% Get nirseeg numerical markers
eeg_fnirs_markers=zeros(1,nirs_data_raw.sampleinfo(2));
for i=1:length(nirs_events)
    if ~isempty(nirs_events(i).value)
        s=nirs_events(i).value;
        s=s(end-4:end);
        eeg_fnirs_markers(nirs_events(i).sample)=str2double(s);
    end
end

% Get eeg only numerical markers
eeg_events=EEG.event;
eeg_markers=zeros(1,EEG.pnts);
for i=2:length(EEG.event)
    a=EEG.event(i).type;
    eeg_markers(EEG.event(i).latency)=str2double(a(end-3:end));
end

nonzero=find(eeg_fnirs_markers~=0);
start_eeg_fnirs=nonzero(1);
stop_eeg_fnirs=nonzero(end);
eeg_fnirs_markers=eeg_fnirs_markers(start_eeg_fnirs:stop_eeg_fnirs);

nonzero=find(eeg_markers~=0);
start_eeg=nonzero(1);
stop_eeg=nonzero(end);
eeg_markers=eeg_markers(start_eeg:stop_eeg);

%% Get column of markers
tmp=struct2table(nirs_events);
value=tmp.value;
eeg_fnirs_markers_only=value(~cellfun('isempty',value));

tmp=struct2table(eeg_events);
eeg_markers_only=tmp.type;

%% Get labels 
eeg_fnirs_lbl=zeros(size(eeg_fnirs_markers_only));

for i=1:length(eeg_fnirs_markers_only)
    a=eeg_fnirs_markers_only(i);
    a=a{:};
    eeg_fnirs_lbl(i)=str2double(a(end-3:end));
end

eeg_lbl=zeros(size(eeg_markers_only));
for i=1:length(eeg_markers_only)
    a=eeg_markers_only(i);
    a=a{:};
    eeg_lbl(i)=str2double(a(end-3:end));
end

%% Get count of each marker
clear marker_summary
lbl=1698:1717;
lbl(end+1)=1777;
lbl(end+1)=1255;
lbl(end+1)=1500;
lbl(end+1)=1555;
lbl(end+1:end+3)=1797:1799;
lbl_numeric=lbl;
n=0;
for i=lbl_numeric
    n=n+1;
    eeg_no=sum(eeg_lbl==i);
    eeg_fnirs_no=sum(eeg_fnirs_lbl==i);
    marker_summary(1,n)=i;
    marker_summary(2,n)=eeg_no;
    marker_summary(3,n)=eeg_fnirs_no;
end

lbl = sprintfc('%d',lbl_numeric);
marker_table=array2table(marker_summary);
marker_table.Properties.VariableNames={'StartMetronome','StopMetronome','StartCue','StopCue','StartAutoCue','StartAutoNoCue','StartNonAutoCue','StartNonAutoNoCue','StartAutoDualCue','StartAutoDualNoCue','StartNonAutoDualCue','StartNonAutoDualNoCue','StopAutoCue','StopAutoNoCue','StopNonAutoCue','StopNonAutoNoCue','StopAutoDualCue','StopAutoDualNoCue','StopNonAutoDualCue','StopNonAutoDualNoCue','Key','CheckFlip','CheckStop','CheckStart','Letter','StartMovie','StopMovie'};
end
