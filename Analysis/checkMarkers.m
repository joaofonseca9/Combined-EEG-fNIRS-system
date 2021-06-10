function [marker_table,eeg_length,nirs_length]=checkMarkers(EEG,nirs_raw, nirs_events)
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
%     eeg_no=sum(eeg_lbl==i);
    eeg_no=length(find(strcmp({EEG.event.type}, sprintf('s%d',i))));
    eeg_fnirs_no=length(find(strcmp({nirs_events.value}, sprintf('LSL %d',i))));
%     eeg_fnirs_no=sum(eeg_fnirs_lbl==i);
    marker_summary(1,n)=i;
    marker_summary(2,n)=eeg_no;
    marker_summary(3,n)=eeg_fnirs_no;
end

lbl = sprintfc('%d',lbl_numeric);
marker_table=array2table(marker_summary);
marker_table.Properties.VariableNames={'StartMetronome','StopMetronome','StartCue','StopCue','StartAutoCue','StartAutoNoCue','StartNonAutoCue','StartNonAutoNoCue','StartAutoDualCue','StartAutoDualNoCue','StartNonAutoDualCue','StartNonAutoDualNoCue','StopAutoCue','StopAutoNoCue','StopNonAutoCue','StopNonAutoNoCue','StopAutoDualCue','StopAutoDualNoCue','StopNonAutoDualCue','StopNonAutoDualNoCue','Key','CheckFlip','CheckStop','CheckStart','Letter','StartMovie','StopMovie'};


%% Check signal length
eeg_event_samp=[EEG.event.latency];
nirs_event_samp=[nirs_events.sample];

%Find first and last event (excluding test markers) 1600
nirs_event_idx=find(strncmp({nirs_events.value},'LSL',3));
test=find(strcmp({nirs_events.value},'LSL 1600'));
last_test=test(end);%sample of last test marker
first_marker=nirs_event_idx(nirs_event_idx==last_test)+1;

last_marker=find(strcmp({nirs_events.value},'LSL 1500'));
nirs_length=(nirs_event_samp(last_marker)-nirs_event_samp(first_marker))/nirs_raw.fsample;

eeg_event_idx=find(strncmp({EEG.event.type}, 's',1));
test=find(strcmp({EEG.event.type},'s1600'));
last_test=test(end);%sample of last test marker
first_marker=eeg_event_idx(eeg_event_idx==last_test)+1;

last_marker=find(strcmp({EEG.event.type},'s1500'));
eeg_length=(eeg_event_samp(last_marker)-eeg_event_samp(first_marker))/EEG.srate;

if nirs_length~=eeg_length
    warning('NIRS and EEG signals seem to have different lengths')
    fprintf('Time difference is %d' ,nirs_length-eeg_length)
end
 
% %Check time in hh:mm:ss
% nirs_length = seconds(nirs_length);  
% nirs_length.Format = 'hh:mm:ss'; 
% eeg_length = seconds(eeg_length);  
% eeg_length.Format = 'hh:mm:ss'; 
end
