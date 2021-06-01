function [marker_table]=checkMarkers(EEG,nirs_data_raw, nirs_events)
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
end
