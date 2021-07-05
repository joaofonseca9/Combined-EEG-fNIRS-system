labrec='C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\Experiment\Data\Exp\source\sub-02\stim\sub-02_rec-02_triggers.xdf';
nirs_events='C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\Experiment\Data\Exp\pre-processed\sub-02\nirs\sub-02_rec-02_nirs_events.mat';
addpath('C:\Users\joaop\Downloads\xdf-Matlab-master');

load(nirs_events);
labrec_events=load_xdf(labrec);
labrec_events=labrec_events{1};

start_labrec=find(labrec_events.time_series==1700);
start_labrec=labrec_events.time_stamps(start_labrec(1));
keys=find(labrec_events.time_series==1777);
keys_time=labrec_events.time_stamps(keys);

nirs_samples=[nirs_events.sample];
start_nirs=find(strcmp({nirs_events.value}, sprintf('LSL 1700')));
start_nirs=nirs_samples(start_nirs(1));
keys_nirs=find(strcmp({nirs_events.value}, sprintf('LSL 1777')));
keys_time_nirs=nirs_samples(keys_nirs);

dif_first_key_labrec=keys_time(end)-start_labrec;
dif_first_key_nirs=(keys_time_nirs(end)-start_nirs)/50;

delay_labrecorder_nirs=abs(dif_first_key_labrec-dif_first_key_nirs)