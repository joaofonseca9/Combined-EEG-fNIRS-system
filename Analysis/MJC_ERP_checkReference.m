%% Initialize EEGLAB and paths´
laptop='laptopMariana';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);

eeglab;
ft_defaults;

sub='03';
rec='02';

file = getFileNames(mainpath_out, sub, rec);

%% Open the file after postICA
load(file.pstICA, 'EEG');

sub_path = fullfile(mainpath_in, 'incoming', ['sub-', sub]);
results = load(fullfile(sub_path, 'stim', ['results_sub-', sub,...
    '_rec-', rec]));
%%
marker_table = load('C:\Users\maria\Universidade do Porto\João Pedro Barbosa Fonseca - Internship\Experiment\Data\Pilots\pre-processed\sub-03\stim\sub-03_rec-02_marker_table').marker_table;


%% Set reference
% Re-reference the system to Cz 
[EEG_Cz] = pop_reref(EEG, 'Cz');

% Re-reference the system to M1
[EEG_M1] = pop_reref(EEG, 'M1');

%% Extract task data
% From system referenced to Cz 
[EEG_Cz_divided, file] = extractTaskData_EEG(EEG_Cz, marker_table,...
    results, file, mainpath_out);

% From system referenced to M1
[EEG_M1_divided, file] = extractTaskData_EEG(EEG_M1, marker_table,...
    results, file, mainpath_out);