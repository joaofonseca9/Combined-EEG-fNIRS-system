%% Comparison of the performance under the dual-tasking condition.

clear; clc; close all;
addpath('C:\Users\maria\OneDrive\Documentos\GitHub\Combined-EEG-fNIRS-system\Analysis');

laptop = 'laptopMariana';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);

subrec = ["02" "02"; "03" "02"; "04" "01"];
% Sub 02
% seq_auto = '434141243212';
% seq_nonauto = '212321324241';
% Sub 03
% seq_auto = '414212322433';
% seq_nonauto= '424423113441';
% Sub 04
seq_auto = '212321324241';
seq_nonauto = '413241423213';

autodual_finalAverageMistakes_cued = 0;
autodual_finalAverageMistakes_uncued = 0;
nonautodual_finalAverageMistakes_cued = 0;
nonautodual_finalAverageMistakes_uncued = 0;
autodual_finalIncorrectSequences_cued = 0;
autodual_finalIncorrectSequences_uncued = 0;
nonautodual_finalIncorrectSequences_cued = 0;
nonautodual_finalIncorrectSequences_uncued = 0;
autodual_finalDelay_cued = 0;
autodual_finalDelay_uncued = 0;
nonautodual_finalDelay_cued = 0;
nonautodual_finalDelay_uncued = 0;

% Go through every subject.
for subject = 3:3%1:size(subrec, 1)
    sub = subrec(subject, 1);
    rec = subrec(subject, 2);
    
    subvar = genvarname(sub);

    % Load the subject's results.
    load([mainpath_in, '\source\sub-', char(sub), '\stim\results_sub-',...
        char(sub), '_rec-', char(rec), '.mat']);
    
    %% Automatic sequence.
    % Check if counting answers from automatic sequence were correct.
    autodual_correct = checkCorrectCountingPerTrial(events_autodual);

    % Check average number of mistakes per condition (cued and uncued).
    [autodual_averageMistakes_cued, autodual_averageMistakes_uncued] = ...
        averageMistakesPerCondition(events_autodual, autodual_correct);
    
    % Check number of incorrectly typed sequences per number of trials.
    [autodual_incorrectSequences_cued, autodual_incorrectSequences_uncued] = ...
        checkIncorrectSequencePerTrial(events_autodual, seq_auto);
       
    % Get the delay in the fingertapping performance.
    [autodual_delay_cued, autodual_delay_uncued] = ...
        calculateDelayPerformance(events_autodual, seq_auto);
    
    % Add values to final average (all subjects).
    autodual_finalAverageMistakes_cued =...
        autodual_finalAverageMistakes_cued + autodual_averageMistakes_cued;
    autodual_finalAverageMistakes_uncued =...
        autodual_finalAverageMistakes_uncued + autodual_averageMistakes_uncued;
    autodual_finalIncorrectSequences_cued =...
        autodual_finalIncorrectSequences_cued + autodual_incorrectSequences_cued;
    autodual_finalIncorrectSequences_uncued =...
        autodual_finalIncorrectSequences_uncued + autodual_incorrectSequences_uncued;
    autodual_finalDelay_cued =...
        autodual_finalDelay_cued + autodual_delay_cued;
    autodual_finalDelay_uncued =...
        autodual_finalDelay_uncued + autodual_delay_uncued;
    
    %% Non-automatic sequence.
    % Check if counting answers from non-automatic sequence were correct.
    nonautodual_correct = checkCorrectCountingPerTrial(events_nonautodual);

    % Check average number of mistakes per condition (cued and uncued).
    [nonautodual_averageMistakes_cued, nonautodual_averageMistakes_uncued] = ...
        averageMistakesPerCondition(events_nonautodual, nonautodual_correct);
    
    % Check number of incorrectly typed sequences per number of trials.
    [nonautodual_incorrectSequences_cued, nonautodual_incorrectSequences_uncued] = ...
        checkIncorrectSequencePerTrial(events_nonautodual, seq_nonauto);
    
    % Get the delay in the fingertapping performance.
    [nonautodual_delay_cued, nonautodual_delay_uncued] = ...
        calculateDelayPerformance(events_nonautodual, seq_nonauto);
    
    % Add values to final average (all subjects).
    nonautodual_finalAverageMistakes_cued =...
        nonautodual_finalAverageMistakes_cued + nonautodual_averageMistakes_cued;
    nonautodual_finalAverageMistakes_uncued =...
        nonautodual_finalAverageMistakes_uncued + nonautodual_averageMistakes_uncued;
    nonautodual_finalIncorrectSequences_cued =...
        nonautodual_finalIncorrectSequences_cued + nonautodual_incorrectSequences_cued;
    nonautodual_finalIncorrectSequences_uncued =...
        nonautodual_finalIncorrectSequences_uncued + nonautodual_incorrectSequences_uncued;
    nonautodual_finalDelay_cued =...
        nonautodual_finalDelay_cued + nonautodual_delay_cued;
    nonautodual_finalDelay_uncued =...
        nonautodual_finalDelay_uncued + nonautodual_delay_uncued;
    
    %% Put values of error into final struct.
    % Values for the current subject.
    s.autodual_avgMistakes_cued = autodual_averageMistakes_cued;
    s.autodual_avgMistakes_uncued = autodual_averageMistakes_uncued;
    s.nonautodual_avgMistakes_cued = nonautodual_averageMistakes_cued;
    s.nonautodual_avgMistakes_uncued = nonautodual_averageMistakes_uncued;
    s.autodual_incorrectSeqs_cued = autodual_incorrectSequences_cued;
    s.autodual_incorrectSeqs_uncued = autodual_incorrectSequences_uncued;
    s.nonautodual_incorrectSeqs_cued = nonautodual_incorrectSequences_cued;
    s.nonautodual_incorrectSeqs_uncued = nonautodual_incorrectSequences_uncued;
    s.autodual_delay_cued = autodual_delay_cued;
    s.autodual_delay_uncued = autodual_delay_uncued;
    s.nonautodual_delay_cued = nonautodual_delay_cued;
    s.nonautodual_delay_uncued = nonautodual_delay_uncued;
    
    % Add struct of current subject to all subjects struct.
    allsubs.(genvarname(strcat('sub', char(sub)))) = s;
    
    % Add average values of all subjects to final struct.
    average.autodual_avgMistakes_cued = autodual_finalAverageMistakes_cued/size(subrec, 1);
    average.autodual_avgMistakes_uncued = autodual_finalAverageMistakes_uncued/size(subrec, 1);
    average.nonautodual_avgMistakes_cued = nonautodual_finalAverageMistakes_cued/size(subrec, 1);
    average.nonautodual_avgMistakes_uncued = nonautodual_finalAverageMistakes_uncued/size(subrec, 1);
    average.autodual_incorrectSeqs_cued = autodual_finalIncorrectSequences_cued/size(subrec, 1);
    average.autodual_incorrectSeqs_uncued = autodual_finalIncorrectSequences_uncued/size(subrec, 1);
    average.nonautodual_incorrectSeqs_cued = nonautodual_finalIncorrectSequences_cued/size(subrec, 1);
    average.nonautodual_incorrectSeqs_uncued = nonautodual_finalIncorrectSequences_uncued/size(subrec, 1);
    average.autodual_avgDelay_cued = autodual_finalDelay_cued/size(subrec, 1);
    average.autodual_avgDelay_uncued = autodual_finalDelay_uncued/size(subrec, 1);
    average.nonautodual_avgDelay_cued = nonautodual_finalDelay_cued/size(subrec, 1);
    average.nonautodual_avgDelay_uncued = nonautodual_finalDelay_uncued/size(subrec, 1);
    allsubs.avg = average;
   
end

% Save the struct from all subs.
save(strcat(pwd, '\allsubs.mat'), 'allsubs')

%% Functions necessary

function dual_correct = checkCorrectCountingPerTrial(events_dual)
% For each trial see if the counting response was equal to the actual
% number of G's.

dual = events_dual.trial;
dual_correct = zeros(1, length(dual));

for trial = 1:length(dual)
    stimuli = dual(trial).stimuli;
    if count(cell2mat(stimuli.value),'G') == str2num(cell2mat(stimuli.response))
        dual_correct(trial) = 1;
    end
end

end

function [dual_averageMistakes_cued, dual_averageMistakes_uncued] = ...
    averageMistakesPerCondition(events_dual, dual_correct)
% Check average number of mistakes per condition (cued and uncued).

dual = events_dual.trial;
dual_mistakes_cued = 0;
dual_mistakes_uncued = 0;

for trial = 1:length(dual)
    if dual_correct(trial) == 0 && dual(trial).cue == 1
        dual_mistakes_cued = dual_mistakes_cued + 1;
    elseif dual_correct(trial) == 0 && dual(trial).cue == 0
        dual_mistakes_uncued = dual_mistakes_uncued + 1;
    end
end
dual_averageMistakes_cued = dual_mistakes_cued/(length(dual)/2);
dual_averageMistakes_uncued = dual_mistakes_uncued/(length(dual)/2);

end

function [dual_incorrectSequences_cued, dual_incorrectSequences_uncued] = ...
    checkIncorrectSequencePerTrial(events_dual, real_sequence)

dual = events_dual.trial;
dual_incorrectSequences_cued = 0;
dual_incorrectSequences_uncued = 0;
    
for trial = 1:length(dual)
    typed_sequence = events_dual.trial(trial).responses.value;
    typed_sequence = typed_sequence(~cellfun('isempty', typed_sequence));
    typed_sequence = cell2mat(typed_sequence)';
    correctSequence = checkCorrectSequence(typed_sequence, real_sequence);
    if correctSequence == false && dual(trial).cue == 1
        dual_incorrectSequences_cued = dual_incorrectSequences_cued + 1;
    elseif correctSequence == false && dual(trial).cue == 0
        dual_incorrectSequences_uncued = dual_incorrectSequences_uncued + 1;
    end    
end
dual_incorrectSequences_cued = dual_incorrectSequences_cued/(length(dual)/2);
dual_incorrectSequences_uncued = dual_incorrectSequences_uncued/(length(dual)/2);

end

function [delay_cued, delay_uncued] = calculateDelayPerformance(events_dual, real_sequence)

dual = events_dual.trial;
delay_cued = 0;
delay_uncued = 0;

len_cued = length(dual)/2;
len_uncued = length(dual)/2;
    
for trial = 1:length(dual)
    typed_sequence = events_dual.trial(trial).responses.value;
    typed_sequence = typed_sequence(~cellfun('isempty', typed_sequence));
    typed_sequence = cell2mat(typed_sequence)';
    correctSequence = checkCorrectSequence(typed_sequence, real_sequence);

    if dual(trial).cue == 1 && correctSequence == true
        delay_cued = delay_cued + abs(mean(diff(dual(trial).responses.onset) - 1/1.50));
    elseif dual(trial).cue == 1 && correctSequence == false
        len_cued = len_cued - 1;
    elseif dual(trial).cue == 0 && correctSequence == true
        delay_uncued = delay_uncued + abs(mean(diff(dual(trial).responses.onset) - 1/1.50));
    elseif dual(trial).cue == 0 && correctSequence == false
        len_uncued = len_uncued - 1;
    end
    
end
    
delay_cued = delay_cued / len_cued;
delay_uncued = delay_uncued / len_uncued;
    
end

function correctSequence = checkCorrectSequence(typed_sequence, real_sequence)

if all(strcmp(typed_sequence, real_sequence))
    correctSequence = true;
else
    correctSequence = false;
end

end