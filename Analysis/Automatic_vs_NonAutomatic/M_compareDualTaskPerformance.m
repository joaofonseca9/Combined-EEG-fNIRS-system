%% Comparison of the performance under the dual-tasking condition.

clear; clc; close all;
addpath('C:\Users\maria\OneDrive\Documentos\GitHub\Combined-EEG-fNIRS-system\Analysis');

laptop = 'laptopMariana';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);

subrec = ["02" "02"; "03" "02"];

autodual_finalAverageMistakes_cued = 0;
autodual_finalAverageMistakes_uncued = 0;
nonautodual_finalAverageMistakes_cued = 0;
nonautodual_finalAverageMistakes_uncued = 0;
autodual_finalDelay_cued = 0;
autodual_finalDelay_uncued = 0;
nonautodual_finalDelay_cued = 0;
nonautodual_finalDelay_uncued = 0;

% Go through every subject.
for subject = 1:size(subrec, 1)
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
    
    % Get the delay in the fingertapping performance.
    [autodual_delay_cued, autodual_delay_uncued] = ...
        calculateDelayPerformance(events_autodual);
    
    % Add values to final average (all subjects).
    autodual_finalAverageMistakes_cued =...
        autodual_finalAverageMistakes_cued + autodual_averageMistakes_cued;
    autodual_finalAverageMistakes_uncued =...
        autodual_finalAverageMistakes_uncued + autodual_averageMistakes_uncued;
    autodual_finalDelay_cued =...
        autodual_finalDelay_cued + abs(autodual_delay_cued);
    autodual_finalDelay_uncued =...
        autodual_finalDelay_uncued + abs(autodual_delay_uncued);
    
    %% Non-automatic sequence.
    % Check if counting answers from non-automatic sequence were correct.
    nonautodual_correct = checkCorrectCountingPerTrial(events_nonautodual);

    % Check average number of mistakes per condition (cued and uncued).
    [nonautodual_averageMistakes_cued, nonautodual_averageMistakes_uncued] = ...
        averageMistakesPerCondition(events_nonautodual, nonautodual_correct);
    
    % Get the delay in the fingertapping performance.
    [nonautodual_delay_cued, nonautodual_delay_uncued] = ...
        calculateDelayPerformance(events_nonautodual);
    
    % Add values to final average (all subjects).
    nonautodual_finalAverageMistakes_cued =...
        nonautodual_finalAverageMistakes_cued + nonautodual_averageMistakes_cued;
    nonautodual_finalAverageMistakes_uncued =...
        nonautodual_finalAverageMistakes_uncued + nonautodual_averageMistakes_uncued;
    nonautodual_finalDelay_cued =...
        nonautodual_finalDelay_cued + abs(nonautodual_delay_cued);
    nonautodual_finalDelay_uncued =...
        nonautodual_finalDelay_uncued + abs(nonautodual_delay_uncued);
    
    %% Put values of error into final struct.
    % Values for the current subject.
    s.autodual_avgMistakes_cued = autodual_averageMistakes_cued;
    s.autodual_avgMistakes_uncued = autodual_averageMistakes_uncued;
    s.nonautodual_avgMistakes_cued = nonautodual_averageMistakes_cued;
    s.nonautodual_avgMistakes_uncued = nonautodual_averageMistakes_uncued;
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

function [delay_cued, delay_uncued] = calculateDelayPerformance(events_dual)

    margin = 0.25; % Margin of error: think about what is most convenient.
    dual = events_dual.trial;
    delay_cued = 0;
    delay_uncued = 0;
    
    for trial = 1:length(dual)
        if dual(trial).cue == 1
            delay_cued = delay_cued + mean(diff(dual(trial).responses.onset) - 1/1.50);
        elseif dual(trial).cue == 0
            delay_uncued = delay_uncued + mean(diff(dual(trial).responses.onset) - 1/1.50);
        end
    end
    
    delay_cued = delay_cued / (length(dual)/2);
    delay_uncued = delay_uncued / (length(dual)/2);
    
end


% 
%     
%%
% Show if the tapping tempo was correct.
%     margin=0.25; % margin of error: think about what is most convenient
%     delay=mean(diff(events_autodual.trial(h).responses.onset)-1/1.50);
%     tempo=sprintf('The tempo was off with on average %f seconds \n', delay);
%     % Show if the tapped sequence was correct
%     correctSequence=sequences(1);
%     if all(strcmp(events_autodual.trial(h).responses.value,correctSequence{:}'))
%         seq='Sequence is correct \n \n';
%     else
%         seq='Sequence is incorrect \n \n';
%     end