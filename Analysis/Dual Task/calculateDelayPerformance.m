function [delay_cued, delay_uncued] = calculateDelayPerformance(events_dual, real_sequence)
% For each trial see if the sequence was performed correctly and if so
% check the delay in the performance.
% Check average delay on correctly performed sequences per condition
% (cued and uncued).

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