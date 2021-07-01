function [dual_incorrectSequences_cued, dual_incorrectSequences_uncued] = ...
    checkIncorrectSequencePerTrial(events_dual, real_sequence)
% For each trial see if the sequence was performed correctly.
% Check average number of incorrectly performed sequences per condition
% (cued and uncued).

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