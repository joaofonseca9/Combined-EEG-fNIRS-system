function [dual_averageMistakes_cued, dual_averageMistakes_uncued] = ...
    averageMistakesPerCondition(events_dual, dual_correct)
% Check average number of mistakes per condition (cued and uncued)

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