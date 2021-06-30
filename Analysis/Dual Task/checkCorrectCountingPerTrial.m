function dual_correct = checkCorrectCountingPerTrial(events_dual)
% For each trial see if the counting response was equal to the actual
% number of G's

dual = events_dual.trial;
dual_correct = zeros(1, length(dual));

for trial = 1:length(dual)
    stimuli = dual(trial).stimuli;
    if count(cell2mat(stimuli.value),'G') == str2num(cell2mat(stimuli.response))
        dual_correct(trial) = 1;
    end
end

end