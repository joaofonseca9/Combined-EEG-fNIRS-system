function nirs_out = keepTrialsInterest(nirs_in)
% Go through the different trials and keep only the ones of interest:
% 3 - Dual Cued
% 4 - Single Cued
% 7 - Dual Uncued
% 8 - Single Uncued

nirs_out = nirs_in;
numRemovedTrials = 0;

% Go through different trials
for i=1:length(nirs_in.trialinfo)
    % If it is not a trial of interest:
    if nirs_in.trialinfo(i) ~= 3 && nirs_in.trialinfo(i) ~= 4 ...
            && nirs_in.trialinfo(i) ~= 7 ...
            && nirs_in.trialinfo(i) ~= 8
        % Eliminate that trial and all the unecessary information
        nirs_out.trial(i-numRemovedTrials) = [];
        nirs_out.time(i-numRemovedTrials) = [];
        nirs_out.trialinfo(i-numRemovedTrials) = [];
        nirs_out.sampleinfo(i-numRemovedTrials, :) = [];
        % Increase number of removed trials
        numRemovedTrials = numRemovedTrials+1;
    end 
end

end