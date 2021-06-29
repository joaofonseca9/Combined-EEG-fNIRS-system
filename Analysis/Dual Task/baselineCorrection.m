function [databc] = baselineCorrection (data, time)
% perform baseline correction by subtracting mean of channel over baseline 
% from channels signal 
idx_bl = find(time<0);
for iTrial = 1:size(data,3)
    for iChan = 1:size(data,1)
        databc(iChan,:,iTrial) = data(iChan,:,iTrial) - mean(data(iChan,idx_bl,iTrial));
    end
end
end