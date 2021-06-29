function [data_filt] = bandpassFilter(data, fs, fc)
% filter data using specified cut-off frequency
band = 2* (1/fs) * fc;
[B,A] = butter(2, band, 'bandpass');
data_filt = zeros(size(data));
for iTrial = 1:size(data,3)
    for iChan = 1:size(data,1)
        data_filt(iChan,:,iTrial) = filtfilt(B,A, data(iChan,:,iTrial));
    end
end
end