function [filt_data] = filterNoise(data, EEG, order)
% Filter the data with a 2nd-order Butterworth bandpass filter [0-48Hz],
% to extract frequencies of interest.
% Filter the data with a 2nd-order Butterworth bandstop filter [48-52Hz],
% to eliminate 50-Hz line noise.

% Settings 2nd-order Butterworth bandpass/stop filter
bandpass = 2*(1/EEG.srate)*[1 48];
bandstop = 2*(1/EEG.srate)*[48 52];  

[B,A] = butter(order, bandpass, 'bandpass');
[b,a] = butter(order, bandstop, 'stop');

% Preallocating...
f_data = zeros(size(EEG.data));
filt_data = zeros(size(EEG.data));

% Filter data
for iChan = 1:EEG.nbchan
    f_data(iChan,:) = filtfilt(B, A, data(iChan, :));
   filt_data(iChan,:)  = filtfilt(b, a, f_data(iChan, :));
end

end