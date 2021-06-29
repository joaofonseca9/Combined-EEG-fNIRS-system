function [badtrial] = badChansTrial (data, EEG)

% AMPLITUDE...
% calculate channel average per trial, calculate channel threshold
% over all trials (mean+2xSD) and find channels above threshold

chanavg = squeeze(mean(abs(data),2));
chanth  = mean(abs(data),[2 3]) + 2*std(abs(data),[],[2 3]);
amp = (chanavg > chanth);

% POWER...
% calculate the power and power threshold (mean+2xSD) per band, 
% and find channels with power above threshold in 2/4 frequency bands 

freqs = [{'delta'};{'theta'};{'alpha'};{'beta'}];
fcut  = [1 4; 4 8; 8 14; 14 30];
% calculate power
for iTrial = 1:size(data,3)
    Tlength = size(data,2);
    [P(:,:,iTrial),f] = periodogram(data(:,:,iTrial)', hann(Tlength), 2^(2+nextpow2(Tlength)),EEG.srate);
end
for iF = 1:length(freqs)
    % get power per freq band and threshold, find channels above threshold
    chanpow.(freqs{iF})     = P((f(:,1)>=fcut(iF,1) & f(:,1)<=fcut(iF,2)),:,:);
    chanpowth.(freqs{iF})   = mean(chanpow.(freqs{iF}),[1 3]) + 2*std(chanpow.(freqs{iF}),[],[1 3]);
    pow(iF,:,:)             = mean(chanpow.(freqs{iF}),1) > chanpowth.(freqs{iF}); 
end
pow = squeeze(sum(pow,1));

% check for channels within trials with amp>th or pow>th (in 2/4 freq bands)
for iTrial = 1:size(data,3)
    for iChan = 1:size(data,1)
        if amp(iChan,iTrial) == 1 || pow(iChan,iTrial)>=2
            badtrial(iChan,iTrial) = 1;
        else
            badtrial(iChan,iTrial) = 0;
        end
    end
end

end