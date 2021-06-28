function [badchans] = badChansAmpPow (data, EEG, checkBadChan)
% AMPLITUDE...
% calculate channel average over trials, calculate channel threshold
% (mean+2xSD) and find channels above threshold

chanavg = mean(abs(data),[2 3]);
chanth  = mean(abs(data),'all') + 2*std(abs(data),[],'all');
amp = (chanavg > chanth)';
% figure
if checkBadChan == true
    figure; subplot(1,6,1); hold on;
    barh([1:EEG.nbchan], flip(chanavg)); line([chanth chanth], [-0.2 EEG.nbchan+1.1], 'color','k');
    xlim([0 chanth+0.1*chanth]); yticks([1:EEG.nbchan]); yticklabels(flip([{EEG.chanlocs.labels}]));
    xlabel('amplitude'); title('amp>mean+2SD');
end

% POWER...
% calculate the power and power threshold (mean+2xSD) per band, 
% and find channels with power above threshold in 2/5 frequency bands 

freqs = [{'delta'};{'theta'};{'alpha'};{'beta'};{'gamma'}];
fcut  = [1 4; 4 8; 8 14; 14 30; 30 48];

% calculate power
for iTrial = 1:size(data,3)
    Tlength = size(data,2);
    [P(:,:,iTrial),f] = periodogram(data(:,:,iTrial)', hann(Tlength), 2^(2+nextpow2(Tlength)),EEG.srate);
end
for iF = 1:length(freqs)        
    % get power per freq band and threshold, find channels above threshold
    chanpow.(freqs{iF})     = mean( P((f(:,1)>=fcut(iF,1) & f(:,1)<=fcut(iF,2)),:,:), 3);
    chanpowth.(freqs{iF})   = mean(chanpow.(freqs{iF}),'all') + 2*std(chanpow.(freqs{iF}),[],'all');
    pow(iF,:)               = mean(chanpow.(freqs{iF}),1) > chanpowth.(freqs{iF});

    if checkBadChan == true
        subplot(1,6,iF+1);  hold on;
        barh([1:EEG.nbchan], flip(mean(chanpow.(freqs{iF}),1)));
        line([chanpowth.(freqs{iF}) chanpowth.(freqs{iF})], [-0.2 EEG.nbchan+1.1], 'color','k');
        xlim([0 chanpowth.(freqs{iF})+0.1*chanpowth.(freqs{iF})]);
        yticks([1:EEG.nbchan]); yticklabels(flip([{EEG.chanlocs.labels}])); xlabel(freqs{iF}); title('pow>mean+2SD');
    end
end

bad = [];
for iChan = 1:EEG.nbchan
    if sum(amp(:,iChan)) >=1 || sum(pow(:,iChan)) >= 2
        bad = [bad; iChan];
    end
end

badchans = [];
if checkBadChan == true
    % check bad channels by looking into raw data and power spectrum    
    if isempty(bad)
        yesno = input('BADCHAN - no bad channels >mean+2SD [yes/no]: ','s');
    else
        for i = 1:length(bad)
            disp(['BADCHAN - bad channel: ', EEG.chanlocs(bad(i)).labels]);
            yesno = input('BADCHAN - remove bad channel [yes/no]: ','s');
            if strcmp(yesno, 'yes')
                badchans = [badchans; {EEG.chanlocs(bad(i)).labels}];
            end
        end
    end
else
    if isempty(bad)
        disp('BADCHAN - no bad channels >mean+2SD'); 
    else
        disp(['BADCHAN - remove bad channel: ', EEG.chanlocs(bad).labels]);
        badchans = [badchans; {EEG.chanlocs(bad).labels}];
    end
end
end
