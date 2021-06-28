function [amp, badchans] = removeBadChansAmp(EEG, checkBadChan)
% Remove bad channels based on amplitude threshold: if amplitude > 50mV or < 1mV

badchans = [];
% average absolute amplitude per channel
chanavg = mean(abs(EEG.data),2);
% get channels with amplitude above threshold
amp = find((chanavg>50) | (chanavg<1));

if checkBadChan == true
    figure(2); hold on;
    barh([1:EEG.nbchan], flip(chanavg)); line([50 50], [-0.2 EEG.nbchan+1.1], 'color','k'); 
    xlim([0 55]); yticks([1:EEG.nbchan]); yticklabels(flip([{EEG.chanlocs.labels}]));
    xlabel('amplitude'); title('amp>50 | amp<1');
end

% check bad channels by looking into raw data
if checkBadChan == true    
    if isempty(amp)
        yesno = input('BADCHAN - no bad channels >50mV [yes/no]: ','s');
    else
        for i = 1:length(amp)
            disp(['BADCHAN - bad channel: ', EEG.chanlocs(amp(i)).labels]);
            yesno = input('BADCHAN - remove bad channel [yes/no]: ','s');
            if strcmp(yesno, 'yes')
                badchans = [badchans; {EEG.chanlocs(amp(i)).labels}];
            end
        end
    end 
else
    if isempty(amp)
        disp('BADCHAN - no bad channels > 50mV'); 
    else
        disp(['BADCHAN - remove bad channel: ', EEG.chanlocs(amp).labels]);
        badchans = [badchans; {EEG.chanlocs(amp).labels}];
    end
end
end