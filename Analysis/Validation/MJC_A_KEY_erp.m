clear all; close all; clc;

%% Settings

% add path to correct folders and open eeglab
laptop = 'laptopJoao';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);

% select ID number and cap
subject=[{'64','28'}];


for iSub = 1:length(subject)
    sub = char(subject(iSub));
    switch sub
        case '28'
            rec='04';
        case '64'
            rec='01';
        case '02'
            rec='02';
    end
    eeglab;
    
    % set settings
    runBadTrial = false;
    checkBadChan = false;

    % create filenames to save/load data 
    mainpath_in    = 'C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\Experiment\Data\Exp\pre-processed\';
    mainpath_out     = 'C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\Experiment\Data\Exp\processed\';
    file.pstICA    = fullfile(mainpath_in,['sub-',sub],'eeg',['sub-',sub,'_rec-',rec,'_eeg_pstICA.mat']);
    file.results     = fullfile([mainpath_out,'sub-',sub,'\ResultsMatrix\sub-',sub,'_','keyVEP.mat']);
    file.results_bad = fullfile([mainpath_out,'sub-',sub,'\ResultsMatrix\sub-',sub,'_','keyVEP_bad.mat']);
    
    %% Load EEG data

    % add paths to subject folders and files
    sub_path = fullfile([mainpath_in,'sub-',sub, '\']);
    addpath(sub_path);

    % load EEG-files and create a new dataset in eeglab
    load(file.pstICA);
%     [EEG]  = pop_loadset(['sub-',sub,'_rec-',rec,'_eeg_filtData.set'],fullfile(sub_path,'eeg'));
    EEG_key=extractKeyTrials(EEG);
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG_key, 1,'setname','eegdata','gui','off');
    EEG.data = double(EEG.data);
    
    %% Bad channels with amp>50 | amp<1
    % before filtering the data, bad channels will be removed
    % if amplitude is >50 mV (noise) or <1 mV (flat line)
    if ~isfield(EEG,'badchan')
        EEG.badchan=[];
    end
    
    [~, badchans] = BadChansAmp50(EEG, checkBadChan);
    % save bad channels in struct and remove channels from data
    EEG.badchan = [EEG.badchan; badchans']; 
    if ~isempty(badchans)
        [EEG] = pop_select(EEG, 'nochannel', badchans);
    end
    clear badchans
    
    %% Extract trials 
    % create data-matrix with trials [channel x time x trial] 
    % create time vector [sec]

    [data, time] = extractEpochs (EEG);
    %% Remove bad channels within condition
    % remove bad channels if:
    % * amplitude is higher than mean+2xSD
    % * power is higher than mean+2xSD for 2/5 frequency bands
    %   [delta/theta/alpha/beta/gamma] 
    
    [badchans] = BadChansAmpPow(data, EEG, checkBadChan);
    EEG.badchan = [EEG.badchan; badchans'];
    if ~isempty(badchans)
        for i = 1:size(badchans,1)
            idx(i,1) = find(strcmp({EEG.chanlocs.labels},badchans(i)));
        end
        data(idx,:,:) = [];
        [EEG] = pop_select(EEG, 'nochannel', badchans);
    end
    clear badchans idx data
    
    %% Frequencies of Interest: 1-30Hz
    % FOI for ERP is 1-30 Hz using bandpass filter 

    [data_filt] = BandpassFilter (EEG.data, EEG.srate, [1 30]);

    % save data in EEG structure and visualize raw signal
    EEG.data = data_filt;   clear data_filt;
    if checkBadChan == true
        pop_eegplot(EEG,1,0,0);
    end
    
    %% Extract filtered (FOI) trials 
    % create data-matrix with trials [channel x time x trial] 
    % (data is filtered for frequencies of interest)

    [data, ~] = extractEpochs (EEG);
    
    %% Baseline correction
    % extract average of channel's baseline from channel's signal per trial

    [data] = BaselineCorrection (data, time);

    %% Remove bad channels from trials
    % remove bad channels from trial if: 
    % * amplitude is higher than mean+2xSD of channel
    % * power is higher than mean+2xSD of channel in 2/4 frequency bands
    %  [delta/theta/alpha/beta] 
    % replace bad channels within trial with NaN

    if runBadTrial == true
        [badtrial] = BadChansTrial (data, EEG);
        datanan = data;
        for iTrial = 1:size(datanan,3)
            for iChan = 1:size(datanan,1)
                if badtrial(iChan,iTrial) == 1
                    datanan(iChan,:,iTrial) = NaN;
                end 
            end
        end
        data = datanan; clear datanan;
    end

    %% Average VEP over trials 
    % average over trials, use omit NaNs (bad channels within trials)

    ERP = mean(data,3,'omitnan');
    
    %% Global field power in time domain (GFPt)
    % calculate difference between two channels, per time point 
    % take the sum of the differences, per time point
    % gives the GFPt over time

    [GFPt]  = GlobalFieldPotential(ERP);
    KEYerp.data      = data;
    KEYerp.ERP       = ERP;
    KEYerp.GFPt      = GFPt; 
    KEYerp.time      = time;
    KEYerp.chans     = {EEG.chanlocs.labels};
    
%     LETTERerpBAD.chans     = {EEG.chanlocs.labels};
%     LETTERerpBAD.badchans  = EEG.badchan;
%     LETTERerpBAD.badtrial  = badtrial;

    save(file.results,'KEYerp');
%     save(file.results_bad,'CHECKerpBAD');

    %% Plotting
%     time(end+1)=time(end);
    h=figure;
    subplot(1,2,1); plot(time,mean(ERP,3,'omitnan'),'b');
    subplot(1,2,2); hold on; h1 = fill([time,fliplr(time)], [GFPt,fliplr(GFPt)],'b','LineStyle','none');
    set(h1,'FaceAlpha',0.4); plot(time,GFPt,'r','LineWidth',1.5);
    
    subplot(1,2,1); set(gca,'FontSize',11); box on;
    ylabel('Potential (\muV)','FontSize',14); title('ERP','FontSize',14); ylim([-6 6])
    xticks([-0.1:0.1:0.4]); yticks([-5:5:5]); line([0 0],[-6 6],'Color','k');xlabel('Time (s)','FontSize',14);
    subplot(1,2,2); set(gca,'FontSize',11); box on;
    ylabel('GFPt (\muV)','FontSize',14); xlabel('Time (s)','FontSize',14); ylim([-2 26])
    xticks([-0.1:0.1:0.4]); yticks([0:10:20]); line([0 0],[-2 25],'Color','k');  title('GFPt');
    
    
    saveas(gcf,['sub-',sub,'_KEY_erp.jpg'])
end


%% HELPER FUNCTIONS

function [EEG_key]=extractKeyTrials(EEG)
event_samp  = [EEG.event.latency];

%LETTER TRIALS
Key=event_samp((strcmp({EEG.event.type}, 's1777'))==1);
starts=Key-floor(0.5*EEG.srate);
stops=Key+ceil(0.5*EEG.srate);


rej=[starts(1)-floor(0.1*EEG.srate) starts(1)];
for ii=1:length(stops)-1
    rej=[rej;stops(ii) starts(ii+1)];
    if ii==length(stops)-1
        rej=[rej;stops(ii+1) event_samp(end)];
    end
end

[EEG_key] = eeg_eegrej(EEG, rej);
end

%%
function [amp, badchans] = BadChansAmp50(EEG, checkBadChan)
% remove bad channels based on amplitude threshold:
% if amplitude >50mV or <1mV

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
            disp('BADCHAN - no bad channels >50mV'); 
        else
            disp(['BADCHAN - remove bad channel: ', EEG.chanlocs(amp).labels]);
            badchans = [badchans; {EEG.chanlocs(amp).labels}];
        end
    end
end


%%
function [data, time] = extractEpochs (EEG)
% extract epochs (trials) using the event-markers [channels x time x trial]
% ERP-trial: -0.1s to 0.4s
% create time vector in seconds

    % task-events: start (s1555) / stop (s1500) / flip (s1255)
    event_samp  = [EEG.event.latency];
    event_key  = event_samp((strcmp({EEG.event.type}, 's1777'))==1);
    
    % extract trials ERP
    for iTrial = 1:length(event_key)
        data(:,:,iTrial) = EEG.data(:, [event_key(iTrial)-floor(0.1*EEG.srate) : event_key(iTrial)+ceil(0.4*EEG.srate)]);
    end
    % create time vector [s]
    time = -0.1 : 1/EEG.srate : 0.4;
end

%%
function [badchans] = BadChansAmpPow (data, EEG, checkBadChan)

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

%%
function [data_filt] = BandpassFilter (data, fs, fc)
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

%% 
function [databc] = BaselineCorrection (data, time)
% perform baseline correction by subtracting mean of channel over baseline 
% from channels signal 

    idx_bl = find(time<0);
    for iTrial = 1:size(data,3)
        for iChan = 1:size(data,1)
            databc(iChan,:,iTrial) = data(iChan,:,iTrial) - mean(data(iChan,idx_bl,iTrial));
        end
    end
end

%%
function [badtrial] = BadChansTrial (data, EEG)

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


%%
function [GFPt] = GlobalFieldPotential(ERP)
% calculate difference between two channels, per time point 
% take the sum of the differences, per time point
% gives the GFPt over time

    for iTime = 1:size(ERP,2)
        Udiff = zeros(size(ERP,1), size(ERP,1));
        for i = 1:size(ERP,1)
            for j = 1:size(ERP,1)
                Udiff(i,j) = (ERP(i,iTime)-ERP(j,iTime)).^2;
            end
        end
        GFPt(1,iTime) = sqrt( (1/(2*size(ERP,1))) * sum(Udiff,'all') );
    end
end



