%% Analysis of the EEG signals - PSD.

clear; clc; close all;
addpath('C:\Users\maria\OneDrive\Documentos\GitHub\Combined-EEG-fNIRS-system\Analysis');
addpath('C:\Users\maria\OneDrive\Documentos\GitHub\Combined-EEG-fNIRS-system\Analysis\Automatic_vs_NonAutomatic');

laptop = 'laptopMariana';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);
eeglab;
ft_defaults;
results_path = 'C:\Users\maria\OneDrive\Ambiente de Trabalho\Automaticity Results\Separate Channels';

subrec = ["28" "04"; "02" "02"; "76" "01"];

% List of the 30 channels present in the cap.
list_channels = ["Fp1"; "Fpz"; "Fp2"; "F7"; "F3"; "AFFz"; "F4"; "F8";...
    "FC5"; "FC1"; "FC2"; "FC6"; "T7"; "C3"; "Cz"; "C4"; "T8"; "CP5";...
    "CP1"; "CP2"; "CP6"; "P7"; "P3"; "Pz"; "P4"; "P8"; "POz"; "O1";...
    "Oz"; "O2"];

% Loop through every subject.
for subject = 1:size(subrec, 1)
    sub = subrec(subject, 1);
    rec = subrec(subject, 2);
    
    % Load the subject's EEG signals.
    load([mainpath_in, '\pre-processed\sub-', char(sub), '\eeg\sub-',...
        char(sub), '_rec-', char(rec), '_eeg_divided.mat']);
    
    % Separate into the four different tasks.
    EEG_AutoUncued = EEG_divided.EEG_AutoNoCue;
    EEG_NonAutoUncued = EEG_divided.EEG_NonAutoNoCue;
    EEG_AutoCued = EEG_divided.EEG_AutoCue;
    EEG_NonAutoCued = EEG_divided.EEG_NonAutoCue;
    
    % Get the locations of the channels of interest for the first subject
    % since he has all the channels.
    if subject==1
        locs = {EEG_AutoCued.chanlocs.labels};
        % DLPFC.
        F7_loc = find(contains(locs, 'F7'));
        F8_loc = find(contains(locs, 'F8'));
        % SMA.
        FC1_loc = find(contains(locs, 'FC1'));
        FC2_loc = find(contains(locs, 'FC2'));
        Cz_loc = find(contains(locs, 'Cz'));
        % M1.
        C3_loc = find(contains(locs, 'C3'));
        % PPC.
        P3_loc = find(contains(locs, 'P3'));
        P4_loc = find(contains(locs, 'P4'));
    end
    
    % Calculate threshold to eliminate noisy epochs.
    th = calculateThreshold(EEG_divided);
    
    %% Auto Uncued.
    
    event_samp  = [EEG_AutoUncued.event.latency];
    startTask = find(strcmp({EEG_AutoUncued.event.type}, 's1703')==1);
    endTask = find(strcmp({EEG_AutoUncued.event.type}, 's1711')==1);
    keypresses = find(strcmp({EEG_AutoUncued.event.type}, 's1777')==1);
    
    % Get the power spectrum density (PSD) averaged over all trials.
    [power, freq] = calculateAveragePowerAllTrials(EEG_AutoUncued,...
        event_samp, startTask, endTask, keypresses, th);
   
    % Compensate for removed channels.
    power = compensateRemovedChannels(power, EEG_AutoUncued, list_channels);
    
    % Save the values onto a allSubjects variable.
    autouncued_power_allSubjects(:, :, subject) = power;
    
    %% Non-Auto Uncued.
    
    event_samp  = [EEG_NonAutoUncued.event.latency];
    startTask = find(strcmp({EEG_NonAutoUncued.event.type}, 's1705')==1);
    endTask = find(strcmp({EEG_NonAutoUncued.event.type}, 's1713')==1);
    keypresses = find(strcmp({EEG_NonAutoUncued.event.type}, 's1777')==1);
    
    % Get the power spectrum density (PSD) averaged over all trials.
    [power, freq] = calculateAveragePowerAllTrials(EEG_NonAutoUncued,...
        event_samp, startTask, endTask, keypresses, th);
    
    % Compensate for removed channels.
    power = compensateRemovedChannels(power, EEG_NonAutoUncued, list_channels);
    
    % Save the values onto a allSubjects variable.
    nonautouncued_power_allSubjects(:, :, subject) = power;
    
    %% Auto Cued.
    
    event_samp  = [EEG_AutoCued.event.latency];
    startTask = find(strcmp({EEG_AutoCued.event.type}, 's1702')==1);
    endTask = find(strcmp({EEG_AutoCued.event.type}, 's1710')==1);
    keypresses = find(strcmp({EEG_AutoUncued.event.type}, 's1777')==1);
    
    % Get the power spectrum density (PSD) averaged over all trials.
    [power, freq] = calculateAveragePowerAllTrials(EEG_AutoCued,...
        event_samp, startTask, endTask, keypresses, th);
    
    % Compensate for removed channels.
    power = compensateRemovedChannels(power, EEG_AutoCued, list_channels);
    
    % Save the values onto a allSubjects variable.
    autocued_power_allSubjects(:, :, subject) = power;
    
    %% Non-Auto Cued.
    
    event_samp  = [EEG_NonAutoCued.event.latency];
    startTask = find(strcmp({EEG_NonAutoCued.event.type}, 's1704')==1);
    endTask = find(strcmp({EEG_NonAutoCued.event.type}, 's1712')==1);
    keypresses = find(strcmp({EEG_AutoUncued.event.type}, 's1777')==1);
    
    % Get the power spectrum density (PSD) averaged over all trials.
    [power, freq] = calculateAveragePowerAllTrials(EEG_NonAutoCued,...
        event_samp, startTask, endTask, keypresses, th);
    
    % Compensate for removed channels.
    power = compensateRemovedChannels(power, EEG_NonAutoCued, list_channels);
    
    % Save the values onto a allSubjects variable.
    nonautocued_power_allSubjects(:, :, subject) = power;
    
    %% Get the values of power for the regions of interest and plot them.
    
    % DLPFC.
    
    autouncued_power_F7 = autouncued_power_allSubjects(:, F7_loc, subject);
    autocued_power_F7 = autocued_power_allSubjects(:, F7_loc, subject);
    nonautouncued_power_F7 = nonautouncued_power_allSubjects(:, F7_loc, subject);
    nonautocued_power_F7 = nonautocued_power_allSubjects(:, F7_loc, subject);
    
    autouncued_power_F8 = autouncued_power_allSubjects(:, F8_loc, subject);
    autocued_power_F8 = autocued_power_allSubjects(:, F8_loc, subject);
    nonautouncued_power_F8 = nonautouncued_power_allSubjects(:, F8_loc, subject);
    nonautocued_power_F8 = nonautocued_power_allSubjects(:, F8_loc, subject);
    
    autouncued_power_DLPFC = mean([autouncued_power_F7 autouncued_power_F8], 2);
    autocued_power_DLPFC = mean([autocued_power_F7 autocued_power_F8], 2);
    nonautouncued_power_DLPFC = mean([nonautouncued_power_F7 nonautouncued_power_F8], 2);
    nonautocued_power_DLPFC = mean([nonautocued_power_F7 nonautocued_power_F8], 2);
    
    figure; title('DLPFC');
    plot(freq, autouncued_power_DLPFC, '-b'); hold on;
    plot(freq, autocued_power_DLPFC, '-m'); hold on;
    plot(freq, nonautouncued_power_DLPFC, '-g'); hold on;
    plot(freq, nonautocued_power_DLPFC, '-r'); hold on;
    xline(4); hold on;
    xline(8); hold on;
    xline(13); hold on;
    xline(32); hold off;
    legend('Auto Uncued','Auto Cued', 'Non-Auto Uncued', 'Non-Auto Cued');
    
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, fullfile(results_path, ['Sub-', char(sub)], 'DLPFC'),'png');
    
    % Save the values onto a subject struct.
    s.autouncued_power_DLPFC = autouncued_power_DLPFC;
    s.autocued_power_DLPFC = autocued_power_DLPFC;
    s.nonautouncued_power_DLPFC = nonautouncued_power_DLPFC;
    s.nonautocued_power_FDLPFC= nonautocued_power_DLPFC;
    
    % SMA.
    
    autouncued_power_FC1 = autouncued_power_allSubjects(:, FC1_loc, subject);
    autocued_power_FC1 = autocued_power_allSubjects(:, FC1_loc, subject);
    nonautouncued_power_FC1 = nonautouncued_power_allSubjects(:, FC1_loc, subject);
    nonautocued_power_FC1 = nonautocued_power_allSubjects(:, FC1_loc, subject);
    
    autouncued_power_FC2 = autouncued_power_allSubjects(:, FC2_loc, subject);
    autocued_power_FC2 = autocued_power_allSubjects(:, FC2_loc, subject);
    nonautouncued_power_FC2 = nonautouncued_power_allSubjects(:, FC2_loc, subject);
    nonautocued_power_FC2 = nonautocued_power_allSubjects(:, FC2_loc, subject);
    
    autouncued_power_Cz = autouncued_power_allSubjects(:, Cz_loc, subject);
    autocued_power_Cz = autocued_power_allSubjects(:, Cz_loc, subject);
    nonautouncued_power_Cz = nonautouncued_power_allSubjects(:, Cz_loc, subject);
    nonautocued_power_Cz = nonautocued_power_allSubjects(:, Cz_loc, subject);
    
    autouncued_power_SMA = mean([autouncued_power_FC1 autouncued_power_FC2...
        autouncued_power_Cz], 2);
    autocued_power_SMA = mean([autocued_power_FC1 autocued_power_FC2...
        autocued_power_Cz], 2);
    nonautouncued_power_SMA = mean([nonautouncued_power_FC1 nonautouncued_power_FC2...
        nonautouncued_power_Cz], 2);
    nonautocued_power_SMA = mean([nonautocued_power_FC1 nonautocued_power_FC2...
        nonautocued_power_Cz], 2);
    
    figure; title('SMA');
    plot(freq, autouncued_power_SMA, '-b'); hold on;
    plot(freq, autocued_power_SMA, '-m'); hold on;
    plot(freq, nonautouncued_power_SMA, '-g'); hold on;
    plot(freq, nonautocued_power_SMA, '-r'); hold on;
    xline(4); hold on;
    xline(8); hold on;
    xline(13); hold on;
    xline(32); hold off;
    legend('Auto Uncued','Auto Cued', 'Non-Auto Uncued', 'Non-Auto Cued');
    
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, fullfile(results_path, ['Sub-', char(sub)], 'SMA'),'png');
    
    % Save the values onto a subject struct.
    s.autouncued_power_SMA = autouncued_power_SMA;
    s.autocued_power_SMA = autocued_power_SMA;
    s.nonautouncued_power_SMA = nonautouncued_power_SMA;
    s.nonautocued_power_SMA = nonautocued_power_SMA;
    
    % M1.
    
    autouncued_power_M1 = autouncued_power_allSubjects(:, C3_loc, subject);
    autocued_power_M1 = autocued_power_allSubjects(:, C3_loc, subject);
    nonautouncued_power_M1 = nonautouncued_power_allSubjects(:, C3_loc, subject);
    nonautocued_power_M1 = nonautocued_power_allSubjects(:, C3_loc, subject);
    
    figure; title('M1');
    plot(freq, autouncued_power_M1, '-b'); hold on;
    plot(freq, autocued_power_M1, '-m'); hold on;
    plot(freq, nonautouncued_power_M1, '-g'); hold on;
    plot(freq, nonautocued_power_M1, '-r'); hold on;
    xline(4); hold on;
    xline(8); hold on;
    xline(13); hold on;
    xline(32); hold off;
    legend('Auto Uncued','Auto Cued', 'Non-Auto Uncued', 'Non-Auto Cued');
    
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, fullfile(results_path, ['Sub-', char(sub)], 'M1'),'png');
    
    % Save the values onto a subject struct.
    s.autouncued_power_M1 = autouncued_power_M1;
    s.autocued_power_M1 = autocued_power_M1;
    s.nonautouncued_power_M1 = nonautouncued_power_M1;
    s.nonautocued_power_M1 = nonautocued_power_M1;
    
    % PPC.
    
    autouncued_power_P3 = autouncued_power_allSubjects(:, P3_loc, subject);
    autocued_power_P3 = autocued_power_allSubjects(:, P3_loc, subject);
    nonautouncued_power_P3 = nonautouncued_power_allSubjects(:, P3_loc, subject);
    nonautocued_power_P3 = nonautocued_power_allSubjects(:, P3_loc, subject);
    
    autouncued_power_P4 = autouncued_power_allSubjects(:, P4_loc, subject);
    autocued_power_P4 = autocued_power_allSubjects(:, P4_loc, subject);
    nonautouncued_power_P4 = nonautouncued_power_allSubjects(:, P4_loc, subject);
    nonautocued_power_P4 = nonautocued_power_allSubjects(:, P4_loc, subject);
    
    autouncued_power_PPC = mean([autouncued_power_P3 autouncued_power_P4], 2);
    autocued_power_PPC = mean([autocued_power_P3 autocued_power_P4], 2);
    nonautouncued_power_PPC = mean([nonautouncued_power_P3 nonautouncued_power_P4], 2);
    nonautocued_power_PPC = mean([nonautocued_power_P3 nonautocued_power_P4], 2);
    
    figure; title('PPC');
    plot(freq, autouncued_power_PPC, '-b'); hold on;
    plot(freq, autocued_power_PPC, '-m'); hold on;
    plot(freq, nonautouncued_power_PPC, '-g'); hold on;
    plot(freq, nonautocued_power_PPC, '-r'); hold on;
    xline(4); hold on;
    xline(8); hold on;
    xline(13); hold on;
    xline(32); hold off;
    legend('Auto Uncued','Auto Cued', 'Non-Auto Uncued', 'Non-Auto Cued');
    
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, fullfile(results_path, ['Sub-', char(sub)], 'PPC'),'png');
    
    % Save the values onto a subject struct.
    s.autouncued_power_PPC = autouncued_power_PPC;
    s.autocued_power_PPC = autocued_power_PPC;
    s.nonautouncued_power_PPC = nonautouncued_power_PPC;
    s.nonautocued_power_PPC = nonautocued_power_PPC;
    
    % Add struct of current subject to all subjects struct.
    allsubs.(genvarname(strcat('sub', char(sub)))) = s;
    
    disp(['These are the results for subject ', char(sub), '.']);
    disp('Press any key to move onto the next subject.');
    pause;
    close all;
    
end

% Get the power spectrum density (PSD) averaged over all subjects.
% Auto Uncued.
autouncued_power = mean(autouncued_power_allSubjects, 3, 'omitnan');
% Non-Auto Uncued.
nonautouncued_power = mean(nonautouncued_power_allSubjects, 3, 'omitnan');
% Auto Cued.
autocued_power = mean(autocued_power_allSubjects, 3, 'omitnan');
% Non-Auto Cued.
nonautocued_power = mean(nonautocued_power_allSubjects, 3, 'omitnan');

% Get the standard deviation over all subjects.
% Auto Uncued.
std_autouncued_power = std(autouncued_power_allSubjects, 1, 3, 'omitnan');
% Non-Auto Uncued.
std_nonautouncued_power = std(nonautouncued_power_allSubjects, 1, 3, 'omitnan');
% Auto Cued.
std_autocued_power = std(autocued_power_allSubjects, 1, 3, 'omitnan');
% Non-Auto Cued.
std_nonautocued_power = std(nonautocued_power_allSubjects, 1, 3, 'omitnan');

%% Plot the PSD for specific regions of interest.

%% DLPFC.
autouncued_power_F7 = autouncued_power(:, F7_loc);
autocued_power_F7 = autocued_power(:, F7_loc);
nonautouncued_power_F7 = nonautouncued_power(:, F7_loc);
nonautocued_power_F7 = nonautocued_power(:, F7_loc);
std_autouncued_power_F7 = std_autouncued_power(:, F7_loc);
std_autocued_power_F7 = std_autocued_power(:, F7_loc);
std_nonautouncued_power_F7 = std_nonautouncued_power(:, F7_loc);
std_nonautocued_power_F7 = std_nonautocued_power(:, F7_loc);

autouncued_power_F8 = autouncued_power(:, F8_loc);
autocued_power_F8 = autocued_power(:, F8_loc);
nonautouncued_power_F8 = nonautouncued_power(:, F8_loc);
nonautocued_power_F8 = nonautocued_power(:, F8_loc);
std_autouncued_power_F8 = std_autouncued_power(:, F8_loc);
std_autocued_power_F8 = std_autocued_power(:, F8_loc);
std_nonautouncued_power_F8 = std_nonautouncued_power(:, F8_loc);
std_nonautocued_power_F8 = std_nonautocued_power(:, F8_loc);

autouncued_power_DLPFC = mean([autouncued_power_F7 autouncued_power_F8], 2);
autocued_power_DLPFC = mean([autocued_power_F7 autocued_power_F8], 2);
nonautouncued_power_DLPFC = mean([nonautouncued_power_F7 nonautouncued_power_F8], 2);
nonautocued_power_DLPFC = mean([nonautocued_power_F7 nonautocued_power_F8], 2);
std_autouncued_power_DLPFC = mean([std_autouncued_power_F7 std_autouncued_power_F8], 2);
std_autocued_power_DLPFC = mean([std_autocued_power_F7 std_autocued_power_F8], 2);
std_nonautouncued_power_DLPFC = mean([std_nonautouncued_power_F7 std_nonautouncued_power_F8], 2);
std_nonautocued_power_DLPFC = mean([std_nonautocued_power_F7 std_nonautocued_power_F8], 2);

figure; title('DLPFC');
tiledlayout(1, 2);
% Automatic sequence.
ax1 = nexttile;
plot(freq, autouncued_power_DLPFC, '-b', 'LineWidth', 1.5); 
hold on;
h1 = patch([freq; flipud(freq)],...
    [(autouncued_power_DLPFC - std_autouncued_power_DLPFC);...
    flipud(autouncued_power_DLPFC + std_autouncued_power_DLPFC)], 'b',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, autocued_power_DLPFC, '-m', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(autocued_power_DLPFC - std_autocued_power_DLPFC);...
    flipud(autocued_power_DLPFC + std_autocued_power_DLPFC)], 'm',...
    'LineStyle', 'none'); set(h2,'FaceAlpha',0.2);
hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
xlim([4 48]);
legend('Auto Uncued', '', 'Auto Cued', '');
xlabel('Frequency (Hz)');
ylabel('PSD');
% Non-automatic sequence.
ax2 = nexttile;
plot(freq, nonautouncued_power_DLPFC, '-g', 'LineWidth', 1.5);
hold on;
h1 = patch([freq; flipud(freq)],...
    [(nonautouncued_power_DLPFC - std_nonautouncued_power_DLPFC);...
    flipud(nonautouncued_power_DLPFC + std_nonautouncued_power_DLPFC)], 'g',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, nonautocued_power_DLPFC, '-r', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(nonautocued_power_DLPFC - std_nonautocued_power_DLPFC);...
    flipud(nonautocued_power_DLPFC + std_nonautocued_power_DLPFC)], 'r',...
    'LineStyle', 'none'); set(h2,'FaceAlpha',0.2);
hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
xlim([4 48]);
legend('Non-Auto Uncued', '', 'Non-Auto Cued', '');
xlabel('Frequency (Hz)');
ylabel('PSD');
linkaxes([ax1 ax2],'y')

set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, fullfile(results_path, 'DLPFC_PowervsFreq'),'png');

%% SMA.
autouncued_power_FC1 = autouncued_power(:, FC1_loc);
autocued_power_FC1 = autocued_power(:, FC1_loc);
nonautouncued_power_FC1 = nonautouncued_power(:, FC1_loc);
nonautocued_power_FC1 = nonautocued_power(:, FC1_loc);
std_autouncued_power_FC1 = std_autouncued_power(:, FC1_loc);
std_autocued_power_FC1 = std_autocued_power(:, FC1_loc);
std_nonautouncued_power_FC1 = std_nonautouncued_power(:, FC1_loc);
std_nonautocued_power_FC1 = std_nonautocued_power(:, FC1_loc);

autouncued_power_FC2 = autouncued_power(:, FC2_loc);
autocued_power_FC2 = autocued_power(:, FC2_loc);
nonautouncued_power_FC2 = nonautouncued_power(:, FC2_loc);
nonautocued_power_FC2 = nonautocued_power(:, FC2_loc);
std_autouncued_power_FC2 = std_autouncued_power(:, FC2_loc);
std_autocued_power_FC2 = std_autocued_power(:, FC2_loc);
std_nonautouncued_power_FC2 = std_nonautouncued_power(:, FC2_loc);
std_nonautocued_power_FC2 = std_nonautocued_power(:, FC2_loc);

autouncued_power_Cz = autouncued_power(:, Cz_loc);
autocued_power_Cz = autocued_power(:, Cz_loc);
nonautouncued_power_Cz = nonautouncued_power(:, Cz_loc);
nonautocued_power_Cz = nonautocued_power(:, Cz_loc);
std_autouncued_power_Cz = std_autouncued_power(:, Cz_loc);
std_autocued_power_Cz = std_autocued_power(:, Cz_loc);
std_nonautouncued_power_Cz = std_nonautouncued_power(:, Cz_loc);
std_nonautocued_power_Cz = std_nonautocued_power(:, Cz_loc);

autouncued_power_SMA = mean([autouncued_power_FC1 autouncued_power_FC2...
    autouncued_power_Cz], 2);
autocued_power_SMA = mean([autocued_power_FC1 autocued_power_FC2...
    autocued_power_Cz], 2);
nonautouncued_power_SMA = mean([nonautouncued_power_FC1 nonautouncued_power_FC2...
    nonautouncued_power_Cz], 2);
nonautocued_power_SMA = mean([nonautocued_power_FC1 nonautocued_power_FC2...
    nonautocued_power_Cz], 2);
std_autouncued_power_SMA = mean([std_autouncued_power_FC1 std_autouncued_power_FC2...
    std_autouncued_power_Cz], 2);
std_autocued_power_SMA = mean([std_autocued_power_FC1 std_autocued_power_FC2...
    std_autocued_power_Cz], 2);
std_nonautouncued_power_SMA = mean([std_nonautouncued_power_FC1 std_nonautouncued_power_FC2...
    std_nonautouncued_power_Cz], 2);
std_nonautocued_power_SMA = mean([std_nonautocued_power_FC1 std_nonautocued_power_FC2...
    std_nonautocued_power_Cz], 2);

figure; title('SMA');
tiledlayout(1, 2);
% Automatic sequence.
ax1 = nexttile;
plot(freq, autouncued_power_SMA, '-b', 'LineWidth', 1.5); 
hold on;
h1 = patch([freq; flipud(freq)],...
    [(autouncued_power_SMA - std_autouncued_power_SMA);...
    flipud(autouncued_power_SMA + std_autouncued_power_SMA)], 'b',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, autocued_power_SMA, '-m', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(autocued_power_SMA - std_autocued_power_SMA);...
    flipud(autocued_power_SMA + std_autocued_power_SMA)], 'm',...
    'LineStyle', 'none'); set(h2,'FaceAlpha',0.2);
hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
xlim([4 48]);
legend('Auto Uncued', '', 'Auto Cued', '');
xlabel('Frequency (Hz)');
ylabel('PSD');
% Non-automatic sequence.
ax2 = nexttile;
plot(freq, nonautouncued_power_SMA, '-g', 'LineWidth', 1.5);
hold on;
h1 = patch([freq; flipud(freq)],...
    [(nonautouncued_power_SMA - std_nonautouncued_power_SMA);...
    flipud(nonautouncued_power_SMA + std_nonautouncued_power_SMA)], 'g',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, nonautocued_power_SMA, '-r', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(nonautocued_power_SMA - std_nonautocued_power_SMA);...
    flipud(nonautocued_power_SMA + std_nonautocued_power_SMA)], 'r',...
    'LineStyle', 'none'); set(h2,'FaceAlpha',0.2);
hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
xlim([4 48]);
legend('Non-Auto Uncued', '', 'Non-Auto Cued', '');
xlabel('Frequency (Hz)');
ylabel('PSD');
linkaxes([ax1 ax2],'y')

set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, fullfile(results_path, 'SMA_PowervsFreq'),'png');

%% M1.
autouncued_power_M1 = autouncued_power(:, C3_loc);
autocued_power_M1 = autocued_power(:, C3_loc);
nonautouncued_power_M1 = nonautouncued_power(:, C3_loc);
nonautocued_power_M1 = nonautocued_power(:, C3_loc);
std_autouncued_power_M1 = std_autouncued_power(:, C3_loc);
std_autocued_power_M1 = std_autocued_power(:, C3_loc);
std_nonautouncued_power_M1 = std_nonautouncued_power(:, C3_loc);
std_nonautocued_power_M1 = std_nonautocued_power(:, C3_loc);

figure; title('M1');
tiledlayout(1, 2);
% Automatic sequence.
ax1 = nexttile;
plot(freq, autouncued_power_M1, '-b', 'LineWidth', 1.5); 
hold on;
h1 = patch([freq; flipud(freq)],...
    [(autouncued_power_M1 - std_autouncued_power_M1);...
    flipud(autouncued_power_M1 + std_autouncued_power_M1)], 'b',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, autocued_power_M1, '-m', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(autocued_power_M1 - std_autocued_power_M1);...
    flipud(autocued_power_M1 + std_autocued_power_M1)], 'm',...
    'LineStyle', 'none'); set(h2,'FaceAlpha',0.2);
hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
xlim([4 48]);
legend('Auto Uncued', '', 'Auto Cued', '');
xlabel('Frequency (Hz)');
ylabel('PSD');
% Non-automatic sequence.
ax2 = nexttile;
plot(freq, nonautouncued_power_M1, '-g', 'LineWidth', 1.5);
hold on;
h1 = patch([freq; flipud(freq)],...
    [(nonautouncued_power_M1 - std_nonautouncued_power_M1);...
    flipud(nonautouncued_power_M1 + std_nonautouncued_power_M1)], 'g',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, nonautocued_power_M1, '-r', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(nonautocued_power_M1 - std_nonautocued_power_M1);...
    flipud(nonautocued_power_M1 + std_nonautocued_power_M1)], 'r',...
    'LineStyle', 'none'); set(h2,'FaceAlpha',0.2);
hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
xlim([4 48]);
legend('Non-Auto Uncued', '', 'Non-Auto Cued', '');
xlabel('Frequency (Hz)');
ylabel('PSD');
linkaxes([ax1 ax2],'y')

set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, fullfile(results_path, 'M1_PowervsFreq'),'png');

%% PPC.
autouncued_power_P3 = autouncued_power(:, P3_loc);
autocued_power_P3 = autocued_power(:, P3_loc);
nonautouncued_power_P3 = nonautouncued_power(:, P3_loc);
nonautocued_power_P3 = nonautocued_power(:, P3_loc);
std_autouncued_power_P3 = std_autouncued_power(:, P3_loc);
std_autocued_power_P3 = std_autocued_power(:, P3_loc);
std_nonautouncued_power_P3 = std_nonautouncued_power(:, P3_loc);
std_nonautocued_power_P3 = std_nonautocued_power(:, P3_loc);

autouncued_power_P4 = autouncued_power(:, P4_loc);
autocued_power_P4 = autocued_power(:, P4_loc);
nonautouncued_power_P4 = nonautouncued_power(:, P4_loc);
nonautocued_power_P4 = nonautocued_power(:, P4_loc);
std_autouncued_power_P4 = std_autouncued_power(:, P4_loc);
std_autocued_power_P4 = std_autocued_power(:, P4_loc);
std_nonautouncued_power_P4 = std_nonautouncued_power(:, P4_loc);
std_nonautocued_power_P4 = std_nonautocued_power(:, P4_loc);

autouncued_power_PPC = mean([autouncued_power_P3 autouncued_power_P4], 2);
autocued_power_PPC = mean([autocued_power_P3 autocued_power_P4], 2);
nonautouncued_power_PPC = mean([nonautouncued_power_P3 nonautouncued_power_P4], 2);
nonautocued_power_PPC = mean([nonautocued_power_P3 nonautocued_power_P4], 2);
std_autouncued_power_PPC = mean([std_autouncued_power_P3 std_autouncued_power_P4], 2);
std_autocued_power_PPC = mean([std_autocued_power_P3 std_autocued_power_P4], 2);
std_nonautouncued_power_PPC = mean([std_nonautouncued_power_P3 std_nonautouncued_power_P4], 2);
std_nonautocued_power_PPC = mean([std_nonautocued_power_P3 std_nonautocued_power_P4], 2);

figure; title('PPC');
tiledlayout(1, 2);
% Automatic sequence.
ax1 = nexttile;
plot(freq, autouncued_power_PPC, '-b', 'LineWidth', 1.5); 
hold on;
h1 = patch([freq; flipud(freq)],...
    [(autouncued_power_PPC - std_autouncued_power_PPC);...
    flipud(autouncued_power_PPC + std_autouncued_power_PPC)], 'b',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, autocued_power_PPC, '-m', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(autocued_power_PPC - std_autocued_power_PPC);...
    flipud(autocued_power_PPC + std_autocued_power_PPC)], 'm',...
    'LineStyle', 'none'); set(h2,'FaceAlpha',0.2);
hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
xlim([4 48]);
legend('Auto Uncued', '', 'Auto Cued', '');
xlabel('Frequency (Hz)');
ylabel('PSD');
% Non-automatic sequence.
ax2 = nexttile;
plot(freq, nonautouncued_power_PPC, '-g', 'LineWidth', 1.5);
hold on;
h1 = patch([freq; flipud(freq)],...
    [(nonautouncued_power_PPC - std_nonautouncued_power_PPC);...
    flipud(nonautouncued_power_PPC + std_nonautouncued_power_PPC)], 'g',...
    'LineStyle', 'none'); set(h1,'FaceAlpha',0.2);
hold on;
plot(freq, nonautocued_power_PPC, '-r', 'LineWidth', 1.5);
hold on;
h2 = patch([freq; flipud(freq)],...
    [(nonautocued_power_PPC - std_nonautocued_power_PPC);...
    flipud(nonautocued_power_PPC + std_nonautocued_power_PPC)], 'r',...
    'LineStyle', 'none'); set(h2,'FaceAlpha',0.2);
hold on;
xline(8); hold on;
xline(13); hold on;
xline(32); hold off;
xlim([4 48]);
legend('Non-Auto Uncued', '', 'Non-Auto Cued', '');
xlabel('Frequency (Hz)');
ylabel('PSD');
linkaxes([ax1 ax2],'y')

set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, fullfile(results_path, 'PPC_PowervsFreq'),'png');

disp('This was the end of individual subjects.');
disp('These are the results for the average of all subjects.');

%% Save the values onto all subs struct.
avg.autouncued_power_DLPFC = autouncued_power_DLPFC;
avg.autocued_power_DLPFC = autocued_power_DLPFC;
avg.nonautouncued_power_DLPFC = nonautouncued_power_DLPFC;
avg.nonautocued_power_DLPFC = nonautocued_power_DLPFC;
avg.autouncued_power_SMA = autouncued_power_SMA;
avg.autocued_power_SMA = autocued_power_SMA;
avg.nonautouncued_power_SMA = nonautouncued_power_SMA;
avg.nonautocued_power_SMA = nonautocued_power_SMA;
avg.autouncued_power_M1 = autouncued_power_M1;
avg.autocued_power_M1 = autocued_power_M1;
avg.nonautouncued_power_M1 = nonautouncued_power_M1;
avg.nonautocued_power_M1 = nonautocued_power_M1;
avg.autouncued_power_PPC = autouncued_power_PPC;
avg.autocued_power_PPC = autocued_power_PPC;
avg.nonautouncued_power_PPC = nonautouncued_power_PPC;
avg.nonautocued_power_PPC = nonautocued_power_PPC;
allsubs.avg = avg;

std.autouncued_power_DLPFC = std_autouncued_power_DLPFC;
std.autocued_power_DLPFC = std_autocued_power_DLPFC;
std.nonautouncued_power_DLPFC = std_nonautouncued_power_DLPFC;
std.nonautocued_power_DLPFC = std_nonautocued_power_DLPFC;
std.autouncued_power_SMA = std_autouncued_power_SMA;
std.autocued_power_SMA = std_autocued_power_SMA;
std.nonautouncued_power_SMA = std_nonautouncued_power_SMA;
std.nonautocued_power_SMA = std_nonautocued_power_SMA;
std.autouncued_power_M1 = std_autouncued_power_M1;
std.autocued_power_M1 = std_autocued_power_M1;
std.nonautouncued_power_M1 = std_nonautouncued_power_M1;
std.nonautocued_power_M1 = std_nonautocued_power_M1;
std.autouncued_power_PPC = std_autouncued_power_PPC;
std.autocued_power_PPC = std_autocued_power_PPC;
std.nonautouncued_power_PPC = std_nonautouncued_power_PPC;
std.nonautocued_power_PPC = std_nonautocued_power_PPC;
allsubs.std = std;

% Save the struct from all subs.
save(strcat(results_path, '\psd_allsubs.mat'), 'allsubs')

%% Functions

% Loop through the power from the individual trials and average them.
function [power, freq] = calculateAveragePowerAllTrials(EEG, event_samp,...
    startTask, endTask, keypresses, th)

for trial=1:length(startTask)
    
    if trial==1
        size_power_allEpochs = 1;
    end
    
    % Get the keypresses within that trial.
    keypresses_trial = keypresses(keypresses > startTask(trial)...
        & keypresses < endTask(trial));
    keypresses_times = event_samp(keypresses_trial);
    
    % Epoch the data into the different keypresses.
    for epoch = 1:length(keypresses_times)
        
        EEG_epoch = pop_select(EEG, 'point',...
            [keypresses_times(epoch)-floor(0.4*EEG.srate)...
            keypresses_times(epoch)+ceil(0.4*EEG.srate)]);
        epoch_data = EEG_epoch.data;
        
        [power_oneEpoch, freq] = calculatePowerPerEpoch(EEG_epoch,...
            epoch_data, th);
        
        power_allEpochs(:, :, size_power_allEpochs) = power_oneEpoch;
        size_power_allEpochs = size_power_allEpochs+1;
        
    end
    
end

% Take the average of every epoch.
power = mean(power_allEpochs, 3, 'omitnan');

end

% From the trial data, calculate the power over the frequencies of the signal
% for all electrodes.
function [power, freq] =...
    calculatePowerPerEpoch(EEG_epoch, epoch_data, th)

% Using a Hann window.
window = 1:0.8*EEG_epoch.srate;
% Select the data of this specific window [channel x time].
data_window = epoch_data(:, window);

% Channel loop.
for channel = 1:size(data_window, 1)
    % Calculate PSD
    [P, f] = periodogram(data_window(channel, :),...
        hann(size(data_window, 2)),...
        2^(2 + nextpow2(size(data_window, 2))), EEG_epoch.srate);
    
    % Save the power for the frequencies of  the signal.
    pow(:, channel) = P((f(:,1)>=4 & f(:,1)<=48),1);
    
end

% Change frequency variable for frequencies of the signal.
freq = f(f(:,1)>=4 & f(:,1)<=48);

% For the bad channels, give NaN value.
power = pow;
for channel = 1:size(data_window, 1)
    if  pow(channel) > th(channel)
        power(:, channel) = NaN;
    end
end

end

% Add NaN in the lines where channels were removed during pre-processing.
function power_array_out = compensateRemovedChannels(power_array_in, EEG, list_channels)

% Array to see which channels are missing.
list_present = zeros(30, 1);

% If there are less than 30 channels.
if size(power_array_in, 1)~=30
    
    % Initialize new power array.
    power_array_out = zeros(size(power_array_in));
    power_array_out(:, 1:size(power_array_in, 2)) = power_array_in;
    % Get the channels present in the signal.
    channels_present = EEG.chanlocs;
    % Go through the lists of channels supposed to be present and channels
    % actually present and mark them as 1 if present and as 0 if not
    % present.
    for i=1:size(list_channels, 1)
        for j=1:size(channels_present, 2)
            if convertCharsToStrings(channels_present(j).labels) == list_channels(i)
                list_present(i) = 1;
                break;
            end
        end
    end
    
    numMissing=0;
    % For the non-present channels, put NaN value in them and update the
    % rest of the values.
    for k=1:size(list_present, 1)
        if list_present(k)==0
            numMissing = numMissing+1;
            power_array_out(:, k) = NaN;
            for x=k+1:size(power_array_in, 2)
                power_array_out(:, x)=power_array_in(:, x-numMissing);
            end
        end
    end
    power_array_out(:, 30:-1:30-numMissing+1) =...
        power_array_in(:, size(power_array_in,2):-1:size(power_array_in, 2)-numMissing+1);
else
    power_array_out = power_array_in;
end

end