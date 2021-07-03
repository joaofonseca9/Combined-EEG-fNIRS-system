%% Statistic for ERD/ERS.

clear; clc; close all;

laptop = 'laptopMariana';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);
results_path = 'C:\Users\maria\OneDrive\Ambiente de Trabalho\Automaticity Results\Topoplots';
statistics_path = 'C:\Users\maria\OneDrive\Ambiente de Trabalho\Statistics\Topoplots';

load(fullfile(results_path, 'erders_allsubs.mat'), 'allsubs');

subrec = ["28" "04"; "02" "02"; "76" "01"];

for subject=1:size(subrec, 1)
    sub = subrec(subject, 1);
    rec = subrec(subject, 2);
    
    % Subject 1 had no channels eliminated.
    if subject==1
        load([mainpath_in, '\pre-processed\sub-', char(sub), '\eeg\sub-',...
            char(sub), '_rec-', char(rec), '_eeg_divided.mat']);
    end
    
    % Get the locations of the channels of interest.
    locs = {EEG_divided.EEG_task.chanlocs.labels};
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

    sub_struct = allsubs.(genvarname(strcat('sub', char(sub))));
    
    %% Get average value over regions for each frequency band - automatic
    % sequence.
    
    % Theta.
    autouncued_theta = sub_struct.autouncued_ERD_ERS_theta;
    autocued_theta = sub_struct.autocued_ERD_ERS_theta;
    
    auto_theta_DLPFC(subject, 1) = mean([autouncued_theta(F7_loc)...
        autouncued_theta(F8_loc)], 'omitnan');
    auto_theta_DLPFC(subject, 2) = mean([autocued_theta(F7_loc)...
        autocued_theta(F8_loc)], 'omitnan');
    auto_theta_SMA(subject, 1) = mean([autouncued_theta(FC1_loc)...
        autouncued_theta(FC2_loc) autouncued_theta(Cz_loc)], 'omitnan');
    auto_theta_SMA(subject, 2) = mean([autocued_theta(FC1_loc)...
        autocued_theta(FC2_loc) autocued_theta(Cz_loc)], 'omitnan');
    auto_theta_M1(subject, 1) = autouncued_theta(C3_loc);
    auto_theta_M1(subject, 2) = autocued_theta(C3_loc);
    auto_theta_PPC(subject, 1) = mean([autouncued_theta(P3_loc)...
        autouncued_theta(P4_loc)], 'omitnan');
    auto_theta_PPC(subject, 2) = mean([autocued_theta(P3_loc)...
        autocued_theta(P4_loc)], 'omitnan');
    
    % Alpha.
    autouncued_alpha = sub_struct.autouncued_ERD_ERS_alpha;
    autocued_alpha = sub_struct.autocued_ERD_ERS_alpha;
    
    auto_alpha_DLPFC(subject, 1) = mean([autouncued_alpha(F7_loc)...
        autouncued_alpha(F8_loc)], 'omitnan');
    auto_alpha_DLPFC(subject, 2) = mean([autocued_alpha(F7_loc)...
        autocued_alpha(F8_loc)], 'omitnan');
    auto_alpha_SMA(subject, 1) = mean([autouncued_alpha(FC1_loc)...
        autouncued_alpha(FC2_loc) autouncued_alpha(Cz_loc)], 'omitnan');
    auto_alpha_SMA(subject, 2) = mean([autocued_alpha(FC1_loc)...
        autocued_alpha(FC2_loc) autocued_alpha(Cz_loc)], 'omitnan');
    auto_alpha_M1(subject, 1) = autouncued_alpha(C3_loc);
    auto_alpha_M1(subject, 2) = autocued_alpha(C3_loc);
    auto_alpha_PPC(subject, 1) = mean([autouncued_alpha(P3_loc)...
        autouncued_alpha(P4_loc)], 'omitnan');
    auto_alpha_PPC(subject, 2) = mean([autocued_alpha(P3_loc)...
        autocued_alpha(P4_loc)], 'omitnan');
    
    % Beta.
    autouncued_beta = sub_struct.autouncued_ERD_ERS_beta;
    autocued_beta = sub_struct.autocued_ERD_ERS_beta;
    
    auto_beta_DLPFC(subject, 1) = mean([autouncued_beta(F7_loc)...
        autouncued_beta(F8_loc)], 'omitnan');
    auto_beta_DLPFC(subject, 2) = mean([autocued_beta(F7_loc)...
        autocued_beta(F8_loc)], 'omitnan');
    auto_beta_SMA(subject, 1) = mean([autouncued_beta(FC1_loc)...
        autouncued_beta(FC2_loc) autouncued_beta(Cz_loc)], 'omitnan');
    auto_beta_SMA(subject, 2) = mean([autocued_beta(FC1_loc)...
        autocued_beta(FC2_loc) autocued_beta(Cz_loc)], 'omitnan');
    auto_beta_M1(subject, 1) = autouncued_beta(C3_loc);
    auto_beta_M1(subject, 2) = autocued_beta(C3_loc);
    auto_beta_PPC(subject, 1) = mean([autouncued_beta(P3_loc)...
        autouncued_beta(P4_loc)], 'omitnan');
    auto_beta_PPC(subject, 2) = mean([autocued_beta(P3_loc)...
        autocued_beta(P4_loc)], 'omitnan');
    
    % Gamma.
    autouncued_gamma = sub_struct.autouncued_ERD_ERS_gamma;
    autocued_gamma = sub_struct.autocued_ERD_ERS_gamma;
    
    auto_gamma_DLPFC(subject, 1) = mean([autouncued_gamma(F7_loc)...
        autouncued_gamma(F8_loc)], 'omitnan');
    auto_gamma_DLPFC(subject, 2) = mean([autocued_gamma(F7_loc)...
        autocued_gamma(F8_loc)], 'omitnan');
    auto_gamma_SMA(subject, 1) = mean([autouncued_gamma(FC1_loc)...
        autouncued_gamma(FC2_loc) autouncued_gamma(Cz_loc)], 'omitnan');
    auto_gamma_SMA(subject, 2) = mean([autocued_gamma(FC1_loc)...
        autocued_gamma(FC2_loc) autocued_gamma(Cz_loc)], 'omitnan');
    auto_gamma_M1(subject, 1) = autouncued_gamma(C3_loc);
    auto_gamma_M1(subject, 2) = autocued_gamma(C3_loc);
    auto_gamma_PPC(subject, 1) = mean([autouncued_gamma(P3_loc)...
        autouncued_gamma(P4_loc)], 'omitnan');
    auto_gamma_PPC(subject, 2) = mean([autocued_gamma(P3_loc)...
        autocued_gamma(P4_loc)], 'omitnan');
    
    %% Get average value over regions for each frequency band - non-automatic
    % sequence.
    
    % Theta.
    nonautouncued_theta = sub_struct.nonautouncued_ERD_ERS_theta;
    nonautocued_theta = sub_struct.nonautocued_ERD_ERS_theta;
    
    nonauto_theta_DLPFC(subject, 1) = mean([nonautouncued_theta(F7_loc)...
        nonautouncued_theta(F8_loc)], 'omitnan');
    nonauto_theta_DLPFC(subject, 2) = mean([nonautocued_theta(F7_loc)...
        nonautocued_theta(F8_loc)], 'omitnan');
    nonauto_theta_SMA(subject, 1) = mean([nonautouncued_theta(FC1_loc)...
        nonautouncued_theta(FC2_loc) nonautouncued_theta(Cz_loc)], 'omitnan');
    nonauto_theta_SMA(subject, 2) = mean([nonautocued_theta(FC1_loc)...
        nonautocued_theta(FC2_loc) nonautocued_theta(Cz_loc)], 'omitnan');
    nonauto_theta_M1(subject, 1) = nonautouncued_theta(C3_loc);
    nonauto_theta_M1(subject, 2) = nonautocued_theta(C3_loc);
    nonauto_theta_PPC(subject, 1) = mean([nonautouncued_theta(P3_loc)...
        nonautouncued_theta(P4_loc)], 'omitnan');
    nonauto_theta_PPC(subject, 2) = mean([nonautocued_theta(P3_loc)...
        nonautocued_theta(P4_loc)], 'omitnan');
    
    % Alpha.
    nonautouncued_alpha = sub_struct.nonautouncued_ERD_ERS_alpha;
    nonautocued_alpha = sub_struct.nonautocued_ERD_ERS_alpha;
    
    nonauto_alpha_DLPFC(subject, 1) = mean([nonautouncued_alpha(F7_loc)...
        nonautouncued_alpha(F8_loc)], 'omitnan');
    nonauto_alpha_DLPFC(subject, 2) = mean([nonautocued_alpha(F7_loc)...
        nonautocued_alpha(F8_loc)], 'omitnan');
    nonauto_alpha_SMA(subject, 1) = mean([nonautouncued_alpha(FC1_loc)...
        nonautouncued_alpha(FC2_loc) nonautouncued_alpha(Cz_loc)], 'omitnan');
    nonauto_alpha_SMA(subject, 2) = mean([nonautocued_alpha(FC1_loc)...
        nonautocued_alpha(FC2_loc) nonautocued_alpha(Cz_loc)], 'omitnan');
    nonauto_alpha_M1(subject, 1) = nonautouncued_alpha(C3_loc);
    nonauto_alpha_M1(subject, 2) = nonautocued_alpha(C3_loc);
    nonauto_alpha_PPC(subject, 1) = mean([nonautouncued_alpha(P3_loc)...
        nonautouncued_alpha(P4_loc)], 'omitnan');
    nonauto_alpha_PPC(subject, 2) = mean([nonautocued_alpha(P3_loc)...
        nonautocued_alpha(P4_loc)], 'omitnan');
    
    % Beta.
    nonautouncued_beta = sub_struct.nonautouncued_ERD_ERS_beta;
    nonautocued_beta = sub_struct.nonautocued_ERD_ERS_beta;
    
    nonauto_beta_DLPFC(subject, 1) = mean([nonautouncued_beta(F7_loc)...
        nonautouncued_beta(F8_loc)], 'omitnan');
    nonauto_beta_DLPFC(subject, 2) = mean([nonautocued_beta(F7_loc)...
        nonautocued_beta(F8_loc)], 'omitnan');
    nonauto_beta_SMA(subject, 1) = mean([nonautouncued_beta(FC1_loc)...
        nonautouncued_beta(FC2_loc) nonautouncued_beta(Cz_loc)], 'omitnan');
    nonauto_beta_SMA(subject, 2) = mean([nonautocued_beta(FC1_loc)...
        nonautocued_beta(FC2_loc) nonautocued_beta(Cz_loc)], 'omitnan');
    nonauto_beta_M1(subject, 1) = nonautouncued_beta(C3_loc);
    nonauto_beta_M1(subject, 2) = nonautocued_beta(C3_loc);
    nonauto_beta_PPC(subject, 1) = mean([nonautouncued_beta(P3_loc)...
        nonautouncued_beta(P4_loc)], 'omitnan');
    nonauto_beta_PPC(subject, 2) = mean([nonautocued_beta(P3_loc)...
        nonautocued_beta(P4_loc)], 'omitnan');
    
    % Gamma.
    nonautouncued_gamma = sub_struct.nonautouncued_ERD_ERS_gamma;
    nonautocued_gamma = sub_struct.nonautocued_ERD_ERS_gamma;
    
    nonauto_gamma_DLPFC(subject, 1) = mean([nonautouncued_gamma(F7_loc)...
        nonautouncued_gamma(F8_loc)], 'omitnan');
    nonauto_gamma_DLPFC(subject, 2) = mean([nonautocued_gamma(F7_loc)...
        nonautocued_gamma(F8_loc)], 'omitnan');
    nonauto_gamma_SMA(subject, 1) = mean([nonautouncued_gamma(FC1_loc)...
        nonautouncued_gamma(FC2_loc) nonautouncued_gamma(Cz_loc)], 'omitnan');
    nonauto_gamma_SMA(subject, 2) = mean([nonautocued_gamma(FC1_loc)...
        nonautocued_gamma(FC2_loc) nonautocued_gamma(Cz_loc)], 'omitnan');
    nonauto_gamma_M1(subject, 1) = nonautouncued_gamma(C3_loc);
    nonauto_gamma_M1(subject, 2) = nonautocued_gamma(C3_loc);
    nonauto_gamma_PPC(subject, 1) = mean([nonautouncued_gamma(P3_loc)...
        nonautouncued_gamma(P4_loc)], 'omitnan');
    nonauto_gamma_PPC(subject, 2) = mean([nonautocued_gamma(P3_loc)...
        nonautocued_gamma(P4_loc)], 'omitnan');
    
end

%% Automatic sequence.
% Wilcoxon signed-rank test if not normally distributed.
% T-test if normally distributed.

% Theta band.

h = kstest(auto_theta_DLPFC)
if h==0 % Normally distributed data.
    [h, p] = ttest(auto_theta_DLPFC(:, 1), auto_theta_DLPFC(:, 2));
    stats_auto.auto_theta_DLPFC = p;
else 
    p = signrank(auto_theta_DLPFC(:, 1), auto_theta_DLPFC(:, 2));
    stats_auto.auto_theta_DLPFC = p;
end

h = kstest(auto_theta_SMA)
if h==0 % Normally distributed data.
    [h, p] = ttest(auto_theta_SMA(:, 1), auto_theta_SMA(:, 2));
    stats_auto.auto_theta_SMA = p;
else
    p = signrank(auto_theta_SMA(:, 1), auto_theta_SMA(:, 2));
    stats_auto.auto_theta_SMA = p;
end

h = kstest(auto_theta_M1)
if h==0 % Normally distributed data.
    [h, p] = ttest(auto_theta_M1(:, 1), auto_theta_M1(:, 2));
    stats_auto.auto_theta_M1 = p;
else
    p = signrank(auto_theta_M1(:, 1), auto_theta_M1(:, 2));
    stats_auto.auto_theta_M1 = p;
end

h = kstest(auto_theta_PPC)
if h==0 % Normally distributed data.
    [h, p] = ttest(auto_theta_PPC(:, 1), auto_theta_PPC(:, 2));
    stats_auto.auto_theta_PPC = p;
else
    p = signrank(auto_theta_PPC(:, 1), auto_theta_PPC(:, 2));
    stats_auto.auto_theta_PPC = p;
end

%%

% Alpha band.

h = kstest(auto_alpha_DLPFC)
if h==0 % Normally distributed data.
    [h, p] = ttest(auto_alpha_DLPFC(:, 1), auto_alpha_DLPFC(:, 2));
    stats_auto.auto_alpha_DLPFC = p;
else 
    p = signrank(auto_alpha_DLPFC(:, 1), auto_alpha_DLPFC(:, 2));
    stats_auto.auto_alpha_DLPFC = p;
end

h = kstest(auto_alpha_SMA)
if h==0 % Normally distributed data.
    [h, p] = ttest(auto_alpha_SMA(:, 1), auto_alpha_SMA(:, 2));
    stats_auto.auto_alpha_SMA = p;
else
    p = signrank(auto_alpha_SMA(:, 1), auto_alpha_SMA(:, 2));
    stats_auto.auto_alpha_SMA = p;
end

h = kstest(auto_alpha_M1)
if h==0 % Normally distributed data.
    [h, p] = ttest(auto_alpha_M1(:, 1), auto_alpha_M1(:, 2));
    stats_auto.auto_alpha_M1 = p;
else
    p = signrank(auto_alpha_M1(:, 1), auto_alpha_M1(:, 2));
    stats_auto.auto_alpha_M1 = p;
end

h = kstest(auto_alpha_PPC)
if h==0 % Normally distributed data.
    [h, p] = ttest(auto_alpha_PPC(:, 1), auto_alpha_PPC(:, 2));
    stats_auto.auto_alpha_PPC = p;
else
    p = signrank(auto_alpha_PPC(:, 1), auto_alpha_PPC(:, 2));
    stats_auto.auto_alpha_PPC = p;
end

%%

% Beta band.

h = kstest(auto_beta_DLPFC)
if h==0 % Normally distributed data.
    [h, p] = ttest(auto_beta_DLPFC(:, 1), auto_beta_DLPFC(:, 2));
    stats_auto.auto_beta_DLPFC = p;
else 
    p = signrank(auto_beta_DLPFC(:, 1), auto_beta_DLPFC(:, 2));
    stats_auto.auto_beta_DLPFC = p;
end

h = kstest(auto_beta_SMA)
if h==0 % Normally distributed data.
    [h, p] = ttest(auto_beta_SMA(:, 1), auto_beta_SMA(:, 2));
    stats_auto.auto_beta_SMA = p;
else
    p = signrank(auto_beta_SMA(:, 1), auto_beta_SMA(:, 2));
    stats_auto.auto_beta_SMA = p;
end

h = kstest(auto_beta_M1)
if h==0 % Normally distributed data.
    [h, p] = ttest(auto_beta_M1(:, 1), auto_beta_M1(:, 2));
    stats_auto.auto_beta_M1 = p;
else
    p = signrank(auto_beta_M1(:, 1), auto_beta_M1(:, 2));
    stats_auto.auto_beta_M1 = p;
end

h = kstest(auto_beta_PPC)
if h==0 % Normally distributed data.
    [h, p] = ttest(auto_beta_PPC(:, 1), auto_beta_PPC(:, 2));
    stats_auto.auto_beta_PPC = p;
else
    p = signrank(auto_beta_PPC(:, 1), auto_beta_PPC(:, 2));
    stats_auto.auto_beta_PPC = p;
end

%% 

% Gamma band.

h = kstest(auto_gamma_DLPFC)
if h==0 % Normally distributed data.
    [h, p] = ttest(auto_gamma_DLPFC(:, 1), auto_gamma_DLPFC(:, 2));
    stats_auto.auto_gamma_DLPFC = p;
else 
    p = signrank(auto_gamma_DLPFC(:, 1), auto_gamma_DLPFC(:, 2));
    stats_auto.auto_gamma_DLPFC = p;
end

h = kstest(auto_gamma_SMA)
if h==0 % Normally distributed data.
    [h, p] = ttest(auto_gamma_SMA(:, 1), auto_gamma_SMA(:, 2));
    stats_auto.auto_gamma_SMA = p;
else
    p = signrank(auto_gamma_SMA(:, 1), auto_gamma_SMA(:, 2));
    stats_auto.auto_gamma_SMA = p;
end

h = kstest(auto_gamma_M1)
if h==0 % Normally distributed data.
    [h, p] = ttest(auto_gamma_M1(:, 1), auto_gamma_M1(:, 2));
    stats_auto.auto_gamma_M1 = p;
else
    p = signrank(auto_gamma_M1(:, 1), auto_gamma_M1(:, 2));
    stats_auto.auto_gamma_M1 = p;
end

h = kstest(auto_gamma_PPC)
if h==0 % Normally distributed data.
    [h, p] = ttest(auto_gamma_PPC(:, 1), auto_gamma_PPC(:, 2));
    stats_auto.auto_gamma_PPC = p;
else
    p = signrank(auto_gamma_PPC(:, 1), auto_gamma_PPC(:, 2));
    stats_auto.auto_gamma_PPC = p;
end

stats.stats_auto = stats_auto;

%% Non-automatic sequence.
% Wilcoxon signed-rank test if not normally distributed.
% T-test if normally distributed.

% Theta band.

h = kstest(nonauto_theta_DLPFC)
if h==0 % Normally distributed data.
    [h, p] = ttest(nonauto_theta_DLPFC(:, 1), nonauto_theta_DLPFC(:, 2));
    stats_nonauto.nonauto_theta_DLPFC = p;
else 
    p = signrank(nonauto_theta_DLPFC(:, 1), nonauto_theta_DLPFC(:, 2));
    stats_nonauto.nonauto_theta_DLPFC = p;
end

h = kstest(nonauto_theta_SMA)
if h==0 % Normally distributed data.
    [h, p] = ttest(nonauto_theta_SMA(:, 1), nonauto_theta_SMA(:, 2));
    stats_nonauto.nonauto_theta_SMA = p;
else
    p = signrank(nonauto_theta_SMA(:, 1), nonauto_theta_SMA(:, 2));
    stats_nonauto.nonauto_theta_SMA = p;
end

h = kstest(nonauto_theta_M1)
if h==0 % Normally distributed data.
    [h, p] = ttest(nonauto_theta_M1(:, 1), nonauto_theta_M1(:, 2));
    stats_nonauto.nonauto_theta_M1 = p;
else
    p = signrank(nonauto_theta_M1(:, 1), nonauto_theta_M1(:, 2));
    stats_nonauto.nonauto_theta_M1 = p;
end

h = kstest(nonauto_theta_PPC)
if h==0 % Normally distributed data.
    [h, p] = ttest(nonauto_theta_PPC(:, 1), nonauto_theta_PPC(:, 2));
    stats_nonauto.nonauto_theta_PPC = p;
else
    p = signrank(nonauto_theta_PPC(:, 1), nonauto_theta_PPC(:, 2));
    stats_nonauto.nonauto_theta_PPC = p;
end

%%

% Alpha band.

h = kstest(nonauto_alpha_DLPFC)
if h==0 % Normally distributed data.
    [h, p] = ttest(nonauto_alpha_DLPFC(:, 1), nonauto_alpha_DLPFC(:, 2));
    stats_nonauto.nonauto_alpha_DLPFC = p;
else 
    p = signrank(nonauto_alpha_DLPFC(:, 1), nonauto_alpha_DLPFC(:, 2));
    stats_nonauto.nonauto_alpha_DLPFC = p;
end

h = kstest(nonauto_alpha_SMA)
if h==0 % Normally distributed data.
    [h, p] = ttest(nonauto_alpha_SMA(:, 1), nonauto_alpha_SMA(:, 2));
    stats_nonauto.nonauto_alpha_SMA = p;
else
    p = signrank(nonauto_alpha_SMA(:, 1), nonauto_alpha_SMA(:, 2));
    stats_nonauto.nonauto_alpha_SMA = p;
end

h = kstest(nonauto_alpha_M1)
if h==0 % Normally distributed data.
    [h, p] = ttest(nonauto_alpha_M1(:, 1), nonauto_alpha_M1(:, 2));
    stats_nonauto.nonauto_alpha_M1 = p;
else
    p = signrank(nonauto_alpha_M1(:, 1), nonauto_alpha_M1(:, 2));
    stats_nonauto.nonauto_alpha_M1 = p;
end

h = kstest(nonauto_alpha_PPC)
if h==0 % Normally distributed data.
    [h, p] = ttest(nonauto_alpha_PPC(:, 1), nonauto_alpha_PPC(:, 2));
    stats_nonauto.nonauto_alpha_PPC = p;
else
    p = signrank(nonauto_alpha_PPC(:, 1), nonauto_alpha_PPC(:, 2));
    stats_nonauto.nonauto_alpha_PPC = p;
end

%%

% Beta band.

h = kstest(nonauto_beta_DLPFC)
if h==0 % Normally distributed data.
    [h, p] = ttest(nonauto_beta_DLPFC(:, 1), nonauto_beta_DLPFC(:, 2));
    stats_nonauto.nonauto_beta_DLPFC = p;
else 
    p = signrank(nonauto_beta_DLPFC(:, 1), nonauto_beta_DLPFC(:, 2));
    stats_nonauto.nonauto_beta_DLPFC = p;
end

h = kstest(nonauto_beta_SMA)
if h==0 % Normally distributed data.
    [h, p] = ttest(nonauto_beta_SMA(:, 1), nonauto_beta_SMA(:, 2));
    stats_nonauto.nonauto_beta_SMA = p;
else
    p = signrank(nonauto_beta_SMA(:, 1), nonauto_beta_SMA(:, 2));
    stats_nonauto.nonauto_beta_SMA = p;
end

h = kstest(nonauto_beta_M1)
if h==0 % Normally distributed data.
    [h, p] = ttest(nonauto_beta_M1(:, 1), nonauto_beta_M1(:, 2));
    stats_nonauto.nonauto_beta_M1 = p;
else
    p = signrank(nonauto_beta_M1(:, 1), nonauto_beta_M1(:, 2));
    stats_nonauto.nonauto_beta_M1 = p;
end

h = kstest(nonauto_beta_PPC)
if h==0 % Normally distributed data.
    [h, p] = ttest(nonauto_beta_PPC(:, 1), nonauto_beta_PPC(:, 2));
    stats_nonauto.nonauto_beta_PPC = p;
else
    p = signrank(nonauto_beta_PPC(:, 1), nonauto_beta_PPC(:, 2));
    stats_nonauto.nonauto_beta_PPC = p;
end

%% 

% Gamma band.

h = kstest(nonauto_gamma_DLPFC)
if h==0 % Normally distributed data.
    [h, p] = ttest(nonauto_gamma_DLPFC(:, 1), nonauto_gamma_DLPFC(:, 2));
    stats_nonauto.nonauto_gamma_DLPFC = p;
else 
    p = signrank(nonauto_gamma_DLPFC(:, 1), nonauto_gamma_DLPFC(:, 2));
    stats_nonauto.nonauto_gamma_DLPFC = p;
end

h = kstest(nonauto_gamma_SMA)
if h==0 % Normally distributed data.
    [h, p] = ttest(nonauto_gamma_SMA(:, 1), nonauto_gamma_SMA(:, 2));
    stats_nonauto.nonauto_gamma_SMA = p;
else
    p = signrank(nonauto_gamma_SMA(:, 1), nonauto_gamma_SMA(:, 2));
    stats_nonauto.nonauto_gamma_SMA = p;
end

h = kstest(nonauto_gamma_M1)
if h==0 % Normally distributed data.
    [h, p] = ttest(nonauto_gamma_M1(:, 1), nonauto_gamma_M1(:, 2));
    stats_nonauto.nonauto_gamma_M1 = p;
else
    p = signrank(nonauto_gamma_M1(:, 1), nonauto_gamma_M1(:, 2));
    stats_nonauto.nonauto_gamma_M1 = p;
end

h = kstest(nonauto_gamma_PPC)
if h==0 % Normally distributed data.
    [h, p] = ttest(nonauto_gamma_PPC(:, 1), nonauto_gamma_PPC(:, 2));
    stats_nonauto.nonauto_gamma_PPC = p;
else
    p = signrank(nonauto_gamma_PPC(:, 1), nonauto_gamma_PPC(:, 2));
    stats_nonauto.nonauto_gamma_PPC = p;
end

stats.stats_nonauto = stats_nonauto;

save(strcat(statistics_path, '\stats_erders.mat'), 'stats');