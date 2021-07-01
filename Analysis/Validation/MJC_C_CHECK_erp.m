%% VAL_Analysis_comb

close all; clear all; clc;

%% Settings
addpath (fullfile(pwd,'..')) %add the path with analysis scripts
% add path to correct folders and open eeglab
laptop = 'laptopJoao';
[mainpath_in, ~, eeglab_path] = addFolders(laptop);

mainpath_out = 'C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\Experiment\Combined-EEG-fNIRS-system\Analysis\Validation';
if ~isfolder(fullfile(mainpath_out,'Fig_CHECK_erp'))
    mkdir(fullfile(mainpath_out,'Fig_CHECK_erp'));
end

mainpath_outfig=fullfile(mainpath_out,'Fig_CHECK_erp');
% select ID number and cap
subject=[{'02','28','64'}];
subject_eeg_only=[{'03','04','10'}];%,'11','12',}];

%% General variables

chans       = {'Fp1';'Fpz';'Fp2';'F7';'F3';'AFFz';'F4';'F8';'FC5';'FC1';'FC2';'FC6';'T7';'C3';'Cz';'C4';'T8';'CP5';'CP1';'CP2';'CP6';'P7';'P3';'Pz';'P4';'P8';'POz';'O1';'Oz';'O2'};
time        = -0.1 : 1/1024 : 0.4;

fig_erp     = fullfile(mainpath_outfig,'ALL_VEP.jpg');
fig_topo    = fullfile(mainpath_outfig,'ALL_VEPtopo.jpg');

%% CHECK - ERP

for iSub = 1:length(subject)
    %% Load and collect subject data
    
    % load data
    ID = char(subject(iSub));
    load([mainpath_in,'\processed\sub-',ID,'\ResultsMatrix\sub-',ID,'_checkVEP.mat']); 
    comb = CHECKerp; 
    clear CHECKerp;
    
    ID_eeg_only= char(subject_eeg_only(iSub));
    load(['C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\Experiment\Data\Validation Data\EEG\ID',ID_eeg_only,'checkVEP_eegonly.mat']); 
    eeg_only = CHECKerp; 
    clear CHECKerp;
    
    % collect VEP of all subjects [chan x time x subject]
    for iChan = 1:length(chans)
        gelloc = find(strcmp(chans(iChan), eeg_only.chans));
        if ~isempty(gelloc);        allVEP.eeg_only(iChan,:,iSub) = eeg_only.ERP(gelloc,:);
        elseif isempty(gelloc);     allVEP.eeg_only(iChan,:,iSub) = NaN;
        end
        dryloc = find(strcmp(chans(iChan), comb.chans));
        if ~isempty(dryloc);        allVEP.comb(iChan,:,iSub) = comb.ERP(dryloc,:);
        elseif isempty(dryloc);     allVEP.comb(iChan,:,iSub) = NaN;
        end
        clear eeg_onlyloc combloc
    end
    
    % collect GFPt of all subjects [subject x time]
    allGFP.eeg_only(iSub,:) = eeg_only.GFPt;
    allGFP.comb(iSub,:) = comb.GFPt;
    
    
    
    clear D Dg Dd comb eeg_only
end


%% plot mean per channel over subjects + mean and std of GFPt

eeg_onlyGFP1 = mean(allGFP.eeg_only,1) + 2*std(allGFP.eeg_only,[],1);
eeg_onlyGFP2 = mean(allGFP.eeg_only,1) - 2*std(allGFP.eeg_only,[],1);
combGFP1 = mean(allGFP.comb,1) + 2*std(allGFP.comb,[],1);
combGFP2 = mean(allGFP.comb,1) - 2*std(allGFP.comb,[],1);
figure; 
subplot(2,2,1); plot(time,mean(allVEP.eeg_only,3,'omitnan'),'b');
subplot(2,2,2); plot(time,mean(allVEP.comb,3,'omitnan'),'r');
subplot(2,2,3); hold on; h1 = fill([time,fliplr(time)], [eeg_onlyGFP1,fliplr(eeg_onlyGFP2)],'b','LineStyle','none');
set(h1,'FaceAlpha',0.4); plot(time,mean(allGFP.eeg_only,1),'b','LineWidth',1.5); 
subplot(2,2,4); hold on; h2 = fill([time,fliplr(time)], [combGFP1,fliplr(combGFP2)],'r','LineStyle','none');
set(h2,'FaceAlpha',0.4); plot(time,mean(allGFP.comb,1),'r','LineWidth',1.5); 

subplot(2,2,1); set(gca,'FontSize',11); box on;
ylabel('Potential (\muV)','FontSize',14); title('eeg only cap','FontSize',14); ylim([-6 6])
xticks([-0.1:0.1:0.4]); yticks([-5:5:5]); line([0 0],[-6 6],'Color','k'); 
subplot(2,2,2); set(gca,'FontSize',11); box on;
ylabel('Potential (\muV)','FontSize',14); title('combined cap','FontSize',14); ylim([-6 6])
xticks([-0.1:0.1:0.4]); yticks([-5:5:5]); line([0 0],[-6 6],'Color','k'); 
subplot(2,2,3); set(gca,'FontSize',11); box on;
ylabel('GFPt (\muV)','FontSize',14); xlabel('Time (s)','FontSize',14); ylim([-2 26])
xticks([-0.1:0.1:0.4]); yticks([0:10:20]); line([0 0],[-2 25],'Color','k'); 
subplot(2,2,4); set(gca,'FontSize',11); box on;
ylabel('GFPt (\muV)','FontSize',14); xlabel('Time (s)','FontSize',14); ylim([-2 26])
xticks([-0.1:0.1:0.4]); yticks([0:10:20]); line([0 0],[-2 25],'Color','k'); 

clear eeg_onlyGFP1 eeg_onlyGFP2 combGFP1 combGFP2 h1 h2

saveas(gcf,fig_erp);

%% find main VEP components 
% find the main VEP components in the average GFPt 

% average GFPt and VEP over subjects
avgGFP.eeg_only = mean(allGFP.eeg_only,1);
avgGFP.comb = mean(allGFP.comb,1);
avgVEP.eeg_only = mean(allVEP.eeg_only,3, 'omitnan');
avgVEP.comb = mean(allVEP.comb,3, 'omitnan');

% get main component N1: 0.11-0.13s (N1 peak)
[eeg_onlyN1.loc, eeg_onlyN1.time, eeg_onlyN1.amp] = get_mainVEPcomponents (time, avgGFP.eeg_only, [0.11 0.13]);
[combN1.loc, combN1.time, combN1.amp] = get_mainVEPcomponents (time, avgGFP.comb, [0.11 0.13]);

% get main component P1: 0.14-0.16s (P100 peak)
[eeg_onlyP1.loc, eeg_onlyP1.time, eeg_onlyP1.amp] = get_mainVEPcomponents (time, avgGFP.eeg_only, [0.14 0.16]);
[combP1.loc, combP1.time, combP1.amp] = get_mainVEPcomponents (time, avgGFP.comb, [0.14 0.16]);

%% standardize VEP data within subject per cap

% standardize (x-meanX)/sdX
allVEPstd.eeg_only = (allVEP.eeg_only - mean(allVEP.eeg_only,1,'omitnan')) ./ std(allVEP.eeg_only,[],1,'omitnan');
allVEPstd.comb = (allVEP.comb - mean(allVEP.comb,1,'omitnan')) ./ std(allVEP.comb,[],1,'omitnan');
% > average standardized VEP
avgVEPstd.eeg_only = mean(allVEPstd.eeg_only,3, 'omitnan');
avgVEPstd.comb = mean(allVEPstd.comb,3, 'omitnan');

% normalize data 
allVEPnrm.eeg_only = (allVEP.eeg_only - min(allVEP.eeg_only,1,'omitnan')) ./ (max(allVEP.eeg_only,1,'omitnan') - min(allVEP.eeg_only,1,'omitnan'));
allVEPnrm.comb = (allVEP.comb - min(allVEP.comb,1,'omitnan')) ./ (max(allVEP.comb,1,'omitnan') - min(allVEP.comb,1,'omitnan'));
% > average normalized VEP
avgVEPnrm.eeg_only = mean(allVEPnrm.eeg_only,3, 'omitnan');
avgVEPnrm.comb = mean(allVEPnrm.comb,3, 'omitnan');

%% Root Mean Square Deviation (RMSD) - Fiedler (2015)
for iL = 1:length(avgGFP.comb)
    D(iL) = (avgGFP.eeg_only(iL)-avgGFP.comb(iL)).^2;
end
RMSD = sqrt(sum(D)/length(avgGFP.comb));    

%% Spearman's rank Correlation (COR) - Fiedler (2015)
for iL = 1:length(avgGFP.comb)
    Dg(iL) = avgGFP.eeg_only(iL)-mean(avgGFP.eeg_only);
    Dd(iL) = avgGFP.comb(iL)-mean(avgGFP.comb);
end
COR = sum(Dg.*Dd) / sqrt( sum(Dg.^2).*sum(Dd.^2) ); 

    % %% test for normality of COR and RMSD
% 
% [~,NormalDistr.COR]   = kstest(COR);
% [~,NormalDistr.RMSD]  = kstest(RMSD);
% 
% %% average RMSD and COR
% 
% RES.RMSD(1,1)   = mean(RMSD);
% RES.RMSD(1,2)   = 2*std(RMSD);
% RES.COR(1,1)    = mean(COR);
% RES.COR(1,2)    = 2*std(COR);

%% data for topoplot of main VEP components, average over subjects

% raw: get topoplot data of N1 and P1 peak
eeg_onlyN1.topo = avgVEP.eeg_only(:,eeg_onlyN1.loc);        eeg_onlyP1.topo = avgVEP.eeg_only(:,eeg_onlyP1.loc);
combN1.topo = avgVEP.comb(:,combN1.loc);        combP1.topo = avgVEP.comb(:,combP1.loc);

% std: get topoplot data of N1 and P1 standardized within subject and cap
eeg_onlyN1std.topo = avgVEPstd.eeg_only(:,eeg_onlyN1.loc);    eeg_onlyP1std.topo = avgVEPstd.eeg_only(:,eeg_onlyP1.loc);
combN1std.topo = avgVEPstd.comb(:,combN1.loc);    combP1std.topo = avgVEPstd.comb(:,combP1.loc);

% norm: get topoplot data of N1 and P1 normalized within subject and cap
eeg_onlyN1nrm.topo = avgVEPnrm.eeg_only(:,eeg_onlyN1.loc);    eeg_onlyP1nrm.topo = avgVEPnrm.eeg_only(:,eeg_onlyP1.loc);
combN1nrm.topo = avgVEPnrm.comb(:,combN1.loc);    combP1nrm.topo = avgVEPnrm.comb(:,combP1.loc);

%% permutation testing: raw data

load('chanlocs.mat')

% eeg_only-N1 VS comb-N1
% use topoplot to obtain grid (pixels) per subject 
for iSub = 1:size(subject,2)
    figure(100); subplot(121); title('combined cap'); [~,grid_comb(:,:,iSub),~,xmesh,ymesh] = topoplot(allVEP.comb(:,combN1.loc,iSub), chanlocs);
    figure(100); subplot(122); title('eeg only cap'); [~,grid_eeg_only(:,:,iSub),~,~,~] = topoplot(allVEP.eeg_only(:,eeg_onlyN1.loc,iSub), chanlocs);
    close Figure 100
end
% permutation testing
% figure, title('N1 peak - Raw Data')
[permN1.zmap, permN1.zmapTH, permN1.zmapTHmcc] = Permutation_VEPcomponent (grid_comb, grid_eeg_only, xmesh, ymesh);

% percentage of significant pixels 
nrpix_total     = numel(grid_comb(:,:,1)) - sum(isnan(grid_comb(:,:,1)),'all');
nrpix_sig       = sum(permN1.zmapTHmcc==1,'all');
permN1.sigPix   = 100*nrpix_sig / nrpix_total; 

clear grid_comb grid_eeg_only nrpix_sig nrix_total

% eeg_only-P1 VS comb-P1
% use topoplot to obtain grid (pixels) per subject 
for iSub = 1:size(subject,2)
    figure(100); subplot(121); title('combined cap'); [~,grid_comb(:,:,iSub),~,~,~] = topoplot(allVEP.comb(:,combP1.loc,iSub), chanlocs);
    figure(100); subplot(122); title('eeg only cap'); [~,grid_eeg_only(:,:,iSub),~,~,~] = topoplot(allVEP.eeg_only(:,eeg_onlyP1.loc,iSub), chanlocs);
    close Figure 100
end

% figure, title('P1 peak - Raw Data')
% permutation testing
[permP1.zmap, permP1.zmapTH, permP1.zmapTHmcc] = Permutation_VEPcomponent (grid_comb, grid_eeg_only, xmesh, ymesh);

% percentage of significant pixels 
nrpix_total     = numel(grid_comb(:,:,1)) - sum(isnan(grid_comb(:,:,1)),'all');
nrpix_sig       = sum(permP1.zmapTHmcc==1,'all');
permP1.sigPix   = 100*nrpix_sig / nrpix_total;

clear grid_comb grid_eeg_only nrpix_sig nrix_total

% plot topoplot with significant regions
PermTopoplot (permN1, permP1, eeg_onlyN1, combN1, eeg_onlyP1, combP1, chanlocs, xmesh, ymesh, 15, 'raw');
saveas(gcf, 'Fig_CHECK_erp/Topoplots_raw_data.jpg');

%% permutation testing: standardized data

% eeg_only-N1 VS comb-N1
% use topoplot to obtain grid (pixels) per subject 
for iSub = 1:size(subject,2)
    figure(100); subplot(121); title('combined cap'); [~,grid_comb(:,:,iSub),~,~,~] = topoplot(allVEPstd.comb(:,combN1.loc,iSub), chanlocs);
    figure(100); subplot(122); title('eeg only cap'); [~,grid_eeg_only(:,:,iSub),~,~,~] = topoplot(allVEPstd.eeg_only(:,eeg_onlyN1.loc,iSub), chanlocs);
    close Figure 100
end

% figure, title('N1 peak - Standardized Data')
% permutation testing
[permN1std.zmap, permN1std.zmapTH, permN1std.zmapTHmcc] = Permutation_VEPcomponent (grid_comb, grid_eeg_only, xmesh, ymesh);
% percentage of pixels significant
nrpix_total         = numel(grid_comb(:,:,1)) - sum(isnan(grid_comb(:,:,1)),'all');
nrpix_sig           = sum(permN1std.zmapTHmcc==1,'all');
permN1std.sigPix    = 100*nrpix_sig / nrpix_total; 

clear grid_comb grid_eeg_only nrpix_sig nrix_total 

% eeg_only-P1 VS comb-P1
% use topoplot to obtain grid (pixels) per subject 
for iSub = 1:size(subject,2)
    figure(100); subplot(121); title('combined cap'); [~,grid_comb(:,:,iSub),~,~,~] = topoplot(allVEPstd.comb(:,combP1.loc,iSub), chanlocs);
    figure(100); subplot(122); title('eeg only cap'); [~,grid_eeg_only(:,:,iSub),~,~,~] = topoplot(allVEPstd.eeg_only(:,eeg_onlyP1.loc,iSub), chanlocs);
    close Figure 100
end



% figure, title('P1 peak - Standardized Data')
% permutation testing
[permP1std.zmap, permP1std.zmapTH, permP1std.zmapTHmcc] = Permutation_VEPcomponent (grid_comb, grid_eeg_only, xmesh, ymesh);
% percentage of pixels significant
nrpix_total         = numel(grid_comb(:,:,1)) - sum(isnan(grid_comb(:,:,1)),'all');
nrpix_sig           = sum(permP1std.zmapTHmcc==1,'all');
permP1std.sigPix    = 100*nrpix_sig / nrpix_total; 

clear grid_comb grid_eeg_only nrpix_sig nrix_total

% plot topoplot with significant regions
PermTopoplot (permN1std, permP1std, eeg_onlyN1std, combN1std, eeg_onlyP1std, combP1std, chanlocs, xmesh, ymesh, 16, 'std');

saveas(gcf, 'Fig_CHECK_erp/Topoplots_standard_data.jpg');
%% permutation testing: normalized data

% eeg_only-N1 VS comb-N1
% use topoplot to obtain grid (pixels) per subject 
% for iSub = 1:size(subject,2)
%     figure(100); subplot(121); title('combEEG'); [~,grid_comb(:,:,iSub),~,~,~] = topoplot(allVEPnrm.comb(:,combN1.loc,iSub), chanlocs);
%     figure(100); subplot(122); title('eeg_onlyEEG'); [~,grid_eeg_only(:,:,iSub),~,~,~] = topoplot(allVEPnrm.eeg_only(:,eeg_onlyN1.loc,iSub), chanlocs);
%     close Figure 100
% end
% % permutation testing
% [permN1nrm.zmap, permN1nrm.zmapTH, permN1nrm.zmapTHmcc] = Permutation_VEPcomponent (grid_comb, grid_eeg_only, xmesh, ymesh);
% % percentage of pixels significant
% nrpix_total         = numel(grid_comb(:,:,1)) - sum(isnan(grid_comb(:,:,1)),'all');
% nrpix_sig           = sum(permN1nrm.zmapTHmcc==1,'all');
% permN1nrm.sigPix    = 100*nrpix_sig / nrpix_total; 
% 
% clear grid_comb grid_eeg_only nrpix_sig nrix_total
% 
% % eeg_only-P1 VS comb-P1
% % use topoplot to obtain grid (pixels) per subject 
% for iSub = 1:size(subject,2)
%     figure(100); subplot(121); title('combEEG'); [~,grid_comb(:,:,iSub),~,~,~] = topoplot(allVEPnrm.comb(:,combP1.loc,iSub), chanlocs);
%     figure(100); subplot(122); title('eeg_onlyEEG'); [~,grid_eeg_only(:,:,iSub),~,~,~] = topoplot(allVEPnrm.eeg_only(:,eeg_onlyP1.loc,iSub), chanlocs);
%     close Figure 100
% end
% % permutation testing
% [permP1nrm.zmap, permP1nrm.zmapTH, permP1nrm.zmapTHmcc] = Permutation_VEPcomponent (grid_comb, grid_eeg_only, xmesh, ymesh);
% % percentage of pixels significant
% nrpix_total         = numel(grid_comb(:,:,1)) - sum(isnan(grid_comb(:,:,1)),'all');
% nrpix_sig           = sum(permP1nrm.zmapTHmcc==1,'all');
% permP1nrm.sigPix    = 100*nrpix_sig / nrpix_total; 
% 
% clear grid_comb grid_eeg_only nrpix_sig nrix_total
% 
% % plot topoplot with significant regions
% PermTopoplot (permN1nrm, permP1nrm, eeg_onlyN1nrm, combN1nrm, eeg_onlyP1nrm, combP1nrm, chanlocs, xmesh, ymesh, 17, 'nrm');
% saveas(gcf, 'Fig_CHECK_erp/Topoplots_normalized_data.jpg');


fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
fprintf('Statistics: \n')
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n \n')
fprintf('RMSD: %d \n', RMSD)
fprintf('Spearman´s rank correlation: %d \n', COR)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END SCRIPT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END SCRIPT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [Ploc, Ptime, Pamp] = get_mainVEPcomponents (time, avgGFP, window)

    % find main VEP component
    Pwin = find(time>=window(1,1) & time<=window(1,2));     % time window in which peak should be found
    [Pamp,idx_loc] = max(avgGFP(Pwin));                     % find max peak in average GFPvalues within time window
    Ploc = Pwin(idx_loc); clear idx_loc;                    % sample number of peak
    Ptime = time(Ploc);                                     % time (s) of  peak

end

%%
function [zmap, zmapTH, zmapTHmcc] = Permutation_VEPcomponent (grid_comb, grid_eeg_only, xmesh, ymesh)

    % settings
    n = size(grid_comb,3);
    n_permutes = 2^n; 
    vox_pval = 0.05; 
    mcc_vox_pval = 0.05;
    
    % compute real paired t-test of difference (observed t-value)
    Xdiff   = grid_comb - grid_eeg_only;
    treal   = mean(Xdiff,3) ./ (std(Xdiff,[],3) ./ sqrt(n));
    
    % grid of both conditions (comb/eeg_only)
    grid_all = cat(3,grid_comb,grid_eeg_only);
    real_cond_map = [zeros(1,n) , ones(1,n)];

    % create mapping of all possible permutations with n = 15
    % > make (2^n x n) matrix with all possible combinations of 1s and 0s
    mat = nan(2^n,n);
    for idx = 1:n
        mat(:,idx) = repmat([zeros(2^(n-idx),1);ones(2^(n-idx),1)],2^(idx-1),1);
    end
    perm_cond_map = [mat 1-mat];
    clear mat idx 

    % initialize null hypothesis matrices
    permuted_tvals  = zeros(n_permutes,size(grid_all,1), size(grid_all,2));
    max_pixel_pvals = zeros(n_permutes,2);
    
    % generate pixel-specific null hypothesis parameter distribution
    for iperm = 1:n_permutes
        % paired t-test statistic
        Xdiff   = grid_all(:,:,perm_cond_map(iperm,:)==0) - grid_all(:,:,perm_cond_map(iperm,:)==1);
        tmap    = mean(Xdiff,3) ./ (std(Xdiff,[],3) ./ sqrt(n));
        % save permuted values
        permuted_tvals(iperm,:,:) = tmap;
        % save maximum pixel values
        max_pixel_pvals(iperm,:) = [min(tmap(:)) max(tmap(:))];
    end

    zmap = (treal - squeeze(mean(permuted_tvals,1))) ./ squeeze(std(permuted_tvals,[],1));

    figure;
    subplot(131);
    contourf(xmesh,ymesh,zmap,40,'linecolor','none');
    axis square; title('unthresholded zmap');

    % apply uncorrected threshold
    subplot(132);
    contourf(xmesh,ymesh, zmap,40,'linecolor', 'none');
    zmapTH = zmap;
    % subthreshod = false; suprathreshold = true;
    zmapTH(abs(zmapTH)< norminv(1-vox_pval/2))=false; % 2-tailed, so norminv(0.05/2) = 1.96
    zmapTH(abs(zmapTH)>0 & ~isnan(zmapTH)) = true;
    hold on;
    contour(xmesh,ymesh,zmapTH,1,'linecolor','k');
    axis square; title('unthresholded zmap')

    % apply pixel-level corrected threshold
    low_TH = prctile(max_pixel_pvals(:,1), mcc_vox_pval*100/2);
    upp_TH = prctile(max_pixel_pvals(:,2), 100-mcc_vox_pval*100/2);
    zmapTHmcc = zmap;
    zmapTHmcc(zmapTHmcc>low_TH & zmapTHmcc<upp_TH) = false;
    zmapTHmcc(abs(zmapTHmcc)>0 & ~isnan(zmapTHmcc)) = true;
    subplot(133)
    contourf(xmesh,ymesh, zmap,40,'linecolor', 'none'); hold on;
    contour(xmesh,ymesh, zmapTHmcc,1,'linecolor', 'k');
    axis square; title('pixel-corrected zmap')

end

%%
function PermTopoplot (permN1, permP1, eeg_onlyN1, combN1, eeg_onlyP1, combP1, chanlocs, xmesh, ymesh, fignr, type)

    figure(fignr); 
    subplot(221); title('N1: eeg only cap'); hold on;
    topoplot(eeg_onlyN1.topo, chanlocs, 'electrodes','off','numcontour',0);
    ax(1) = gca; axis on; hold on;
    contour(xmesh,ymesh,permN1.zmapTHmcc,1,'k','LineWidth',1); 
    axis off; c1 = colorbar; caxlim(1,:) = caxis; 
    subplot(222); title('N1: combined cap'); hold on;
    topoplot(combN1.topo, chanlocs, 'electrodes','off','numcontour',0);
    ax(2) = gca; axis on; hold on;
    contour(xmesh,ymesh,permN1.zmapTHmcc,1,'k','LineWidth',1); 
    axis off; c2 = colorbar; caxlim(2,:) = caxis; 

    set(ax,'clim',[-max(caxlim(:,2)) max(caxlim(:,2))]); 
    c1.Ticks = [-floor(max(caxlim(:,2))) 0 floor(max(caxlim(:,2)))]; 
    c2.Ticks = [-floor(max(caxlim(:,2))) 0 floor(max(caxlim(:,2)))]; 
    clear caxlim
    
    subplot(223); title('P1: eeg only cap'); hold on;
    topoplot(eeg_onlyP1.topo, chanlocs, 'electrodes','off','numcontour',0);
    ax(1) = gca; axis on; hold on;
    contour(xmesh,ymesh,permP1.zmapTHmcc,1,'k','LineWidth',1); 
    axis off; c3 = colorbar; caxlim(1,:) = caxis; 
    subplot(224); title('P1: combined cap'); hold on;
    topoplot(combP1.topo, chanlocs, 'electrodes','off','numcontour',0);
    ax(2) = gca; axis on; hold on;
    contour(xmesh,ymesh,permP1.zmapTHmcc,1,'k','LineWidth',1); 
    axis off; c4 = colorbar; caxlim(2,:) = caxis; 

    set(ax,'clim',[-max(caxlim(:,2)) max(caxlim(:,2))]); 
    c3.Ticks = [-floor(max(caxlim(:,2))) 0 floor(max(caxlim(:,2)))]; 
    c4.Ticks = [-floor(max(caxlim(:,2))) 0 floor(max(caxlim(:,2)))]; 
    
    if strcmp(type,'raw')==1
        ylabel(c1,'Potential (\muV)','FontSize',9); 
        ylabel(c2,'Potential (\muV)','FontSize',9);
        ylabel(c3,'Potential (\muV)','FontSize',9); 
        ylabel(c4,'Potential (\muV)','FontSize',9);
    elseif strcmp(type,'std')==1
        ylabel(c1,'Z-score','FontSize',11); 
        ylabel(c2,'Z-score','FontSize',11); 
        ylabel(c3,'Z-score','FontSize',11); 
        ylabel(c4,'Z-score','FontSize',11);
    elseif strcmp(type,'nrm')==1
        ylabel(c1,'min-max','FontSize',9); 
        ylabel(c2,'min-max','FontSize',9); 
        ylabel(c3,'min-max','FontSize',9); 
        ylabel(c4,'min-max','FontSize',9); 
    end

end