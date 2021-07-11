R%% VAL_Analysis_comb

close all; clear all; clc;

%% Settings
% add path to correct folders
% indicate ID number and cap 
% create filenames to save data in later stadium

laptop = 'laptopJanne';
addFolders;

mainpath_out    = 'C:\1. PROMPT\A. Dry-EEG\4. Analysis\AnalysisOutput\';
mainpath_outfig = 'C:\1. PROMPT\A. Dry-EEG\4. Analysis\AnalysisFigures\';

%% General variables

subject     = [{'03'};{'04'};{'10'};{'11'};{'12'};{'13'};{'14'};{'15'};{'16'};{'17'};{'18'};{'19'};{'20'};{'21'};{'22'}];
chans       = {'Fp1';'Fpz';'Fp2';'F7';'F3';'Fz';'F4';'F8';'FC5';'FC1';'FC2';'FC6';'T7';'C3';'Cz';'C4';'T8';'CP5';'CP1';'CP2';'CP6';'P7';'P3';'Pz';'P4';'P8';'O1';'Oz';'O2'};
time        = -0.1 : 1/1024 : 0.4;

fig_erp     = fullfile([mainpath_outfig,'ALL_VEP.jpg']);
fig_topo    = fullfile([mainpath_outfig,'ALL_VEPtopo.jpg']);

%% CHECK - ERP

for iSub = 1:size(subject,1)
    %% Load and collect subject data
    
    % load data
    ID = char(subject(iSub));
    load([mainpath_out,'ID',ID,'\ResultsMatrix\ID',ID,'_dryVEP.mat']); dry = CHECKerp; clear CHECKerp;
    load([mainpath_out,'ID',ID,'\ResultsMatrix\ID',ID,'_gelVEP.mat']); gel = CHECKerp; clear CHECKerp;
    
    % collect VEP of all subjects [chan x time x subject]
    for iChan = 1:length(chans)
        gelloc = find(strcmp(chans(iChan), gel.chans));
        if ~isempty(gelloc);        allVEP.gel(iChan,:,iSub) = gel.ERP(gelloc,:);
        elseif isempty(gelloc);     allVEP.gel(iChan,:,iSub) = NaN;
        end
        dryloc = find(strcmp(chans(iChan), dry.chans));
        if ~isempty(dryloc);        allVEP.dry(iChan,:,iSub) = dry.ERP(dryloc,:);
        elseif isempty(dryloc);     allVEP.dry(iChan,:,iSub) = NaN;
        end
        clear gelloc dryloc
    end
    
    % collect GFPt of all subjects [subject x time]
    allGFP.gel(iSub,:) = gel.GFPt;
    allGFP.dry(iSub,:) = dry.GFPt;
    
    %% Root Mean Square Deviation (RMSD) - Fiedler (2015)
    for iL = 1:length(dry.GFPt)
        D(iL) = (gel.GFPt(iL)-dry.GFPt(iL)).^2;
    end
    RMSD(iSub,1) = sqrt(sum(D)/length(dry.GFPt));    
    
    %% Spearman's rank Correlation (COR) - Fiedler (2015)
    for iL = 1:length(dry.GFPt)
        Dg(iL) = gel.GFPt(iL)-mean(gel.GFPt);
        Dd(iL) = dry.GFPt(iL)-mean(dry.GFPt);
    end
    COR(iSub,1) = sum(Dg.*Dd) / sqrt( sum(Dg.^2).*sum(Dd.^2) ); 
    
    clear D Dg Dd dry gel
end

%% test for normality of COR and RMSD

[~,NormalDistr.COR]   = kstest(COR);
[~,NormalDistr.RMSD]  = kstest(RMSD);

%% average RMSD and COR

RES.RMSD(1,1)   = mean(RMSD);
RES.RMSD(1,2)   = 2*std(RMSD);
RES.COR(1,1)    = mean(COR);
RES.COR(1,2)    = 2*std(COR);

%% plot mean per channel over subjects + mean and std of GFPt

gelGFP1 = mean(allGFP.gel,1) + 2*std(allGFP.gel,[],1);
gelGFP2 = mean(allGFP.gel,1) - 2*std(allGFP.gel,[],1);
dryGFP1 = mean(allGFP.dry,1) + 2*std(allGFP.dry,[],1);
dryGFP2 = mean(allGFP.dry,1) - 2*std(allGFP.dry,[],1);
figure; 
subplot(2,2,1); plot(time,mean(allVEP.gel,3,'omitnan'),'b');
subplot(2,2,2); plot(time,mean(allVEP.dry,3,'omitnan'),'r');
subplot(2,2,3); hold on; h1 = fill([time,fliplr(time)], [gelGFP1,fliplr(gelGFP2)],'b','LineStyle','none');
set(h1,'FaceAlpha',0.4); plot(time,mean(allGFP.gel,1),'b','LineWidth',1.5); 
subplot(2,2,4); hold on; h2 = fill([time,fliplr(time)], [dryGFP1,fliplr(dryGFP2)],'r','LineStyle','none');
set(h2,'FaceAlpha',0.4); plot(time,mean(allGFP.dry,1),'r','LineWidth',1.5); 

subplot(2,2,1); set(gca,'FontSize',11); box on;
ylabel('Potential (\muV)','FontSize',14); title('gelEEG','FontSize',14); ylim([-6 6])
xticks([-0.1:0.1:0.4]); yticks([-5:5:5]); line([0 0],[-6 6],'Color','k'); 
subplot(2,2,2); set(gca,'FontSize',11); box on;
ylabel('Potential (\muV)','FontSize',14); title('dryEEG','FontSize',14); ylim([-6 6])
xticks([-0.1:0.1:0.4]); yticks([-5:5:5]); line([0 0],[-6 6],'Color','k'); 
subplot(2,2,3); set(gca,'FontSize',11); box on;
ylabel('GFPt (\muV)','FontSize',14); xlabel('Time (s)','FontSize',14); ylim([-2 26])
xticks([-0.1:0.1:0.4]); yticks([0:10:20]); line([0 0],[-2 25],'Color','k'); 
subplot(2,2,4); set(gca,'FontSize',11); box on;
ylabel('GFPt (\muV)','FontSize',14); xlabel('Time (s)','FontSize',14); ylim([-2 26])
xticks([-0.1:0.1:0.4]); yticks([0:10:20]); line([0 0],[-2 25],'Color','k'); 
set(gcf,'Position',[400 420 1325 420]);

clear gelGFP1 gelGFP2 dryGFP1 dryGFP2 h1 h2

% saveas(gcf,fig_erp);

%% find main VEP components 
% find the main VEP components in the average GFPt 

% average GFPt and VEP over subjects
avgGFP.gel = mean(allGFP.gel,1);
avgGFP.dry = mean(allGFP.dry,1);
avgVEP.gel = mean(allVEP.gel,3, 'omitnan');
avgVEP.dry = mean(allVEP.dry,3, 'omitnan');

% get main component N1: 0.11-0.13s (N75 peak)
[gelN1.loc, gelN1.time, gelN1.amp] = get_mainVEPcomponents (time, avgGFP.gel, [0.11 0.13]);
[dryN1.loc, dryN1.time, dryN1.amp] = get_mainVEPcomponents (time, avgGFP.dry, [0.11 0.13]);

% get main component P1: 0.14-0.16s (P100 peak)
[gelP1.loc, gelP1.time, gelP1.amp] = get_mainVEPcomponents (time, avgGFP.gel, [0.14 0.16]);
[dryP1.loc, dryP1.time, dryP1.amp] = get_mainVEPcomponents (time, avgGFP.dry, [0.14 0.16]);

%% standardize VEP data within subject per cap

% standardize (x-meanX)/sdX
allVEPstd.gel = (allVEP.gel - mean(allVEP.gel,1,'omitnan')) ./ std(allVEP.gel,[],1,'omitnan');
allVEPstd.dry = (allVEP.dry - mean(allVEP.dry,1,'omitnan')) ./ std(allVEP.dry,[],1,'omitnan');
% > average standardized VEP
avgVEPstd.gel = mean(allVEPstd.gel,3, 'omitnan');
avgVEPstd.dry = mean(allVEPstd.dry,3, 'omitnan');

% normalize data 
allVEPnrm.gel = (allVEP.gel - min(allVEP.gel,1,'omitnan')) ./ (max(allVEP.gel,1,'omitnan') - min(allVEP.gel,1,'omitnan'));
allVEPnrm.dry = (allVEP.dry - min(allVEP.dry,1,'omitnan')) ./ (max(allVEP.dry,1,'omitnan') - min(allVEP.dry,1,'omitnan'));
% > average normalized VEP
avgVEPnrm.gel = mean(allVEPnrm.gel,3, 'omitnan');
avgVEPnrm.dry = mean(allVEPnrm.dry,3, 'omitnan');

%% data for topoplot of main VEP components, average over subjects

% raw: get topoplot data of N1 and P1 peak
gelN1.topo = avgVEP.gel(:,gelN1.loc);        gelP1.topo = avgVEP.gel(:,gelP1.loc);
dryN1.topo = avgVEP.dry(:,dryN1.loc);        dryP1.topo = avgVEP.dry(:,dryP1.loc);

% std: get topoplot data of N1 and P1 standardized within subject and cap
gelN1std.topo = avgVEPstd.gel(:,gelN1.loc);    gelP1std.topo = avgVEPstd.gel(:,gelP1.loc);
dryN1std.topo = avgVEPstd.dry(:,dryN1.loc);    dryP1std.topo = avgVEPstd.dry(:,dryP1.loc);

% norm: get topoplot data of N1 and P1 normalized within subject and cap
gelN1nrm.topo = avgVEPnrm.gel(:,gelN1.loc);    gelP1nrm.topo = avgVEPnrm.gel(:,gelP1.loc);
dryN1nrm.topo = avgVEPnrm.dry(:,dryN1.loc);    dryP1nrm.topo = avgVEPnrm.dry(:,dryP1.loc);

%% permutation testing: raw data

load('chanlocs.mat')

% GEL-N1 VS DRY-N1
% use topoplot to obtain grid (pixels) per subject 
for iSub = 1:size(subject,1)
    figure(100); subplot(121); title('dryEEG'); [~,grid_dry(:,:,iSub),~,xmesh,ymesh] = topoplot(allVEP.dry(:,dryN1.loc,iSub), chanlocs);
    figure(100); subplot(122); title('gelEEG'); [~,grid_gel(:,:,iSub),~,~,~] = topoplot(allVEP.gel(:,gelN1.loc,iSub), chanlocs);
    close Figure 100
end
% permutation testing
[permN1.zmap, permN1.zmapTH, permN1.zmapTHmcc] = Permutation_VEPcomponent (grid_dry, grid_gel, xmesh, ymesh);

% percentage of significant pixels 
nrpix_total     = numel(grid_dry(:,:,1)) - sum(isnan(grid_dry(:,:,1)),'all');
nrpix_sig       = sum(permN1.zmapTHmcc==1,'all');
permN1.sigPix   = 100*nrpix_sig / nrpix_total; 

clear grid_dry grid_gel nrpix_sig nrix_total

% GEL-P1 VS DRY-P1
% use topoplot to obtain grid (pixels) per subject 
for iSub = 1:size(subject,1)
    figure(100); subplot(121); title('dryEEG'); [~,grid_dry(:,:,iSub),~,~,~] = topoplot(allVEP.dry(:,dryP1.loc,iSub), chanlocs);
    figure(100); subplot(122); title('gelEEG'); [~,grid_gel(:,:,iSub),~,~,~] = topoplot(allVEP.gel(:,gelP1.loc,iSub), chanlocs);
    close Figure 100
end
% permutation testing
[permP1.zmap, permP1.zmapTH, permP1.zmapTHmcc] = Permutation_VEPcomponent (grid_dry, grid_gel, xmesh, ymesh);

% percentage of significant pixels 
nrpix_total     = numel(grid_dry(:,:,1)) - sum(isnan(grid_dry(:,:,1)),'all');
nrpix_sig       = sum(permP1.zmapTHmcc==1,'all');
permP1.sigPix   = 100*nrpix_sig / nrpix_total;

clear grid_dry grid_gel nrpix_sig nrix_total

% plot topoplot with significant regions
PermTopoplot (permN1, permP1, gelN1, dryN1, gelP1, dryP1, chanlocs, xmesh, ymesh, 15, 'raw');

%% permutation testing: standardized data

% GEL-N1 VS DRY-N1
% use topoplot to obtain grid (pixels) per subject 
for iSub = 1:size(subject,1)
    figure(100); subplot(121); title('dryEEG'); [~,grid_dry(:,:,iSub),~,~,~] = topoplot(allVEPstd.dry(:,dryN1.loc,iSub), chanlocs);
    figure(100); subplot(122); title('gelEEG'); [~,grid_gel(:,:,iSub),~,~,~] = topoplot(allVEPstd.gel(:,gelN1.loc,iSub), chanlocs);
    close Figure 100
end
% permutation testing
[permN1std.zmap, permN1std.zmapTH, permN1std.zmapTHmcc] = Permutation_VEPcomponent (grid_dry, grid_gel, xmesh, ymesh);
% percentage of pixels significant
nrpix_total         = numel(grid_dry(:,:,1)) - sum(isnan(grid_dry(:,:,1)),'all');
nrpix_sig           = sum(permN1std.zmapTHmcc==1,'all');
permN1std.sigPix    = 100*nrpix_sig / nrpix_total; 

clear grid_dry grid_gel nrpix_sig nrix_total 

% GEL-P1 VS DRY-P1
% use topoplot to obtain grid (pixels) per subject 
for iSub = 1:size(subject,1)
    figure(100); subplot(121); title('dryEEG'); [~,grid_dry(:,:,iSub),~,~,~] = topoplot(allVEPstd.dry(:,dryP1.loc,iSub), chanlocs);
    figure(100); subplot(122); title('gelEEG'); [~,grid_gel(:,:,iSub),~,~,~] = topoplot(allVEPstd.gel(:,gelP1.loc,iSub), chanlocs);
    close Figure 100
end
% permutation testing
[permP1std.zmap, permP1std.zmapTH, permP1std.zmapTHmcc] = Permutation_VEPcomponent (grid_dry, grid_gel, xmesh, ymesh);
% percentage of pixels significant
nrpix_total         = numel(grid_dry(:,:,1)) - sum(isnan(grid_dry(:,:,1)),'all');
nrpix_sig           = sum(permP1std.zmapTHmcc==1,'all');
permP1std.sigPix    = 100*nrpix_sig / nrpix_total; 

clear grid_dry grid_gel nrpix_sig nrix_total

% plot topoplot with significant regions
PermTopoplot (permN1std, permP1std, gelN1std, dryN1std, gelP1std, dryP1std, chanlocs, xmesh, ymesh, 16, 'std');

% %% permutation testing: normalized data
% 
% % GEL-N1 VS DRY-N1
% % use topoplot to obtain grid (pixels) per subject 
% for iSub = 1:size(subject,1)
%     figure(100); subplot(121); title('dryEEG'); [~,grid_dry(:,:,iSub),~,~,~] = topoplot(allVEPnrm.dry(:,dryN1.loc,iSub), chanlocs);
%     figure(100); subplot(122); title('gelEEG'); [~,grid_gel(:,:,iSub),~,~,~] = topoplot(allVEPnrm.gel(:,gelN1.loc,iSub), chanlocs);
%     close Figure 100
% end
% % permutation testing
% [permN1nrm.zmap, permN1nrm.zmapTH, permN1nrm.zmapTHmcc] = Permutation_VEPcomponent (grid_dry, grid_gel, xmesh, ymesh);
% % percentage of pixels significant
% nrpix_total         = numel(grid_dry(:,:,1)) - sum(isnan(grid_dry(:,:,1)),'all');
% nrpix_sig           = sum(permN1nrm.zmapTHmcc==1,'all');
% permN1nrm.sigPix    = 100*nrpix_sig / nrpix_total; 
% 
% clear grid_dry grid_gel nrpix_sig nrix_total
% 
% % GEL-P1 VS DRY-P1
% % use topoplot to obtain grid (pixels) per subject 
% for iSub = 1:size(subject,1)
%     figure(100); subplot(121); title('dryEEG'); [~,grid_dry(:,:,iSub),~,~,~] = topoplot(allVEPnrm.dry(:,dryP1.loc,iSub), chanlocs);
%     figure(100); subplot(122); title('gelEEG'); [~,grid_gel(:,:,iSub),~,~,~] = topoplot(allVEPnrm.gel(:,gelP1.loc,iSub), chanlocs);
%     close Figure 100
% end
% % permutation testing
% [permP1nrm.zmap, permP1nrm.zmapTH, permP1nrm.zmapTHmcc] = Permutation_VEPcomponent (grid_dry, grid_gel, xmesh, ymesh);
% % percentage of pixels significant
% nrpix_total         = numel(grid_dry(:,:,1)) - sum(isnan(grid_dry(:,:,1)),'all');
% nrpix_sig           = sum(permP1nrm.zmapTHmcc==1,'all');
% permP1nrm.sigPix    = 100*nrpix_sig / nrpix_total; 
% 
% clear grid_dry grid_gel nrpix_sig nrix_total
% 
% % plot topoplot with significant regions
% PermTopoplot (permN1nrm, permP1nrm, gelN1nrm, dryN1nrm, gelP1nrm, dryP1nrm, chanlocs, xmesh, ymesh, 17, 'nrm');

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
function [zmap, zmapTH, zmapTHmcc] = Permutation_VEPcomponent (grid_dry, grid_gel, xmesh, ymesh)

    % settings
    n = size(grid_dry,3);
    n_permutes = 2^n; 
    vox_pval = 0.05; 
    mcc_vox_pval = 0.05;
    
    % compute real paired t-test of difference (observed t-value)
    Xdiff   = grid_dry - grid_gel;
    treal   = mean(Xdiff,3) ./ (std(Xdiff,[],3) ./ sqrt(n));
    
    % grid of both conditions (dry/gel)
    grid_all = cat(3,grid_dry,grid_gel);
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
function PermTopoplot (permN1, permP1, gelN1, dryN1, gelP1, dryP1, chanlocs, xmesh, ymesh, fignr, type)

    figure(fignr); 
    subplot(221); title('N1: gel-EEG'); hold on;
    topoplot(gelN1.topo, chanlocs, 'electrodes','off','numcontour',0);
    ax(1) = gca; axis on; hold on;
    contour(xmesh,ymesh,permN1.zmapTHmcc,1,'k','LineWidth',1); 
    axis off; c1 = colorbar; caxlim(1,:) = caxis; 
    subplot(222); title('N1: dry-EEG'); hold on;
    topoplot(dryN1.topo, chanlocs, 'electrodes','off','numcontour',0);
    ax(2) = gca; axis on; hold on;
    contour(xmesh,ymesh,permN1.zmapTHmcc,1,'k','LineWidth',1); 
    axis off; c2 = colorbar; caxlim(2,:) = caxis; 

    set(ax,'clim',[-max(caxlim(:,2)) max(caxlim(:,2))]); 
    c1.Ticks = [-floor(max(caxlim(:,2))) 0 floor(max(caxlim(:,2)))]; 
    c2.Ticks = [-floor(max(caxlim(:,2))) 0 floor(max(caxlim(:,2)))]; 
    clear caxlim
    
    subplot(223); title('P1: gel-EEG'); hold on;
    topoplot(gelP1.topo, chanlocs, 'electrodes','off','numcontour',0);
    ax(1) = gca; axis on; hold on;
    contour(xmesh,ymesh,permP1.zmapTHmcc,1,'k','LineWidth',1); 
    axis off; c3 = colorbar; caxlim(1,:) = caxis; 
    subplot(224); title('P1: dry-EEG'); hold on;
    topoplot(dryP1.topo, chanlocs, 'electrodes','off','numcontour',0);
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