clear;
close all;

%% Specify where 3D model folder is saved
laptop='laptopCatarina';
% laptop='laptopJoao';
% laptop='laptopMariana';
[mainpath_in, mainpath_out] = addFolders(laptop);
ft_defaults;

sub='28';

model_path = fullfile(mainpath_in,'source',['sub-',sub],'3d');
obj_file = fullfile(model_path, 'Model', 'Model.obj');
out_path = fullfile(mainpath_out,['sub-',sub],'3d');
cd(out_path) % folder to save the data

[ftver, ftpath] = ft_version; % specify your path to fieldtrip

%% Load the 3D model
head_surface = ft_read_headshape(obj_file);

% Convert the units to mm
head_surface = ft_convert_units(head_surface, 'mm');

% Visualize the mesh surface
ft_plot_mesh(head_surface)

%% Specify the fiducials + optodes
% note that we call Cz a fiducial (which is not a fiducial) but we will use
% it in a later stage for coregistration
if isfile('opto.mat')
     load(['opto.mat'], 'opto')
else
    cfg = [];
    cfg.channel={};
    cfg.channel{1} = 'Nz';
    cfg.channel{2} = 'LPA';
    cfg.channel{3} = 'RPA';
    cfg.channel{4}= 'Iz';
    cfg.channel{5} = 'Cz';
    tnames = {'Tx1a', 'Tx1b', 'Tx1c', 'Tx1d', 'Tx6a', 'Tx6b', 'Tx6c', 'Tx6d', 'Tx11a', 'Tx11b', 'Tx11c', 'Tx11d', 'Tx2', 'Tx3', 'Tx4', 'Tx5', 'Tx7', 'Tx8', 'Tx9', 'Tx10', 'Tx12', 'Tx13', 'Tx14', 'Tx15'};
    for i =6:29
        cfg.channel{i}=sprintf('%s', tnames{i-5}); % change this according to the number of transmitters your layout has
    end
    rnames = {'Rx1', 'Rx2', 'Rx3', 'Rx4', 'Rx5', 'Rx6', 'Rx7', 'Rx8', 'Rx9', 'Rx10', 'Rx11', 'Rx12'};
    for i =30:41
       cfg.channel{i}=sprintf('%s', rnames{i-29}); % change this according to the number of receivers your layout has
    end
    cfg.method = 'headshape';
    opto = ft_electrodeplacement(cfg, head_surface);
    save(['opto.mat'], 'opto')
end


% Notes:
% Do you want to change the anatomical labels for the axes [Y, n]? --> n
% Use "Tools > Rotate 3D" to rotate the 3D model
% Click/unclick "Colors" to toggle the colors on and off: best is to use 
% the color view for the fiducials, but the structure view for the optodes
% Use the mouse to click on fiducials/optodes and subsequently on the 
% corresponding label to assign the markers (make sure you're not in 
% "Rotate 3D" mode anymore!)
% If an error was made: double click on the label to remove this marker
% If ready --> press Q

%% Align the axes of the coordinate system with the fiducial positions (ctf coordinates)
% For the mesh
clf;
cfg = [];
cfg.method = 'fiducial';
cfg.coordsys = 'ctf';
cfg.fiducial.nas = opto.elecpos(1,:); % position of Nz
cfg.fiducial.lpa = opto.elecpos(2,:); % position of LPA
cfg.fiducial.rpa = opto.elecpos(3,:); % position of RPA
head_surface_aligned = ft_meshrealign(cfg, head_surface);

ft_plot_axes(head_surface_aligned)
ft_plot_mesh(head_surface_aligned)

% For the optodes
fid.chanpos = [110 0 0; 0 90 0; 0 -90 0]; % CTF coordinates of the fiducials
fid.elecpos = [110 0 0; 0 90 0; 0 -90 0]; % like electrode positions
fid.label = {'Nz','LPA','RPA'}; % same labels as in elec
fid.unit = 'mm'; % same units as mri

cfg = [];
cfg.method = 'fiducial';
cfg.coordsys = 'ctf';
cfg.target = fid;
cfg.elec = opto;
cfg.fiducial = {'Nz', 'LPA', 'RPA'}; % labels of fiducials in fid and in elec
opto_aligned = ft_electroderealign(cfg);

% Visualize the optodes on the aligned head surface
% the colorview is not that clean as the structure view
figure;
ft_plot_mesh(head_surface_aligned)
ft_plot_sens(opto_aligned, 'elecsize', 10, 'style', 'b')

% Same visualization without the colors
figure;
ft_plot_mesh(removefields(head_surface_aligned, 'color'), 'tag', 'headshape', 'facecolor', 'skin', 'material', 'dull', 'edgecolor', 'none', 'facealpha', 1);
lighting gouraud
l = lightangle(0, 90); set(l, 'Color', [1 1 1]/2)
l = lightangle(0, 0); set(l, 'Color', [1 1 1]/3)
l = lightangle(90, 0); set(l, 'Color', [1 1 1]/3)
l = lightangle(180, 0); set(l, 'Color', [1 1 1]/3)
l = lightangle(270, 0); set(l, 'Color', [1 1 1]/3)
alpha 0.9
ft_plot_sens(opto_aligned, 'elecsize', 10, 'style', 'b')

save(['opto_aligned.mat'], 'opto_aligned')

%% Move optode inward
cfg = [];
cfg.method = 'moveinward';
cfg.moveinward = 5; % determine distance to skin
cfg.channel = 2:length(opto_aligned.label); % do not move the nasion inward!
cfg.keepchannel = true;
cfg.elec = opto_aligned;
opto_inw = ft_electroderealign(cfg);

% Visualize the optode inward
figure;
ft_plot_mesh(removefields(head_surface_aligned, 'color'), 'tag', 'headshape', 'facecolor', 'skin', 'material', 'dull', 'edgecolor', 'none', 'facealpha', 1);
lighting gouraud
l = lightangle(0, 90); set(l, 'Color', [1 1 1]/2)
l = lightangle(0, 0); set(l, 'Color', [1 1 1]/3)
l = lightangle(90, 0); set(l, 'Color', [1 1 1]/3)
l = lightangle(180, 0); set(l, 'Color', [1 1 1]/3)
l = lightangle(270, 0); set(l, 'Color', [1 1 1]/3)
alpha 0.7
ft_plot_sens(opto_inw, 'elecsize', 10, 'style', 'b')

save(['opto_inw.mat'], 'opto_inw')

%% Coregister optodes to MNI atlas 
% Load skin surface of standard atlas
skin_template = fullfile(ftpath, 'template', 'headmodel', 'skin', 'standard_skin_14038.vol');
skin = ft_read_headshape(skin_template);

% Visualize skin surface of standard atlas
figure;
ft_plot_mesh(skin, 'edgecolor', 'none', 'facecolor', 'skin'); 
camlight

% Visualization with optodes
hold on; ft_plot_sens(opto_inw);
% the optodes (represented in the ctf coordinate system) are rotated 
% with 90 degrees compared to the head model (represented in the MNI 
% coordinate system)

% Use the MNI coordinates of the fiducials to guide the coregistration 
elec1020 = ft_read_sens(fullfile(ftpath, 'template', 'electrode', 'standard_1020.elc'));

% Note:
% Select rpa, lpa, nasion, inion and cz from elec1020

cfg=[];
cfg.method = 'moveinward'; % a hack to select the right channels
cfg.moveinward = 0;
cfg.channel = {'Nz'; 'RPA'; 'LPA'; 'Iz'; 'Cz'};
electarg = ft_electroderealign(cfg, elec1020);

% There are two ways to coregister the optodes to the MNI space:
% (1) interactively/manually
% (2) automatically based on the fiducials we specified above

% (1) interactively/manually (time consuming and needs practice. the
% automatic method is better but only when you have enough fiducials/marker
% points)
% 90 degree turn to go from ctf to MNI coordinate system
%{
cfg = [];
cfg.method = 'interactive';
cfg.headshape = skin;
cfg.target = electarg; % you can use the fiducials to guide your coregistration
opto_MNI = ft_electroderealign(cfg, opto_inw);
%}

% (2) automatically based on the fiducials we specified above 
cfg = [];
cfg.method = 'template';
cfg.warp = 'globalrescale';
cfg.channel = {'LPA', 'RPA', 'Nz', 'Cz', 'Iz'};
cfg.target = electarg;
opto_MNI = ft_electroderealign(cfg, opto_inw);

save('opto_MNI.mat', 'opto_MNI')

% Plot all together
figure;
ft_plot_mesh(skin, 'edgecolor', 'none', 'facecolor', 'skin'); camlight
hold on; 
ft_plot_sens(opto_MNI, 'elecsize', 20,'facecolor', 'k', 'label', 'label');
hold on;
ft_plot_sens(electarg, 'elecsize', 20, 'facecolor', 'b');

%% Calculate channel locations based on optode positions (in MNI space)
% define here the channels you want to create:
channels = {'Rx9-Tx12', 'Rx9-Tx13', 'Rx11-Tx12', 'Rx11-Tx13', 'Rx10-Tx14', 'Rx12-Tx14', 'Rx12-Tx15', 'Rx4-Tx4', 'Rx4-Tx5', 'Rx5-Tx7', 'Rx5-Tx8', 'Rx7-Tx7', 'Rx7-Tx8', 'Rx6-Tx9', 'Rx8-Tx9', 'Rx8-Tx10', 'Rx1-Tx2', 'Rx1-Tx3', 'Rx3-Tx2', 'Rx3-Tx3', 'Rx3-Tx5', 'Rx2-Tx4', 'Rx2-Tx3'};
%channels = {'Rx5-Tx6b','Rx7-Tx6d','Rx6-Tx6a','Rx8-Tx6c','Rx4-Tx1c','Rx2-Tx1a','Rx3-Tx1d','Rx1-Tx1b','Rx9-Tx11a','Rx11-Tx11c','Rx12-Tx11d','Rx10-Tx11b'}; %short channels
[rxnames, rem] = strtok(channels, {'-', ' '});
[txnames, rem] = strtok(rem, {'-', ' '});

% Change the naming
opto_def = opto_MNI;
opto_def.optopos = opto_MNI.chanpos;
opto_def.chantype = cell(length(opto_MNI.chantype),1);
opto_def.chantype(:) = {'nirs'};
opto_def.optolabel = opto_MNI.label;
opto_def.label = channels;
opto_def.tra = zeros(length(channels),length(opto_def.optolabel));
for i=1:length(channels)
    opto_def.tra(i,:)=strcmp(rxnames{i},opto_def.optolabel)+strcmp(txnames{i}, opto_def.optolabel);
end
opto_def = rmfield(opto_def, {'chanpos', 'chantype', 'chanunit', 'elecpos'});

% Calculate channel positions
opto_chan = ft_datatype_sens(opto_def)

% Visualize
figure;
ft_plot_mesh(skin,'edgecolor','none','facealpha',0.8,'facecolor',[0.6 0.6 0.8]);
hold on;
ft_plot_sens(opto_chan, 'opto', true, 'optosize', 10,'facecolor', 'k', 'label', 'label');

%% Determine anatomical labels for the channels
% Load atlas (AAL atlas)
atlas = ft_read_atlas([ftpath filesep 'template/atlas/aal/ROI_MNI_V4.nii']);
channels = {'Rx9-Tx12', 'Rx9-Tx13', 'Rx11-Tx12', 'Rx11-Tx13', 'Rx10-Tx14', 'Rx12-Tx14', 'Rx12-Tx15', 'Rx4-Tx4', 'Rx4-Tx5', 'Rx5-Tx7', 'Rx5-Tx8', 'Rx7-Tx7', 'Rx7-Tx8', 'Rx6-Tx9', 'Rx8-Tx9', 'Rx8-Tx10', 'Rx1-Tx2', 'Rx1-Tx3', 'Rx3-Tx2', 'Rx3-Tx3', 'Rx3-Tx5', 'Rx2-Tx4', 'Rx2-Tx3'};
%channels = {'Rx5-Tx6b','Rx7-Tx6d','Rx6-Tx6a','Rx8-Tx6c','Rx4-Tx1c','Rx2-Tx1a','Rx3-Tx1d','Rx1-Tx1b','Rx9-Tx11a','Rx11-Tx11c','Rx12-Tx11d','Rx10-Tx11b'}; %short channels

% Look up the corresponding anatomical label
cfg = [];
cfg.roi = opto_chan.chanpos(match_str(opto_chan.label,channels),:);
cfg.atlas = atlas;
cfg.output = 'multiple';
cfg.minqueryrange = 1;
cfg.maxqueryrange = 25; 

% Note:
% If no label was found, increase the queryrange

labels = ft_volumelookup(cfg, atlas);

% Select the anatomical label with the highest probability
for i=1:length(channels)
    [~, indx] = max(labels(i).count);
    label{i}=char(labels(i).name(indx));
end

% Show results in table
roi_table = table(channels', label', 'VariableNames', {'channel', 'label'});
save('roi_table.mat', 'roi_table')

%% Create layout based on 3D positions
cfg = [];
cfg.opto = opto_chan;
cfg.rotate = 0;
layout = ft_prepare_layout(cfg);

% Plot layout
cfg = [];
cfg.layout = layout;
ft_layoutplot(cfg);

% Create outline and mask based on image
image = fullfile(mainpath_in, 'plot_layout.png');
cfg = [];
cfg.image = image;
bg = ft_prepare_layout(cfg);

% Notes:
% Specify electrode locations --> only select Cz as electrode (we need
% it as the middlepoint to align the layout with the mask & outline) --> press q
% Create mask: skip this step --> press q
% Create outline: these are the lines that will be shown in the layout
% (outline the circumference with nose + central sulcus) --> press q

% Scale according to layout and place Cz in the middle
for i=1:length(bg.outline)
    bg.outline{i}=[bg.outline{i}]-bg.pos(1,:); % center around Cz
    bg.outline{i}=[bg.outline{i}]/250; % scale
end
bg.pos=bg.pos-bg.pos(1,:);
bg.pos=bg.pos/250;

% Combine all into one layout
layout.outline = bg.outline;

% Visualize
figure;
ft_plot_layout(layout);
save('layout.mat', 'layout');

% Change the mask based on the layout with outline just created
saveas(gcf, 'layout.png')
image = 'layout.png';
cfg = [];
cfg.image = image;
bg = ft_prepare_layout(cfg);

% Notes:
% Specify electrode locations --> only select Cz as electrode --> press q
% Create mask: this area will be used when making a topoplot (select optode
% areas) --> press q
% Create outline: skip this step --> press q

% Center around Cz + rescale
for i=1:length(bg.mask)
    bg.mask{i}=[bg.mask{i}]-bg.pos(1,:); % center around Cz
    bg.mask{i}=[bg.mask{i}]*(1.1/150); % scale
end

% Visualize and save the layout
load('layout.mat')
layout.mask = bg.mask;
figure; ft_plot_layout(layout)
save('layout.mat', 'layout')
saveas(gcf, 'layout.png')

