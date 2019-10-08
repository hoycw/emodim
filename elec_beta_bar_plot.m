%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'hoycw/Apps/fieldtrip/'];
elseif exist('/Users/lapate/','dir');root_dir = '/Users/lapate/knight/';ft_dir = '/Users/lapate/knight/hoycw/Apps/fieldtrip/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%% Set Up Directories
% addpath([root_dir 'emodim/scripts/']);
% addpath([root_dir 'emodim/scripts/utils/']);
% addpath(ft_dir);
% ft_defaults

%% Data Preparation
SBJ = 'IR84';
ch_lab = 'tmp';
% SBJ_vars_cmd = ['run ' root_dir 'emodim/scripts/SBJ_vars/' SBJ '_vars.m'];
% eval(SBJ_vars_cmd);

% Load data
load([root_dir 'emodim/data/predCorr_HGmLogLowpass.mat']);
load([SBJ_vars.dirs.recon SBJ '_elec_main_ft_pat_Dx_gROI.mat']);

emo_lab = {'Adoration', 'Aesthetic Appreciation', 'Amusement', 'Anger',...
    'Anxiety', 'Awe', 'Confusion', 'Disgust', 'Empathic Pain', 'Excitement',...
    'Fear', 'Horror', 'Interest', 'Joy', 'Sadness', 'Surprise'};
data = rand(size(emo_lab))-0.5;

%% Plot bars
fig_name = [SBJ '_' ch_lab '_bar_betas'];
f = figure('Name',fig_name);
ax = gca;

b = bar(numel(emo_lab):-1:1,data,'FaceColor','flat');

cmap = colormap('parula');
cmap_idx = linspace(min(data),max(data),size(cmap,1));
cdata = zeros([numel(data) 3]);
for i = 1:numel(data)
    [~, c_ix] = min(abs(data(i)-cmap_idx));
    cdata(i,:) = cmap(c_ix,:);
end
b.CData = cdata;

% Plot labels
ax.XLim       = [0 numel(data)+1];
ax.XTick      = 1:numel(data);
ax.XTickLabel = emo_lab;
ax.XLabel.String = 'Emotion Category';

ax.YLabel.String   = 'Model Weight';

ax.Title.String = 'Emotional Encoding: One Amygdala Site';

view([90 -90]);
set(ax,'FontSize',16);
set(findall(f,'type','text'),'FontSize',18);

%% Save figure
fig_fname = [root_dir 'emodim/results/' fig_name '.png'];

