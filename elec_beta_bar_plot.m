%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'hoycw/Apps/fieldtrip/'];
elseif exist('/Users/lapate/','dir');root_dir = '/Users/lapate/knight/';ft_dir = '/Users/lapate/knight/hoycw/Apps/fieldtrip/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%% Set Up Directories
addpath([root_dir 'emodim/scripts/']);
addpath([root_dir 'emodim/scripts/utils/']);
addpath([root_dir 'emodim/scripts/Colormaps/']);
addpath(ft_dir);
ft_defaults

%% Data Preparation
SBJ = 'IR84';
stat_lab = 'OrigMean';
thresh = 0.05;
SBJ_vars_cmd = ['run ' root_dir 'emodim/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

% Load data
load([root_dir 'emodim/data/predCorr_HGmLogLowpass.mat']);
load([SBJ_vars.dirs.proc SBJ '_LowPassLog_betas.mat']);
load([SBJ_vars.dirs.recon SBJ '_elec_main_ft_pat_Dx_gROI.mat']);

ch_lab = 'tmp';

emo_lab = ylabelnames;
% emo_lab = {'Adoration', 'Aesthetic Appreciation', 'Amusement', 'Anger',...
%     'Anxiety', 'Awe', 'Confusion', 'Disgust', 'Empathic Pain', 'Excitement',...
%     'Fear', 'Horror', 'Interest', 'Joy', 'Sadness', 'Surprise'};
% data = rand(size(emo_lab))-0.5;

%%
% Select sig elecs
% cors = predCorr.(SBJ).(stat_lab).Correlations;
% boot = predCorr.(SBJ).(stat_lab).Bootstrap;
% 
% % Threshold data
% sig  = false(size(elec.label));
% stat = zeros(size(elec.label));
% fprintf('%s Significant Electrodes:\n\t',stat_lab);
% for e = 1:numel(elec.label)
%     %     if strcmp(sig_method,'prop_boot')
%     % Proportion of bootstrap values above zero
%     stat(e) = 1-sum(boot(:,e)>0)/size(boot,1);
%     if stat(e) <= thresh
%         sig(e) = true;
%         fprintf('%s\t',elec.label{e});
%     end
%     %     elseif strcmp(sig_method,'norminv')
%     %         error('dont use this until sure its appropriate, see email Slides ready to edit on 12/11/18');
%     %         % STDs of bootstrap above 0
%     %         sd = std(boot{sbj_ix,st_ix}(:,e),0,1);
%     %         stat(e) = -norminv(thresh)*sd;
%     %         if cors{sbj_ix,st_ix}(e) >= stat(e)
%     %             sig(e) = true;
%     %             fprintf('%s\t',elec_sbj{sbj_ix,st_ix}.label{e});
%     %         end
%     %     end
% end
% fprintf('\nTotal n_sig = %i / %i\n',sum(sig),numel(cors));

%% Plot betas
imagesc(beta);
colormap(redblue(100));
colorbar;
ax = gca;
ax.XTick = 1:size(beta,2);
ax.XTickLabel = xlabelnames;
ax.XTickLabelRotation = 90;
ax.YTick = 1:numel(emo_lab);
ax.YTickLabel = emo_lab;

%% Plot bars
ch_lab = 'AMG1';
ch_ix  = 1;
data   = beta(2:end,ch_ix);    % don't include offset/intercept
fig_name = [SBJ '_' ch_lab '_bar_betas'];
f = figure('Name',fig_name);
ax = gca;

b = bar(1:numel(data),data,'FaceColor','flat');%numel(emo_lab):-1:1

cmap = redblue(100);%parula
cmap_idx = linspace(-max(abs(data)),max(abs(data)),size(cmap,1));
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

ax.Title.String = 'Emotional Encoding in the Amygdala';

view([90 -90]);
set(ax,'FontSize',14);
set(findall(f,'type','text'),'FontSize',16);

%% Save figure
fig_fname = [root_dir 'emodim/results/' fig_name '.png'];
fprintf('Saving %s\n',fig_fname);
saveas(gcf,fig_fname);

