if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'hoycw/Apps/fieldtrip/'];
elseif exist('/Users/lapate/','dir');root_dir = '/Users/lapate/knight/';ft_dir = '/Users/lapate/knight/hoycw/Apps/fieldtrip/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%% Set Up Directories
addpath([root_dir 'emodim/scripts/']);
addpath([root_dir 'emodim/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Channel Label Assignment
bad_lab = {'LTH1','LTH2','LHH1','LHH2','RTH4','RTH5','RHH3','RHH4'};
sus_lab = {'LAM1','LAM2','LHH3','RAM2','RTH2','RTH3'};
out_lab = {'LOF10','LAC10','RAM10','RHH1','RAC10','ROF10','RAM9','LAM9','LAM10','LHH10'};
jnk_lab = {'---(13)','---(14)','---(17)','---(18)','FZ','CZ','C3','C4','OZ','FPZ','Z',...
            'EKG','REF','DC01','DC02','DC03','DC04'};

%% Load data 
load([root_dir 'emodim/scripts/utils/cfg_plot.mat']);

%% Toss junk

jnk_neg = fn_ch_lab_negate(jnk_lab);
cfgs = [];
cfgs.channel = {'all',jnk_neg{:}};
data = ft_selectdata(cfgs,data);

%% Assign colors
bad_color = [1 0 0];
sus_color = [0.7 0 0.7];
out_color = [0 0.1 1];
all_color = [0.3 0.3 0.4];
colormap  = [bad_color; sus_color; out_color; all_color];
cfg_plot.channelcolormap = colormap;


% Add color groups as first 2 characters of labels
% make cfg_plot.colorgroups to be  'labelchar2'

cgroups = zeros([numel(data.label) 1]);
colors = zeros([numel(data.label) 3]);
for l = 1:numel(data.label)
    if numel(data.label{l})<3
        data.label{l} = [data.label{l} '  '];
    end
    
    if any(strcmp(data.label{l},bad_lab))
        cgroups(l)  = 1;
        colors(l,:) = bad_color;
    elseif any(strcmp(data.label{l},sus_lab))
        cgroups(l)  = 2;
        colors(l,:) = sus_color;
    elseif any(strcmp(data.label{l},out_lab))
        cgroups(l)  = 3;
        colors(l,:) = out_color;
    else
        cgroups(l)  = 4;
        colors(l,:) = all_color;
    end
    fprintf('%s \t- %i %.01f %.01f %.01f\n',data.label{l},cgroups(l),colors(l,:));
end
cfg_plot.colorgroups = cgroups;

% for l = 1:numel(bad_lab)
%     ix = strcmp(data.label,bad_lab{l});
%     data.label{ix} = ['bad' data.label{ix}];
% end
% for l = 1:numel(sus_lab)
%     ix = strcmp(data.label,sus_lab{l});
%     data.label{ix} = ['sus' data.label{ix}];
% end
% % for l = 1:numel(out_lab)
% %     ix = strcmp(data.label,out_lab{l});
% %     data.label{ix} = ['sus' data.label{ix}];
% % end
% 
% % not_real_neg = fn_ch_lab_negate(not_real);
% % eeg_neg = fn_ch_lab_negate({'FPZ','CZ','OZ','C3','C4'});
% % eog_neg = fn_ch_lab_negate({'LUC','LLC','RUC','RLC'});
% % cfgs = [];
% % cfgs.channel = {'all',not_real_neg{:},eeg_neg{:},eog_neg{:}};
% % data = ft_selectdata(cfgs,data);
% 
