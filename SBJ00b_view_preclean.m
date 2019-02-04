function SBJ00b_view_preclean(SBJ)

%% Check which root directory
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'hoycw/Apps/fieldtrip/'];
elseif exist('/Users/lapate/','dir');root_dir = '/Users/lapate/knight/';app_dir = '/Users/lapate/knight/hoycw/Apps/fieldtrip/';
else root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%% Set Up Directories
addpath([root_dir 'emodim/scripts/']);
addpath([root_dir 'emodim/scripts/utils/']);
addpath(genpath([app_dir 'wave_clus/']));
addpath(genpath('/Users/colinhoy/Code/Apps/UR_NLX2MAT_releaseDec2015/'));
addpath([app_dir 'fieldtrip/']);
ft_defaults


%% ========================================================================
SBJ_vars_cmd = ['run ' root_dir 'emodim/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

%% Process channel labels
% Keep any real but "bad" channels
toss = ones(size(SBJ_vars.ch_lab.bad));
for bad_ix = 1:numel(SBJ_vars.ch_lab.bad)
    ch_match = zeros(size(SBJ_vars.ch_lab.probes));
    for probe_ix = 1:numel(SBJ_vars.ch_lab.probes)
        if ~isempty(strfind(SBJ_vars.ch_lab.bad{bad_ix},SBJ_vars.ch_lab.probes{probe_ix}))
            ch_match(probe_ix) = 1;
        end
    end
    if any(ch_match)
        toss(bad_ix) = 0;
    end
end
SBJ_vars.ch_lab.bad = SBJ_vars.ch_lab.bad(logical(toss));

% Handle prefix and suffix
if isfield(SBJ_vars.ch_lab,'prefix')
    for bad_ix = 1:numel(SBJ_vars.ch_lab.bad)        
        SBJ_vars.ch_lab.bad{bad_ix} = [SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.bad{bad_ix}];
    end
    for eeg_ix = 1:numel(SBJ_vars.ch_lab.eeg)
        SBJ_vars.ch_lab.eeg{eeg_ix} = [SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.eeg{eeg_ix}];
    end
    for eog_ix = 1:numel(SBJ_vars.ch_lab.eog)
        SBJ_vars.ch_lab.eog{eog_ix} = [SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.eog{eog_ix}];
    end
    SBJ_vars.ch_lab.photod = {[SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.photod{1}]};
end
if isfield(SBJ_vars.ch_lab,'suffix')
    for bad_ix = 1:numel(SBJ_vars.ch_lab.bad)
        SBJ_vars.ch_lab.bad{bad_ix} = [SBJ_vars.ch_lab.bad{bad_ix} SBJ_vars.ch_lab.suffix];
    end
    for eeg_ix = 1:numel(SBJ_vars.ch_lab.eeg)
        SBJ_vars.ch_lab.eeg{eeg_ix} = [SBJ_vars.ch_lab.eeg{eeg_ix} SBJ_vars.ch_lab.suffix];
    end
    for eog_ix = 1:numel(SBJ_vars.ch_lab.eog)
        SBJ_vars.ch_lab.eog{eog_ix} = [SBJ_vars.ch_lab.eog{eog_ix} SBJ_vars.ch_lab.suffix];
    end
    SBJ_vars.ch_lab.photod = {[SBJ_vars.ch_lab.photod{1} SBJ_vars.ch_lab.suffix]};
end
bad_ch_neg = fn_ch_lab_negate(SBJ_vars.ch_lab.bad);
eeg_ch_neg = fn_ch_lab_negate(SBJ_vars.ch_lab.eeg);
eog_ch_neg = fn_ch_lab_negate(SBJ_vars.ch_lab.eog);
photod_ch_neg = fn_ch_lab_negate(SBJ_vars.ch_lab.photod);

%% Load data
b_ix = 1;

if numel(SBJ_vars.raw_file)>1
    block_suffix = strcat('_',SBJ_vars.block_name{b_ix});
else
    block_suffix = SBJ_vars.block_name{b_ix};   % should just be ''
end

if strcmp(SBJ_vars.raw_file{b_ix}(end-2:end),'mat')
    orig = load(SBJ_vars.dirs.raw_filename{b_ix});
    data = orig.data;
else
    % Load Neural Data
    cfg            = [];
    cfg.dataset    = SBJ_vars.dirs.raw_filename{b_ix};
    cfg.continuous = 'yes';
%     cfg.channel    = {'all',eeg_ch_neg{:},eog_ch_neg{:},photod_ch_neg{:}};%bad_ch_neg{:},
    data = ft_preprocessing(cfg);
end
data_orig = data;

%% Select Neural Data
cfg = [];
cfg.channel = {'all',bad_ch_neg{:},eeg_ch_neg{:},eog_ch_neg{:},photod_ch_neg{:}};
data = ft_selectdata(cfg,data);

%% Reorder
data = fn_reorder_data(data, {});

%% Plot
load([root_dir 'emodim/scripts/utils/cfg_plot.mat']);

out = ft_databrowser(cfg_plot,data);
bad_epochs = out.artfctdef.visual.artifact;

%% Colorize
% edit fn_databrowser_colors.m

%% Save
save([SBJ_vars.dirs.events SBJ '_bob_bad_epochs_preclean.mat'],'-v7.3','bad_epochs');
end