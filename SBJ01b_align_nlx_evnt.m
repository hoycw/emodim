function SBJ01b_align_nlx_evnt(SBJ, pipeline_id, save_it)
% save_it == 0: don't save plots, compare to raw data
%         == 1: save plots and data, compare to import

%% Load, preprocess, and save out photodiode
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'hoycw/Apps/'];
elseif exist('/Users/lapate/','dir');root_dir = '/Users/lapate/knight/';app_dir = '/Users/lapate/knight/hoycw/Apps/';
else root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

%% Paths
addpath([root_dir 'emodim/scripts/']);
addpath([root_dir 'emodim/scripts/utils/']);
addpath(genpath([app_dir 'wave_clus/']));
addpath(genpath([app_dir 'UR_NLX2MAT_releaseDec2015/']));
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% SBJ vars
eval(['run ' root_dir 'emodim/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run ' root_dir 'emodim/scripts/proc_vars/' pipeline_id '_proc_vars.m']);

%% Read photodiode, NLX macro, clinical data
% Neuralynx photodiode
evnt       = ft_read_neuralynx_interp({[SBJ_vars.dirs.nlx 'photo/' ...
                            SBJ_vars.ch_lab.photod{1} SBJ_vars.ch_lab.nlx_suffix '.ncs']});
evnt.label = {'photo'};
evnt_orig  = evnt;

% Neuralynx macro channel
macro_fnames = SBJ_vars.ch_lab.nlx_nk_align;
for m_ix = 1:numel(macro_fnames)
    macro_fnames{m_ix} = [SBJ_vars.dirs.nlx 'macro/' SBJ_vars.ch_lab.nlx_nk_align{m_ix}...
                            SBJ_vars.ch_lab.nlx_suffix '.ncs'];
end
macro = ft_read_neuralynx_interp(macro_fnames);
macro.label = SBJ_vars.ch_lab.nlx_nk_align;
macro_orig  = macro;

% Nihon Kohden clinical channel
if ~save_it
    clin_fname = SBJ_vars.dirs.raw_filename{b_ix};
else
    if any(SBJ_vars.low_srate)
        srate_str = num2str(SBJ_vars.low_srate(b_ix));
    else
        srate_str = num2str(proc_vars.resample_freq);
    end
    if numel(SBJ_vars.raw_file)>1
        data_blocks = {};
        for b_ix = 1:numel(SBJ_vars.block_name)
            load([SBJ_vars.dirs.import SBJ '_' srate_str 'hz_' SBJ_vars.block_name{b_ix} '.mat']);
            data_blocks{b_ix} = data;
        end
        data = fn_concat_blocks(data_blocks);
    else
        block_suffix = SBJ_vars.block_name{b_ix};   % should just be ''
        clin_fname = [SBJ_vars.dirs.import SBJ '_' srate_str 'hz' block_suffix '.mat'];
        load(clin_fname);
    end
end

cfgs         = [];
cfgs.channel = SBJ_vars.ch_lab.nlx_nk_align;
clin         = ft_selectdata(cfgs,data);
clin_orig    = clin;

%% Preprocess
% Cut NLX to nlx_analysis_time
if isfield(SBJ_vars,'nlx_analysis_time')
    cfgs = [];
    cfgs.latency = SBJ_vars.nlx_analysis_time;
    evnt = ft_selectdata(cfgs, evnt);
    macro = ft_selectdata(cfgs, macro);
end
if any(isnan(evnt.trial{1}(:))) || any(isnan(macro.trial{1}(:)))
    error('NaNs detected in NLX data, check for discontinuities!');
end

% Inversion on NLX data
if SBJ_vars.nlx_macro_inverted
    macro.trial{1} = macro.trial{1}*-1;
end
if SBJ_vars.photo_inverted
    evnt.trial{1} = evnt.trial{1}*-1;
end

% Bipolar rereference
if numel(SBJ_vars.ch_lab.nlx_nk_align)>1
    cfg = [];
    cfg.channel = macro.label;
    cfg.montage.labelold = cfg.channel;
    [cfg.montage.labelnew, cfg.montage.tra, ~] = fn_create_ref_scheme_bipolar(SBJ_vars.ch_lab.nlx_nk_align);
    macro = ft_preprocessing(cfg,macro);
    
    cfg = [];
    cfg.channel = clin.label;
    cfg.montage.labelold = cfg.channel;
    [cfg.montage.labelnew, cfg.montage.tra, ~] = fn_create_ref_scheme_bipolar(SBJ_vars.ch_lab.nlx_nk_align);
    clin = ft_preprocessing(cfg,clin);
end

% % Cut clincial macro to analysis time
% if numel(SBJ_vars.analysis_time{b_ix})>1
%     error('havent set up processing for multi block concat!');
% end
% cfgs = []; cfgs.latency = SBJ_vars.analysis_time{b_ix}{1};
% clin = ft_selectdata(cfgs,clin);
% clin.time{1} = clin.time{1}-SBJ_vars.analysis_time{b_ix}{1}(1);

% Match sampling rates
if macro.fsample > clin.fsample
    fprintf('downsampling Neuralynx from %d to %d Hz\n', macro.fsample, clin.fsample)
    cfgr = [];
    cfgr.demean     = 'yes';
    cfgr.resamplefs = clin.fsample;
    macro = ft_resampledata(cfgr, macro);
elseif macro.fsample < clin.fsample
    error('Clinical data has higher sampling rate than NLX macro, whats going on?');
end
if evnt.fsample > clin.fsample
    fprintf('downsampling Neuralynx photodiode from %d to %d Hz\n', evnt.fsample, clin.fsample)
    cfgr = [];
    cfgr.demean     = 'yes';
    cfgr.resamplefs = clin.fsample;
    evnt = ft_resampledata(cfgr, evnt);
end

% % Remove extreme values
% clin_thresh  = proc_vars.nlx_nk_align_std_thresh*std(clin.trial{1});
% macro_thresh = proc_vars.nlx_nk_align_std_thresh*std(macro.trial{1});
% clin.trial{1}((clin.trial{1}>median(clin.trial{1})+clin_thresh)|(clin.trial{1}<median(clin.trial{1})-clin_thresh)) = median(clin.trial{1});
% macro.trial{1}((macro.trial{1}>median(macro.trial{1})+macro_thresh)|(macro.trial{1}<median(macro.trial{1})-macro_thresh)) = median(macro.trial{1});

%% Compare PSDs
fn_plot_PSD_1by1_compare(clin.trial{1},macro.trial{1},clin.label,macro.label,...
    clin.fsample,'clinical','macro');
if save_it
    saveas(gcf,[SBJ_vars.dirs.import SBJ '_nlx_nk_PSD_compare_' macro.label{1} '.png']);
end

%% Compute cross correlation at varying time lags
if numel(clin.trial{1}) <= numel(macro.trial{1})
    warning('Clinical data is smaller than macro data, recut clinical block!');
end

% Arjen's way: synchronize nihon kohden and neuralynx timeseries
%   Find cross-variance and lags by shifting macro along clin
[covar, lags] = xcov(clin.trial{1}', macro.trial{1}');
%   Find the peak (previously: get sharp peak; remove very long trends with
%       smooth (moving average) to find big spike in covariance)
[~, idx]      = max(covar);% - smooth(covar, clin.fsample*10));   
if isfield(SBJ_vars,'nlx_nk_align_force')
    algorithm_idx = idx;
    idx = SBJ_vars.nlx_nk_align_force;
end

%% Plot match
figure; subplot(2,1,1);
hold on; plot(lags,covar);
hold on; plot(lags(idx),covar(idx),'k*');
if exist('algorithm_idx','var')
    plot(lags(algorithm_idx),covar(algorithm_idx),'r*');
    legend('corrected','alogrithm');
else
    legend('algorithm');
end
ylabel('correlation');
xlabel('lag');
subplot(2,1,2);
t = 1:numel(clin.time{1}); 
hold on; plot(t, zscore(clin.trial{1}));
t2 = lags(idx):lags(idx)+numel(macro.trial{1})-1;
hold on; plot(t2, zscore(macro.trial{1})+10);
t3 = lags(idx):lags(idx)+numel(evnt.time{1})-1; % ignore any offset between photo and chan
hold on; plot(t3, zscore(evnt.trial{1})+20);
legend('NK', 'NLX', 'NLX photo');
title(macro.label{1});

%% Save figure and data
if save_it
    saveas(gcf,[SBJ_vars.dirs.import SBJ '_nlx_nk_macro_alignment_' macro.label{1} '.fig']);
    % print([subj(1).datadir 'datafiles/sync_nk-nl_' subjectm(9:end) '_' num2str(tcgver) '_' num2str(d)], '-dpdf');
    
    %% Create photodiode channel matched to clinical data
    evnt_nlx   = evnt;
    evnt       = clin;
    evnt.label = SBJ_vars.ch_lab.photod;
    % Create dummy time series of the median of the photodiode
    evnt.trial{1} = ones(1,numel(clin.trial{1})).*median(evnt.trial{1});
    % Add in photodiode data for segments when NLX overlaps
    evnt.trial{1}(t3(t3>0 & t3<numel(evnt.trial{1}))) = evnt_nlx.trial{1}(t3>0 & t3<numel(evnt.trial{1}));
    
    %% Save data out
    evnt_out_filename = strcat(SBJ_vars.dirs.import,SBJ,'_evnt.mat');
    save(evnt_out_filename, '-v7.3', 'evnt');
end

end