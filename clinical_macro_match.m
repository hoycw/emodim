%% Load, preprocess, and save out photodiode
SBJ = 'IR84';
inverted = 1;

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
b_ix = 1;   %block
eval(['run ' root_dir 'emodim/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run ' root_dir 'emodim/scripts/proc_vars/SU_nlx_proc_vars.m']);

if numel(SBJ_vars.raw_file)>1
    block_suffix = strcat('_',SBJ_vars.block_name{b_ix});
else
    block_suffix = SBJ_vars.block_name{b_ix};   % should just be ''
end

%% Read photodiode, NLX macro, clinical data
photo = ft_read_neuralynx_interp({[SBJ_vars.dirs.nlx 'photo/' ...
                            SBJ_vars.ch_lab.photod{1} SBJ_vars.ch_lab.nlx_suffix '.ncs']});
photo.label = {'photo'};
photo_orig  = photo;

macro = ft_read_neuralynx_interp({[SBJ_vars.dirs.nlx 'macro/' ...
                            SBJ_vars.ch_lab.nlx_nk_align{1} SBJ_vars.ch_lab.nlx_suffix '.ncs']});
macro.label = SBJ_vars.ch_lab.nlx_nk_align;
macro_orig  = macro;

load([SBJ_vars.dirs.raw SBJ '_raw_emodim_clinical.mat']);
cfgs = [];
cfgs.channel = SBJ_vars.ch_lab.nlx_nk_align;
clin         = ft_selectdata(cfgs,data);
clin_orig    = clin;

%% Preprocess
% Inversion on NLX data
if SBJ_vars.nlx_inverted
    photo.trial{1} = photo.trial{1}*-1;
    macro.trial{1} = macro.trial{1}*-1;
end

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
if photo.fsample > clin.fsample
    fprintf('downsampling Neuralynx photodiode from %d to %d Hz\n', photo.fsample, clin.fsample)
    cfgr = [];
    cfgr.demean     = 'yes';
    cfgr.resamplefs = clin.fsample;
    photo = ft_resampledata(cfgr, photo);
end

% Remove extreme values
std_thresh   = 4;
clin_thresh  = std_thresh*std(clin.trial{1});
macro_thresh = std_thresh*std(macro.trial{1});
clin.trial{1}((clin.trial{1}>median(clin.trial{1})+clin_thresh)|(clin.trial{1}<median(clin.trial{1})-clin_thresh)) = median(clin.trial{1});
macro.trial{1}((macro.trial{1}>median(macro.trial{1})+macro_thresh)|(macro.trial{1}<median(macro.trial{1})-macro_thresh)) = median(macro.trial{1});

% if z_transform
%     macro_orig = macro;
%     clin_orig = clin;
%     macro.trial{1} = (macro.trial{1}-mean(macro.trial{1}))/std(macro.trial{1});
%     clin.trial{1} = (clin.trial{1}-mean(clin.trial{1}))/std(clin.trial{1});
% end

%% Compare PSDs
fn_plot_PSD_1by1_compare(clin.trial{1},macro.trial{1},clin.label,macro.label,...
    clin.fsample,'clinical','macro');
saveas(gcf,[SBJ_vars.dirs.import SBJ '_nlx_nk_PSD_compare.png']);

%% Plot TS
plot_len = 10000;
figure;
subplot(2,1,1);
plot(clin.time{1}(1:plot_len),clin.trial{1}(1:plot_len));
subplot(2,1,2);
plot(macro.time{1}(1:plot_len),macro.trial{1}(1:plot_len));

%% Compute cross correlation at varying time lags
n_clin  = numel(clin.trial{1});
n_macro = numel(macro.trial{1});
if n_clin <= n_macro
    error('Clinical data is smaller than macro data, recut clinical block!');
end

% Arjen's way
% synchronize nihon kohden and neuralynx timeseries
% Find cross-variance and lags by shifting macro along clin
[covar, lags] = xcov(clin.trial{1}', macro.trial{1}');
% find the sharp peak; remove very long trends with smooth (moving average) to find big spike in covariance
[~, idx]      = max(covar - smooth(covar, clin.fsample*10));   

%% Plot match
figure; subplot(2,1,1);
hold on; plot(lags,covar);
hold on; plot(lags(idx),covar(idx),'k*');
ylabel('correlation');
xlabel('lag');
subplot(2,1,2);
t = 1:numel(clin.time{1}); 
hold on; plot(t, zscore(clin.trial{1}));
t2 = lags(idx):lags(idx)+numel(macro.trial{1})-1;
hold on; plot(t2, zscore(macro.trial{1})+10);
t3 = lags(idx):lags(idx)+numel(photo.time{1})-1; % ignore any offset between photo and chan
hold on; plot(t3, zscore(photo.trial{1})+20);
legend('NK', 'NLX', 'NLX photo');
saveas(gcf,[SBJ_vars.dirs.import SBJ '_nlx_nk_macro_alignment.fig']);
% print([subj(1).datadir 'datafiles/sync_nk-nl_' subjectm(9:end) '_' num2str(tcgver) '_' num2str(d)], '-dpdf');

%% create nihon kohden photodiode channel
photo_nlx = photo;
photo = clin;
photo.label = {'photodiode'};
% Create dummy time series of the median of the photodiode
photo.trial{1} = ones(1,numel(clin.trial{1})).*median(photo.trial{1});
% Add in photodiode data for segments when NLX overlaps
photo.trial{1}(t3(t3>0 & t3<numel(photo.trial{1}))) = photo_nlx.trial{1}(t3>0 & t3<numel(photo.trial{1}));

% Save data out
photo_fname = [SBJ_vars.dirs.import SBJ '_photodiode_nk_aligned.mat'];
save(photo_fname, 'photo');


%% OLD WAY:
% Colin original way:
% corr_init_ix = floor(1:proc_vars.align_resamp_init*match_srate:n_clin-n_macro);
% corr_init = zeros(size(corr_init_ix));
% for ix = 1:numel(corr_init_ix)
%     tmp_corr = corrcoef(clin.trial{1}(corr_init_ix(ix):corr_init_ix(ix)+n_macro-1),...
%                             macro.trial{1});
%     corr_init(ix) = tmp_corr(1,2);
% end
% 
% [max_init, max_init_ix] = max(abs(corr_init));
% figure; hold on;
% plot(corr_init);
% scatter(max_init_ix,corr_init(max_init_ix),'r');
% 
% % Second level
% corr_full_ix = corr_init_ix(max_init_ix-1):corr_init_ix(max_init_ix+1);
% corr_full = zeros(size(corr_full_ix));
% for ix = 1:numel(corr_full_ix)
%     tmp_corr = corrcoef(clin.trial{1}(corr_full_ix(ix):corr_full_ix(ix)+n_macro-1),...
%                             macro.trial{1});
%     corr_full(ix) = tmp_corr(1,2);
% end
% 
% [max_full, max_full_ix] = max(abs(corr_full));
% figure; hold on;
% plot(corr_full);
% scatter(max_full_ix,corr_full(max_full_ix),'r');


% match_ix = corr_full_ix(max_full_ix);
% match_ts = clin.trial{1}(1,match_ix:match_ix+n_macro-1);
% match_ts_z = (match_ts-mean(match_ts))/std(match_ts);
% macro_z = (macro.trial{1}-mean(macro.trial{1}))/std(macro.trial{1});
% figure; hold on;
% plot(macro.time{1},macro_z*-1,'b');
% plot(macro.time{1},match_ts_z,'r');

