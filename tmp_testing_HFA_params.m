if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'hoycw/Apps/fieldtrip/'];
elseif exist('/Users/lapate/','dir');root_dir = '/Users/lapate/knight/';ft_dir = '/Users/lapate/knight/hoycw/Apps/fieldtrip/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%% Set Up Directories
addpath([root_dir 'emodim/scripts/']);
addpath([root_dir 'emodim/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%%
SBJ = 'IR71';
pipeline_id = 'main_ft';

% Data Preparation
SBJ_vars_cmd = ['run ' root_dir 'emodim/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

% Load Data
load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',pipeline_id,'.mat'));
cfgs = [];
cfgs.channel = {'AM6-7'};
roi = ft_selectdata(cfgs,data);

load([root_dir 'emodim/scripts/utils/cfg_plot.mat']);
cfg_plot.ylim = [-3 3];

%%
an_id = 'HGh_zsc1to8_sm0_fLog';
an_vars_cmd = ['run ' root_dir 'emodim/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
cfghil = cfg_hfa;
cfghil.trials = 'all';

an_id = 'HGm_zsc1to8_sm0_wn100';
an_vars_cmd = ['run ' root_dir 'emodim/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
cfgorig = cfg_hfa;
cfgorig.trials = 'all';

an_id = 'HGm_zsc1to8_sm0_wn250';
an_vars_cmd = ['run ' root_dir 'emodim/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
cfglong = cfg_hfa;
cfglong.trials = 'all';

an_id = 'HGm_zsc1to8_oct12_wn100';
an_vars_cmd = ['run ' root_dir 'emodim/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
cfgoct = cfg_hfa;
cfgoct.trials = 'all';

%%
orig = ft_freqanalysis(cfgorig,roi);
orig = fn_bsln_ft_tfr(orig,bsln_lim,bsln_type,n_boots);

long = ft_freqanalysis(cfglong,roi);
long = fn_bsln_ft_tfr(long,bsln_lim,bsln_type,n_boots);
% has up to 19 tapers!

oct = ft_freqanalysis(cfgoct,roi);
oct = fn_bsln_ft_tfr(oct,bsln_lim,bsln_type,n_boots);

%%
% split = {};
% for f = 1:numel(cfg_hfa.foi)
%     cfgtmp = cfg_hfa;
%     cfgtmp.output = 'fourier';
%     cfgtmp.foi = cfg_hfa.foi(f);
%     cfgtmp.tapsmofrq = cfg_hfa.tapsmofrq(f);
%     cfgtmp.t_ftimwin = cfg_hfa.t_ftimwin(f);
%     split{f} = ft_freqanalysis(cfgtmp,roi);
% end

cfgsimp = cfg_hfa;
cfgsimp.tapsmofrq = repmat(cfg_hfa.tapsmofrq(1),size(cfg_hfa.tapsmofrq));
cfgsimp.keeptapers = 'no';
simp = ft_freqanalysis(cfgsimp,roi);
% Looks the same, with all frequencies looking nearly identical

cfghfasm = cfg_hfa;
cfghfasm.keeptapers = 'no';
cfghfasm.tapsmofrq = cfg_hfa.tapsmofrq/2;
hfasm67 = ft_freqanalysis(cfghfasm, roi);
hfasm_bsln67 = fn_bsln_ft_tfr(hfasm67,bsln_lim,bsln_type,n_boots);
% looks the same as simp, maybe a little messier at lower freq (only 1 taper in 70, 80 Hz)

cfghan = cfg_hfa;
cfghan.keeptapers = 'no';
cfghan.taper = 'hanning';
han67 = ft_freqanalysis(cfghan, roi);
han_bsln67 = fn_bsln_ft_tfr(han67,bsln_lim,bsln_type,n_boots);
% EXCELLENT! This look smuch more like what I want, without double peaks
% and very distinct content across the frequency sub-bands!

cfgwvlt = [];
cfgwvlt.method = 'wavelet';
cfgwvlt.output = 'pow';
cfgwvlt.trials = 'all';
cfgwvlt.channel = 'all';
cfgwvlt.pad = 'maxperlen';
cfgwvlt.padtype = 'zero';
cfgwvlt.foi = [70:10:150];
cfgwvlt.toi = 'all';
cfgwvlt.keeptrials = 'yes';
wvlt67 = ft_freqanalysis(cfgwvlt,roi);
wvlt_bsln67 = fn_bsln_ft_tfr(wvlt67,bsln_lim,bsln_type,n_boots);