function SBJ06_HFA_save(SBJ,pipeline_id,an_id)
% Calculates high frequency activity and saves fieldtrip struct

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'hoycw/Apps/fieldtrip/'];
elseif exist('/Users/lapate/','dir');root_dir = '/Users/lapate/knight/';ft_dir = '/Users/lapate/knight/hoycw/Apps/fieldtrip/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%% Set Up Directories
addpath([root_dir 'emodim/scripts/']);
addpath([root_dir 'emodim/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Data Preparation
SBJ_vars_cmd = ['run ' root_dir 'emodim/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
an_vars_cmd = ['run ' root_dir 'emodim/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);

% Load Data
load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',pipeline_id,'.mat'));
% load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));

%% Select Channel(s)
cfgs = [];
cfgs.channel = SBJ_vars.ch_lab.ROI;
roi = ft_selectdata(cfgs,data);

%% Compute HFA
fprintf('===================================================\n');
fprintf('------------------ HFA Calculations ---------------\n');
fprintf('===================================================\n');
cfg_hfa.trials = 'all';
if strcmp(HFA_type,'multiband')
    hfa = ft_freqanalysis(cfg_hfa, roi);
elseif any(strcmp(HFA_type,{'broadband','hilbert'}))
    % Hilbert method to extract power
    if strcmp(bsln_type,'zboot')
        error('Treating filter-hilbert output as "trials" is incompatible with baslines smapled across trials: zboot');
    end
    hfas = cell(size(fois));
    for f_ix = 1:numel(fois)
        cfg_hfa.bpfreq = bp_lim(f_ix,:);
        fprintf('\n------> %s filtering: %.03f - %.03f\n', HFA_type, bp_lim(f_ix,1), bp_lim(f_ix,2));
        hfas{f_ix} = ft_preprocessing(cfg_hfa,roi);
    end
    % Treat different freqs as trials
    hfa = ft_appenddata([], hfas{:}); clear hfas;
else
    error('Unknown HFA_type provided');
end

%% Baseline Correction
fprintf('===================================================\n');
fprintf('---------------- Baseline Correction --------------\n');
fprintf('===================================================\n');
switch bsln_type
    case {'zboot', 'zscore'}
        if strcmp(HFA_type,'multiband')
            hfa = fn_bsln_ft_tfr(hfa,bsln_lim,bsln_type,n_boots);
        elseif any(strcmp(HFA_type,{'broadband','hilbert'}))
            hfa = fn_bsln_ft_filtered(hfa,bsln_lim,bsln_type,n_boots);
        else
            error('Unknown HFA_type provided');
        end
    case {'relchange', 'demean', 'my_relchange'}
        error(['bsln_type ' bsln_type ' is not compatible with one-sample t test bsln activation stats']);
%         cfgbsln = [];
%         cfgbsln.baseline     = bsln_lim;
%         cfgbsln.baselinetype = bsln_type;
%         cfgbsln.parameter    = 'powspctrm';
%         hfa = ft_freqbaseline(cfgbsln,hfa);
    otherwise
        error(['No baseline implemented for bsln_type: ' bsln_type]);
end

%% Smooth Power Time Series
if smooth_pow_ts
    fprintf('===================================================\n');
    fprintf('----------------- Filtering Power -----------------\n');
    fprintf('===================================================\n');
    if ~strcmp(HFA_type,'multiband')
        error('HFA filtering not implemented for hilbert-filter data yet!');
    end
    for ch_ix = 1:numel(hfa.label)
        for f_ix = 1:numel(fois)
            if strcmp(lp_yn,'yes') && strcmp(hp_yn,'yes')
                hfa.powspctrm(:,ch_ix,f_ix,:) = fn_EEGlab_bandpass(...
                    hfa.powspctrm(:,ch_ix,f_ix,:), roi.fsample, hp_freq, lp_freq);
            elseif strcmp(lp_yn,'yes')
                hfa.powspctrm(:,ch_ix,f_ix,:) = fn_EEGlab_lowpass(...
                    squeeze(hfa.powspctrm(:,ch_ix,f_ix,:)), roi.fsample, lp_freq);
            elseif strcmp(hp_yn,'yes')
                error('Why are you only high passing?');
            else
                error('Why did you say yes smooth but no to both low and high pass?');
            end
        end
    end
end

%% Merge multiple bands
cfg_avg = [];
if strcmp(HFA_type,'multiband')
    cfg_avg.freq = 'all';
    cfg_avg.avgoverfreq = 'yes';
elseif any(strcmp(HFA_type,{'broadband','hilbert'}))
    cfg_avg.trials = 'all';
    cfg_avg.avgoverrpt = 'yes';
end
hfa = ft_selectdata(cfg_avg,hfa);

%% Save Results
data_out_filename = strcat(SBJ_vars.dirs.proc,SBJ,'_',an_id,'.mat');
fprintf('===================================================\n');
fprintf('--- Saving %s ------------------\n',data_out_filename);
fprintf('===================================================\n');
save(data_out_filename,'-v7.3','hfa');

end
