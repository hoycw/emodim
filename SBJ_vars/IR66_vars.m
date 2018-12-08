%% IR66 Processing Variables
[root_dir, app_dir] = fn_get_root_dir();
if isempty(strfind(path,'fieldtrip'))
    addpath([app_dir 'fieldtrip/']);
    ft_defaults
end

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ        = 'IR66';
SBJ_vars.raw_file   = {'2017121913_0025_alan.besa'};
SBJ_vars.block_name = {''};
SBJ_vars.restart    = {1};

SBJ_vars.dirs.SBJ     = [root_dir 'emodim/data/' SBJ_vars.SBJ '/'];
SBJ_vars.dirs.raw     = [SBJ_vars.dirs.SBJ '00_raw/'];
SBJ_vars.dirs.import  = [SBJ_vars.dirs.SBJ '01_import/'];
SBJ_vars.dirs.preproc = [SBJ_vars.dirs.SBJ '02_preproc/'];
SBJ_vars.dirs.events  = [SBJ_vars.dirs.SBJ '03_events/'];
SBJ_vars.dirs.proc    = [SBJ_vars.dirs.SBJ '04_proc/'];
SBJ_vars.dirs.recon   = [SBJ_vars.dirs.SBJ '05_recon/'];
if ~exist(SBJ_vars.dirs.import,'dir')
    mkdir(SBJ_vars.dirs.import);
end
if ~exist(SBJ_vars.dirs.preproc,'dir')
    mkdir(SBJ_vars.dirs.preproc);
end
if ~exist(SBJ_vars.dirs.events,'dir')
    mkdir(SBJ_vars.dirs.events);
end
if ~exist(SBJ_vars.dirs.proc,'dir')
    mkdir(SBJ_vars.dirs.proc);
end
if ~exist(SBJ_vars.dirs.recon,'dir')
    mkdir(SBJ_vars.dirs.recon);
end

SBJ_vars.dirs.raw_filename = strcat(SBJ_vars.dirs.raw,SBJ_vars.raw_file);

SBJ_vars.recon.surf_l     = [SBJ_vars.dirs.recon 'Surfaces/' SBJ_vars.SBJ '_cortex_lh.mat'];
SBJ_vars.recon.surf_r     = [SBJ_vars.dirs.recon 'Surfaces/' SBJ_vars.SBJ '_cortex_rh.mat'];
SBJ_vars.recon.elec_pat   = [SBJ_vars.dirs.recon 'Electrodes/' SBJ_vars.SBJ '_elec_acpc_f.mat'];
SBJ_vars.recon.elec_mni_v = [SBJ_vars.dirs.recon 'Electrodes/' SBJ_vars.SBJ '_elec_mni_frv.mat'];
SBJ_vars.recon.elec_mni_s = [];
SBJ_vars.recon.fs_T1      = [SBJ_vars.dirs.recon 'Scans/' SBJ_vars.SBJ '_fs_postop_T1.mgz'];
SBJ_vars.recon.fs_DK      = [SBJ_vars.dirs.recon 'Scans/' SBJ_vars.SBJ '_fs_postop_aparc+aseg.mgz'];
SBJ_vars.recon.fs_Dx      = [SBJ_vars.dirs.recon 'Scans/' SBJ_vars.SBJ '_fs_postop_aparc.a2009s+aseg.mgz'];

%--------------------------------------
% Channel Selection
%--------------------------------------
%hdr = ft_read_header(SBJ_vars.dirs.raw_filename);
%SBJ_vars.orig_n_ch = length(hdr.label);
%SBJ_vars.orig_n_samples = hdr.nSamples;
%SBJ_vars.orig_srate = hdr.Fs;
%clear hdr;

SBJ_vars.ch_lab.probes     = {'RAM','RHH','RTH','RAC','ROF','RPC','LAM','LHH','LTH','LAC','LOF','LPC'};
SBJ_vars.ch_lab.probe_type = {'SEEG','SEEG','SEEG','SEEG','SEEG','SEEG','SEEG','SEEG','SEEG','SEEG','SEEG','SEEG'};
SBJ_vars.ch_lab.ref_type   = {'BP','BP','BP','BP','BP','BP','BP','BP','BP','BP','BP','BP'};
SBJ_vars.ch_lab.ROI        = {'all'};
SBJ_vars.ch_lab.eeg_ROI    = {};

%SBJ_vars.ch_lab.prefix = 'POL ';    % before every channel except 'EDF Annotations'
%SBJ_vars.ch_lab.suffix = '-Ref';    % after every channel except 'EDF Annotations'
SBJ_vars.ch_lab.mislabel = {{'CZ ','CZ'}};

SBJ_vars.ref_exclude = {}; %exclude from the CAR
SBJ_vars.ch_lab.bad = {...
    'LHH1','LHH2','LHH3','LTH1','LTH2','LTH3','LTH4',...% epileptic
    'RHH1','RHH2','RHH3','RHH4','RTH1','RTH2','RTH3',...% epileptic, trying to save RAM1, RAM2
    'LHH4','LHH10','LTH5','LTH6',...% tossed after looking at preproc
    'LAM1','LAM2','LAM3','LAM4','LAM5','LAM6',...% tossed for spread after looking at preproc
    'RAM1','RAM2',...% tossed after looking at preproc
    'LAC8','LAC9','LAC10',...% out of brain
    'DC01','DC03','DC04',...% empty analogs
    'GND','REF','EKG'...% not real data
    };
SBJ_vars.ch_lab.eeg = {'FZ','CZ ','OZ','C3','C4'};
SBJ_vars.ch_lab.eog = {'LUE','LLE','RUE','RLE'};
SBJ_vars.ch_lab.photod = {'DC02'};

%--------------------------------------
% Line Noise Parameters
%--------------------------------------
SBJ_vars.notch_freqs = [60 120 180 240 300];
SBJ_vars.bs_width    = 2;

%--------------------------------------
% Time Parameters
%--------------------------------------
SBJ_vars.analysis_time = {{[20 1300]}};

%--------------------------------------
% Artifact Rejection Parameters
%--------------------------------------
SBJ_vars.artifact_params.std_limit_raw = 7;
SBJ_vars.artifact_params.hard_threshold_raw = 700;

SBJ_vars.artifact_params.std_limit_diff = 7;
SBJ_vars.artifact_params.hard_threshold_diff = 70;

%--------------------------------------
% Trials to Reject
%--------------------------------------
% SBJ_vars.trial_reject_n = [];
