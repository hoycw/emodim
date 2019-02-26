%% IR78 Processing Variables
[root_dir, app_dir] = fn_get_root_dir();
if isempty(strfind(path,'fieldtrip'))
    addpath([app_dir 'fieldtrip/']);
    ft_defaults
end

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ        = 'IR78';
SBJ_vars.raw_file   = {'2018062715_0046.besa'};   %_0044.besa and _0045.besa have partial B1s
SBJ_vars.block_name = {''};
SBJ_vars.restart    = {1};
SBJ_vars.low_srate  = [0];

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
SBJ_vars.recon.fs_T1      = [SBJ_vars.dirs.recon 'Scans/' SBJ_vars.SBJ '_fs_preop_T1.mgz'];
SBJ_vars.recon.fs_DK      = [SBJ_vars.dirs.recon 'Scans/' SBJ_vars.SBJ '_fs_preop_aparc+aseg.mgz'];
SBJ_vars.recon.fs_Dx      = [SBJ_vars.dirs.recon 'Scans/' SBJ_vars.SBJ '_fs_preop_aparc.a2009s+aseg.mgz'];

%--------------------------------------
% Channel Selection
%--------------------------------------
%hdr = ft_read_header(SBJ_vars.dirs.raw_filename);
%SBJ_vars.orig_n_ch = length(hdr.label);
%SBJ_vars.orig_n_samples = hdr.nSamples;
%SBJ_vars.orig_srate = hdr.Fs;
%clear hdr;

SBJ_vars.ch_lab.probes     = {'LOF','LSM','LPC','LAC','LTH','LHH','LAM','RSM','MEG','SMG','ROF','RPC','RAC','RTH','RHH','RAM'};
SBJ_vars.ch_lab.probe_type = {'SEEG','SEEG','SEEG','SEEG','SEEG','SEEG','SEEG','SEEG','SEEG','SEEG','SEEG','SEEG','SEEG','SEEG','SEEG','SEEG'};
SBJ_vars.ch_lab.ref_type   = {'BP','BP','BP','BP','BP','BP','BP','BP','BP','BP','BP','BP','BP','BP','BP','BP'};
SBJ_vars.ch_lab.ROI        = {'all'};
SBJ_vars.ch_lab.eeg_ROI    = {};

%SBJ_vars.ch_lab.prefix = 'POL ';    % before every channel except 'EDF Annotations'
%SBJ_vars.ch_lab.suffix = '-Ref';    % after every channel except 'EDF Annotations'
SBJ_vars.ch_lab.mislabel = {{'RAM2-1','RSM2'},{'2MEG1','SMG1'},{'2MEG2','SMG2'},{'2MEG3','SMG3'},...% remove # in label of 2nd MEG probe
                            {'2MEG4','SMG4'}};%,{'2MEG5','SMG5'},{'2MEG6','SMG6'},{'2MEG7','SMG7'},...
%                            {'2MEG8','SMG8'},{'2MEG9','SMG9'},{'2MEG10','SMG10'}};

SBJ_vars.ref_exclude = {}; %exclude from the CAR
SBJ_vars.ch_lab.bad = {...
    'RTH4','RTH5','2MEG5','2MEG6','2MEG7','2MEG8','RAC8','RAC9','RAC10',...% spiking 
    'RAM10','ROF10',...% secondary spiking
    'MEG10','2MEG9','2MEG10','LAC9','LAC10','RPC10','RSM10',...% out of brain
    'REF','EKG','E','DC02','DC03','DC04',...% not real data
    'Z','Gnd','EKG'...% not real data
    };
SBJ_vars.ch_lab.eeg = {'FZ','CZ','OZ','C3','C4'};
SBJ_vars.ch_lab.eog = {'LUC','LLC','RUC','RLC'};
SBJ_vars.ch_lab.photod = {'DC01'};

%--------------------------------------
% Line Noise Parameters
%--------------------------------------
SBJ_vars.notch_freqs = [60 120 180 240 300];
SBJ_vars.bs_width    = 2;

%--------------------------------------
% Time Parameters
%--------------------------------------
% first trial is at 36 s (weird blip befor that at 34s
SBJ_vars.analysis_time = {{[24 1836]}};

%--------------------------------------
% Artifact Rejection Parameters
%--------------------------------------
SBJ_vars.artifact_params.std_limit_raw = 10;
SBJ_vars.artifact_params.hard_threshold_raw = 1000;

SBJ_vars.artifact_params.std_limit_diff = 10;
SBJ_vars.artifact_params.hard_threshold_diff = 60;

%--------------------------------------
% Trials to Reject
%--------------------------------------
% SBJ_vars.trial_reject_n = [];
