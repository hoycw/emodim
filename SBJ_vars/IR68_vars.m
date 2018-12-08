%% IR68 Processing Variables
[root_dir, app_dir] = fn_get_root_dir();
if isempty(strfind(path,'fieldtrip'))
    addpath([app_dir 'fieldtrip/']);
    ft_defaults
end

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ        = 'IR68';
SBJ_vars.raw_file   = {'2018012413_0001.besa'};
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
SBJ_vars.recon.elec_mni_s = [];%SBJ_vars.dirs.recon 'Electrodes/' SBJ_vars.SBJ '_elec_mni_s.mat'];
SBJ_vars.recon.fs_T1      = [SBJ_vars.dirs.recon 'Scans/' SBJ_vars.SBJ '_fs_preop_T1.mgz'];
SBJ_vars.recon.fs_DK      = [SBJ_vars.dirs.recon 'Scans/' SBJ_vars.SBJ '_fs_preop_aparc+aseg.mgz'];
SBJ_vars.recon.fs_Dx      = [SBJ_vars.dirs.recon 'Scans/' SBJ_vars.SBJ '_fs_preop_aparc.a2009s+aseg.mgz'];

%--------------------------------------
% Channel Selection
%--------------------------------------
SBJ_vars.ch_lab.probes     = {'LAM','LHH','LTH','AIN','MIN','PIN','LOF','LAC','LPC'};
SBJ_vars.ch_lab.probe_type = {'seeg','seeg','seeg','seeg','seeg','seeg','seeg','seeg','seeg'};
SBJ_vars.ch_lab.ref_type   = {'BP','BP','BP','BP','BP','BP','BP','BP','BP'};
SBJ_vars.ch_lab.ROI        = {'all'};
SBJ_vars.ch_lab.eeg_ROI    = {};

%SBJ_vars.ch_lab.prefix = 'POL ';    % before every channel except 'EDF Annotations'
%SBJ_vars.ch_lab.suffix = '-Ref';    % after every channel except 'EDF Annotations'
%SBJ_vars.ch_lab.mislabel = {{'RLT12','FPG12'},{'IH;L8','IHL8'}};

SBJ_vars.ref_exclude = {}; %exclude from the CAR
SBJ_vars.ch_lab.bad = {...
    '*-1',...%bad duplicate channels
    'LAM7','LAM8','LAM9','LAM10','LHH8','LHH9','LHH10','LTH8','LTH9','LTH10',...%epileptic
    'AIN8','AIN9','AIN10',...%out of brain
    'AIN5',...% loose
    'EKG',...% not using right now
    'LPL*','GRND','XREF','DC02','DC03','DC04'...% not real data
    };
SBJ_vars.ch_lab.eeg = {'FZ','CZ','OZ','C3','C4'};
% SBJ_vars.ch_lab.CZ_lap_ref = {};
SBJ_vars.ch_lab.eog = {'LUE','LLE','RUE','RLE'};
SBJ_vars.ch_lab.photod = {'DC01'};
%SBJ_vars.ch_lab.mic    = {};

%--------------------------------------
% Line Noise Parameters
%--------------------------------------
SBJ_vars.notch_freqs = [60 120 180 240 300];
SBJ_vars.bs_width    = 2;

%--------------------------------------
% Time Parameters
%--------------------------------------
% behavioral log shows 1st 2 videos repeated, task was restarted
% task was interrupted after first 2 trials in block 5, so cut those 2 off
SBJ_vars.analysis_time = {{[25 1330]}};

%--------------------------------------
% Artifact Rejection Parameters
%--------------------------------------
SBJ_vars.artifact_params.std_limit_raw = 7;
SBJ_vars.artifact_params.hard_threshold_raw = 700;

SBJ_vars.artifact_params.std_limit_diff = 7;
SBJ_vars.artifact_params.hard_threshold_diff = 50;

%--------------------------------------
% Trials to Reject
%--------------------------------------
% SBJ_vars.trial_reject_n = [];
