function SU00_extract_evnt(SBJ)
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
b_ix = 1;   %block
eval(['run ' root_dir 'emodim/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run ' root_dir 'emodim/scripts/proc_vars/SU_nlx_proc_vars.m']);

if numel(SBJ_vars.raw_file)>1
    block_suffix = strcat('_',SBJ_vars.block_name{b_ix});
else
    block_suffix = SBJ_vars.block_name{b_ix};   % should just be ''
end

%% Read photodiode
inverted = 1;
photo = ft_read_neuralynx_interp({[SBJ_vars.dirs.SU 'photo/' SBJ_vars.ch_lab.photod{1} SBJ_vars.ch_lab.suffix '.ncs']});
photo.label = {'photo'};

% Preprocess
if inverted
    photo.trial{1} = photo.trial{1}*-1;
end
cfgpp = []; cfgpp.demean = 'yes';
photo = ft_preprocessing(cfgpp,photo);

% Cut to analysis time
if numel(SBJ_vars.analysis_time{b_ix})>1
    error('havent set up processing for multi block concat!');
end
cfgs = []; cfgs.latency = SBJ_vars.analysis_time{b_ix}{1};
photo = ft_selectdata(cfgs,photo);
photo.time{1} = photo.time{1}-SBJ_vars.analysis_time{b_ix}{1}(1);

% Resample if at very high srate
if proc_vars.photo_resample && photo.fsample~=proc_vars.evnt_resample_freq
    cfg_dsmp = [];
    cfg_dsmp.resamplefs = proc_vars.evnt_resample_freq;
    photo = ft_resampledata(cfg_dsmp,photo);
end    

%% Read, downsample, and save mic
mic = ft_read_neuralynx_interp({[SBJ_vars.dirs.SU 'mic/' SBJ_vars.ch_lab.mic{1} SBJ_vars.ch_lab.suffix '.ncs']});
mic.label = {'mic'};

% Cut to analysis time
mic = ft_selectdata(cfgs,mic);
mic.time{1} = mic.time{1}-SBJ_vars.analysis_time{b_ix}{1}(1);

% Downsample
if proc_vars.mic_resample
    cfg_dsmp = [];
    cfg_dsmp.resamplefs = photo.fsample;
    mic_dsmp = ft_resampledata(cfg_dsmp,mic);
end

%% Process Microphone data
mic_data = mic.trial{1};
%rescale to prevent clipping, add 0.05 fudge factor
mic_data_rescale = mic_data./(max(abs(mic_data))+0.05);

%% Concatenate photo and mic
% Check that time vectors are close enough
dif = photo.time{1}-mic_dsmp.time{1};
if max(dif)>0.000001
    error('time vectors not aligned, check that!');
elseif ~isequal(photo.time{1},mic_dsmp.time{1})
    warning(['Mic and photo time vectors still not equal, but differences is small: ' num2str(max(dif))]);
    mic_dsmp.time{1} = photo.time{1};
end

% Append
cfga = [];
cfga.keepsampleinfo = 'no';
evnt = ft_appenddata(cfga,photo,mic_dsmp);
evnt.fsample = photo.fsample;

%% Save data out
evnt_out_filename = strcat(SBJ_vars.dirs.import,SBJ,'_evnt',block_suffix,'.mat');
save(evnt_out_filename, '-v7.3', 'evnt');

% Save Microphone data separately as .wav for listening
mic_data_filename = strcat(SBJ_vars.dirs.import,SBJ,'_mic_recording',block_suffix,'.wav');
audiowrite(mic_data_filename,mic_data_rescale,evnt.fsample);

end