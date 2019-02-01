%% Photodiode Trace Cleaning Parameters: IR84
% very slow drift:
cfg = [];
cfg.hpfilter = 'yes';
cfg.hpfreq   = 0.1;
cfg.hpfiltord= 4;
evnt = ft_preprocessing(cfg,evnt);

% Mark trials to ignore e.g., interruptions
%   screen saver came on and messed up timing of 297 (2091.mp4) and also last trial
%   tossing all 5 of those to be certain
ignore_trials = [297:301]; 

% Set zero/baseline during a block
bsln_val = -5;

% Record epochs (in sec) with fluctuations that should be set to baseline
bsln_times = {...
    [0.0 70.0],... % initial wandering
    [132.0 136.0],...% spike after baseline shift
    [1536.0 1685.0]...% last 5 trials + end of task
    };
% Fix last datapoint
evnt.trial{1}(1,end) = bsln_val;

% Record epochs (in sec) when baseline has shifted
bsln_shift_times = {};% without hp filter: [70.0 134.0]};
% Record amount of shift in baseline for each epoch 
bsln_shift_val = [];% without hp filter: [125];
if length(bsln_shift_times)~=length(bsln_shift_val)
    error('Number of epochs and values for baseline shift periods do not match.');
end

% Record within trial corrections
stim_times = {};
stim_yval = [];
if length(stim_times)~=length(stim_yval)
    error('Number of epochs and values for stimulus correction periods do not match.');
end

