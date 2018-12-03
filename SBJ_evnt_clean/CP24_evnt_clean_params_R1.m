%% Photodiode Trace Cleaning Parameters: CP24 Run 1
[root_dir, ~] = fn_get_root_dir();

% Mark trials to ignore e.g., interruptions
ignore_trials = [277:330]; % all trials from B5, which are in R2

% Photodiode came loose, baseline shift and strength of events decreased
eval(['run ' root_dir 'emodim/scripts/SBJ_evnt_clean/CP24_R1_drift_corrections.m']);
evnt.trial{1} = new;

% Set zero/baseline during a block
bsln_val = -1.5;   % somehow we already demeaned this photodiode channel...
% Record epochs (in sec) with fluctuations that should be set to baseline
bsln_times = {...
    [307.0 308.837],...% weird blips around end of 1st correction
    [314.19 317.0]...% weird blip around end of correction 2
    };
% Record epochs (in sec) when baseline has shifted
bsln_shift_times = {};
% Record amount of shift in baseline for each epoch 
bsln_shift_val = [];
if length(bsln_shift_times)~=length(bsln_shift_val)
    error('Number of epochs and values for baseline shift periods do not match.');
end

% Record within trial corrections
stim_times = {};
stim_yval = [];
if length(stim_times)~=length(stim_yval)
    error('Number of epochs and values for stimulus correction periods do not match.');
end

