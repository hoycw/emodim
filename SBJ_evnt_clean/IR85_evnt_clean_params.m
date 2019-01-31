%% Photodiode Trace Cleaning Parameters: IR85

% Mark trials to ignore e.g., interruptions
ignore_trials = [1 2]; 

% Set zero/baseline during a block
bsln_val = 0;

% Record epochs (in sec) with fluctuations that should be set to baseline
bsln_times = {...
    [0.0 60.0],... % initial wandering
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

