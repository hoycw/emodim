%% Photodiode Trace Cleaning Parameters: IR67

% Mark trials to ignore e.g., interruptions
ignore_trials = []; 

% Set zero/baseline during a block
bsln_val = -10;

% Record epochs (in sec) with fluctuations that should be set to baseline
bsln_times = {...
    [0.0 24.0],... % initial wandering
    [25.0 34.0],...% spike that may or may not be an event...
    [1700.0 1980.0]...% end of task
    };
% Fix last datapoint
evnt.trial{1}(1,end) = bsln_val;
% first several blocks have tiny signal amplitude, make all big size of last blocks
evnt.trial{1}(1,evnt.trial{1}(1,:)>10) = 230;

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

