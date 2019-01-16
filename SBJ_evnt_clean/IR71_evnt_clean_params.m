%% Photodiode Trace Cleaning Parameters: IR71

% Mark trials to ignore e.g., interruptions
ignore_trials = []; 

% Set zero/baseline during a block
bsln_val = 0;   % somehow we already demeaned this photodiode channel...

% Record epochs (in sec) with fluctuations that should be set to baseline
%   beginning until ~1s before first photodiode signal
%   274-304 to remove nurse interruption at end of B1
%   1056-1066s pc for baseline crap
bsln_times = {...
    [0.0 9.0],...  % initial offset
    [274.0 304.0],... % nurse interruption end of B1
    [990.0 1000.0]... % baseline crap
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

