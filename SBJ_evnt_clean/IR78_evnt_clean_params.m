%% Photodiode Trace Cleaning Parameters: IR78
% This photodiode trace requires no fixing! Yay!
%   Actually, it has some weird baseline shifts, but the algorithm handles them fine so onsets are good

% Mark trials to ignore e.g., interruptions
ignore_trials = []; 
% first two trials in logs were before restart, and a repeated
% last two trials were first ones in last block, which was interrupted immediately

% Set zero/baseline during a block
bsln_val = 0;   % somehow we already demeaned this photodiode channel...

% Record epochs (in sec) with fluctuations that should be set to baseline
bsln_times = {...
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

