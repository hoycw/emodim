%% Photodiode Trace Cleaning Parameters: IR78
% This photodiode trace requires no fixing! Yay!
%   Actually, it has some weird baseline shifts, but the algorithm handles them fine so onsets are good

% Mark trials to ignore e.g., interruptions
%   a priest came in to talk to the patient during the last block, from 24:05-24:35 in preclean time
%   1440:1475s in preclean = 1426:1451 in analysis_time, which is videos 1802-1992 in run5.csv
ignore_trials = [287:293]; 
% first two trials in logs were before restart, and a repeated
% last two trials were first ones in last block, which was interrupted immediately

% Set zero/baseline during a block
bsln_val = 0;   % somehow we already demeaned this photodiode channel...

% Record epochs (in sec) with fluctuations that should be set to baseline
bsln_times = {...
    };
% would add this to bsln_times, but it makes the read_photodiode fn break: [1426.0 1451.0]...% priest distraction
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

