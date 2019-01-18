%% Photodiode Trace Cleaning Parameters: IR61

% Mark trials to ignore e.g., interruptions
%   B3T59 interrupt, taking out B3T58 and on
%   B1 + B2 + B3T58 = 72 + 66 + 58 = 196:215
ignore_trials = [196:215]; 

% Set zero/baseline during a block
bsln_val = 900;   % somehow we already demeaned this photodiode channel...

% Record epochs (in sec) with fluctuations that should be set to baseline
%   nurse interruption @ 887s, removing video right before that too (though end of B3)
%   strange drifts in baseline in between blocks and indivdiual trials
% Also, strange extra event at beginning of 3 individually run blocks (B3 restart, B4, B5)
bsln_times = {...
    [880.0 1130.0],...% B3 interruption + between B3-B4 drift + initial B3 restart @ 1124.992
    [1450.0 1485.0],...% bw block drift + initial B4 @ 1483.02
    [1809.588 1862.0]...% bw block shift + initial B5 @ 1860.159
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

