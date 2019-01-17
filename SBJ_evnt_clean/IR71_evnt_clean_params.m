%% Photodiode Trace Cleaning Parameters: IR71

% Mark trials to ignore e.g., interruptions
%   58:63 in photodiode events, but +1 for missing 1st video
ignore_trials = [59:64]; 

% Set zero/baseline during a block
bsln_val = 2830;   % somehow we already demeaned this photodiode channel...

% Record epochs (in sec) with fluctuations that should be set to baseline
%   beginning until ~1s before first photodiode signal
%   274-304 to remove nurse interruption at end of B1
%   1056-1066s pc for baseline crap
bsln_times = {...
    [274.0 304.0],... % nurse interruption end of B1
    [990.0 1000.0],... % baseline crap
    [1359.0 1360.0]... % baseline crap
    };
%    [0.0 9.0],...  % initial offset

    
% Record epochs (in sec) when baseline has shifted
bsln_shift_times = {};
% Record amount of shift in baseline for each epoch 
bsln_shift_val = [];
if length(bsln_shift_times)~=length(bsln_shift_val)
    error('Number of epochs and values for baseline shift periods do not match.');
end

% Record within trial corrections
%   events get small in middle and read_photo can't detect
%   all events are above 9500 after above corrections, so make them all big
stim_times = {};
stim_yval = [];
if length(stim_times)~=length(stim_yval)
    error('Number of epochs and values for stimulus correction periods do not match.');
end
% Corrected manually running this code after normal correction code:
% stim_epochs = find(evnt.trial{1}>9500);
% evnt.trial{1}(stim_epochs) = 70000;

