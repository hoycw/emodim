%% Photodiode Trace Cleaning Parameters: CP24 Run 1
% This photodiode trace requires no fixing! Yay!

% Mark trials to ignore e.g., interruptions
ignore_trials = [277:330]; % all trials from B5, which are in R2

% Photodiode came loose, baseline shift and strength of events decreased
%   x1 = [256204, 307287];
%   y1 = [4.3955 87.912];
%   coeff  = polyfit(x1,y1, 1);
%   line1  = coeff(1)*[x1(1):x1(2)] + coeff(2);
%   x2 = [x1(2) 309015];
%   y2 = [y1(2) -3.5];
%   coeff2 = polyfit(x2,y2,1);
%   line2  = coeff2(1)*[x2(1):x2(2)] + coeff2(2);
%   
%   new = evnt.trial{1};
%   new(x1(1):x1(2)) = new(x1(1):x1(2))-line1;
%   new(x2(1):x2(2)) = new(x2(1):x2(2))-line2;
%   
%   x3 = [313245 315690];
%   y3 = [-3.4 81.32];
%   coeff3  = polyfit(x3,y3, 1);
%   line3  = coeff3(1)*[x3(1):x3(2)] + coeff3(2);
%   x4 = [x3(2) 316165];
%   y4 = [y3(2) -2.198];
%   coeff4 = polyfit(x4,y4,1);
%   line4  = coeff4(1)*[x4(1):x4(2)] + coeff4(2);
%   
%   new(x3(1):x3(2)) = new(x3(1):x3(2))-line3;
%   new(x4(1):x4(2)) = new(x4(1):x4(2))-line4;
%   
%   evnt.trial{1} = new;
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

