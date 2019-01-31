%% Photodiode Trace Cleaning Parameters: IR75
% Big drift down over course of task
% %   Tools --> linear + show equations:
% %   original attempt with one line: y = -0.000019*[1:numel(evnt.time{1})] + 18;
% %   segment 1:
% x1 = 1:43165; % flat
% y1 = zeros(size(x1));
% %   segment 2: 
% x2 = 43166:813620; % sloping down
% y2 = -0.00003693*x2 + 25; % old intercept without x1 values: 12.38;
% %   segment 3:
% x3 = 813621:1127900; % is flat ~5
% y3 = zeros(size(x3));
% %   segment 4: 
% x4 = 1127901:1849700; % sloping down
% y4 = -0.000040744*x4 + 53;% old intercept without x1:x3 values = 2.8949; 
% %   segment 5: 
% x5 = 1849701:numel(evnt.trial{1}); %jumps up and then is flat...
% y5 = zeros(size(x5));
% y = [y1 y2 y3 y4 y5];
% evnt.trial{1} = evnt.trial{1}-y;
% 
% % increase amplitude of photodide events to improve SNR
% evnt.trial{1}(evnt.trial{1}>15) = 60;

% Above manual didn't work, so trying high pass filtering
cfg = [];
cfg.hpfilter = 'yes';
cfg.hpfreq   = 0.1;
cfg.hpfiltord= 4;
evnt = ft_preprocessing(cfg,evnt);

% Mark trials to ignore e.g., interruptions
ignore_trials = []; 

% Set zero/baseline during a block
bsln_val = 0;

% Record epochs (in sec) with fluctuations that should be set to baseline
bsln_times = {...
    [0.0 49.0],... % initial wandering
    [1826.0 1900.0]...% end of task
    };
% Fix last datapoint
evnt.trial{1}(1,end) = bsln_val;

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

