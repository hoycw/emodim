function [trial_info_clean] = SBJ05_reject_artifacts_KLA_report(trials, trial_info, artifact_struct, plot_chans, report)
% trials is a fieldtrip data structure already segmented into trials
% artifact_struct contains the thresholds to reject trials
% plot_chans is a cell array of channels to be plotted
%   {} - empty array indicates skip plotting
%   'worst' - says plot channels that had > 5 bad trials
%   full array - list of channel names to plot
% report [struct] - list of 0/1/2 flags for plotting and reporting options,
%   0 = no report; 1 = concise report; 2 = verbose report
%   .hard_thresh: 1=number of trials; 2=print rejected trial_n
%   .std_thresh: 1=number of trials; 2=print rejected trial_n
%   .std_plot: 1=plot std distribution with cut off

% if(nargin < 3)
%   plot_chans = [];
% end

%% File paths
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'Apps/fieldtrip/'];
elseif exist('/Users/lapate/','dir');root_dir = '/Users/lapate/knight/';ft_dir = '/Users/lapate/knight/hoycw/Apps/fieldtrip/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
helper_function_dir_name = [root_dir 'PRJ_Stroop/scripts/_TOOLBOXES/'];
if(~exist(fullfile(helper_function_dir_name)))
  fprintf('\nERROR: Helper function directory does not exist.\n');
  return;
end

% Import helper functions to Matlab path
addpath(genpath(fullfile(helper_function_dir_name)));

%% Load input data
s_rate = trials.fsample;

plot_struct.plot_chans = [];
if isempty(plot_chans)
    plot_struct.plot_chan_names = {};
else
    if strcmp(plot_chans{1},'worst')  % send all channels, later only plot those with > 5 rejected trials
        plot_struct.plot_chan_names = trials.label;
        plot_struct.bad_trial_thresh = plot_chans{2};
    else
        plot_struct.plot_chan_names = plot_chans;
    end
    for ch_ix = 1:numel(plot_struct.plot_chan_names)
        plot_struct.plot_chans = [plot_struct.plot_chans strmatch(plot_struct.plot_chan_names{ch_ix},trials.label,'exact')];
    end
end

%% Convert to single data matrix
[trial_mat,~] = fn_format_trials_ft2KLA(trials);
n_epochs = size(trial_mat,2);

%% Artifact rejection
% Reject artifacts from raw data
fprintf('==============================================================================================\n');
fprintf('=================================== PHASE 1: ORIGINAL DATA ===================================\n');
fprintf('==============================================================================================\n');
ok_epochs_raw = artifact_reject(trial_mat, artifact_struct.std_limit_raw,...
    artifact_struct.hard_threshold_raw, s_rate, plot_struct, report, trial_info);
% Reject artifacts from diff data (looking for fast changes)
fprintf('==============================================================================================\n');
fprintf('================================== PHASE 2: DERIVATIVE DATA ==================================\n');
fprintf('==============================================================================================\n');
ok_epochs_diff = artifact_reject(diff(trial_mat,1,3), artifact_struct.std_limit_diff,...
    artifact_struct.hard_threshold_diff, s_rate, plot_struct, report, trial_info);
% Get intersection of the two arrays
ok_epochs = intersect(ok_epochs_raw,ok_epochs_diff);

%% Compile output in trial_info_clean
trial_info_clean = trial_info;
trial_info_clean.artifact_params = artifact_struct;

% Document bad trials
var_reject_raw_ix = setdiff(trial_info.trial_n,ok_epochs_raw);
var_reject_diff_ix = setdiff(trial_info.trial_n,ok_epochs_diff);
trial_info_clean.bad_trials.KLA_raw = var_reject_raw_ix;
trial_info_clean.bad_trials.KLA_diff = var_reject_diff_ix;
trial_info_clean.bad_trials.KLA_all = sort(unique([trial_info_clean.bad_trials.KLA_raw; trial_info_clean.bad_trials.KLA_diff]));

% Print results
fprintf('\tTOTAL %d OUT OF %d EPOCHS RETAINED\n', length(ok_epochs), n_epochs);

%% Artifact rejection function
function [ok_epochs] = artifact_reject(current_data, std_limit, hard_threshold, s_rate, plot_struct, report, trial_info)

n_channels = size(current_data,1);
n_epochs = size(current_data,2);

% Remove epochs that exceed a certain threshold
ok_epochs = 1:n_epochs;
threshold_epochs = squeeze(max(max(abs(current_data),[],3),[],1)) > hard_threshold;
if report.hard_thresh > 0
    fprintf('\t----> Hard Threshold rejected %i trials, trial_n=:\n',sum(threshold_epochs));
    if report.hard_thresh > 1
        disp(ok_epochs(threshold_epochs));
    end
end
ok_epochs(threshold_epochs) = [];

% Remove temporal mean
for channel_n = 1:n_channels
  for epoch_n = 1:n_epochs
    current_data(channel_n,epoch_n,:) = detrend(squeeze(current_data(channel_n,epoch_n,:)),'constant');
  end
end
% Keep removing epochs until all are below standard deviation threshold
%   Initialize variables for loop
n_bad_epochs_last = -1;
n_bad_epochs = 0;
exceeded_channels_last1 = -ones(n_channels,1);
exceeded_channels_last2 = zeros(n_channels,1);
%   Start checking EEG Standard Deviations
[eeg_std] = get_std_of_data(current_data(:,ok_epochs,:), [], s_rate);
std_last = eeg_std(:,1,:); clear eeg_std;
iter_n = 0;
while( n_bad_epochs > n_bad_epochs_last )
  iter_n = iter_n+1;
  n_bad_epochs_last = n_bad_epochs;
  % Get standard deviation of data (excluding bad epochs)
  [eeg_std] = get_std_of_data(current_data(:,ok_epochs,:), [], s_rate); clear ok_epochs; clear bad_epochs;
  std_data = eeg_std(:,1,:); clear eeg_std;
  % Replace standard deviation on channels whose number of bad epochs has not increased
  for channel_n = 1:n_channels
    if(exceeded_channels_last1(channel_n) <= exceeded_channels_last2(channel_n))
      % Number of bad epochs for channel has not increased, use old std
      std_data(channel_n,1,:) = std_last(channel_n,1,:);
    end
  end
  % Find bad epochs
  [ok_epochs, bad_epochs, exceeded_channels] = find_exceeded_std(current_data, std_data, std_limit, s_rate);
  n_bad_epochs = numel(bad_epochs);
  std_last = std_data;
  exceeded_channels_last2 = exceeded_channels_last1;
  exceeded_channels_last1 = exceeded_channels;
  
  fprintf('\tITERATION %2.0f\tOK: %3.0f\tREJECTED: %3.0f\tPCT REJ: %3.0f\n',...
      iter_n, numel(ok_epochs), n_bad_epochs, (n_bad_epochs/n_epochs)*100);
end
fprintf('\tOK: %3.0f\tREJECTED: %3.0f\tPCT REJ: %3.0f\n',...
    numel(ok_epochs), numel(bad_epochs), (n_bad_epochs/n_epochs)*100);
if report.std_thresh > 0
    fprintf('\t----> StD Threshold rejected %i trials, trial_n:\n',numel(bad_epochs));
    if report.std_thresh > 1
        disp(trial_info.trial_n(bad_epochs));
    end
end
 
if report.std_plot > 0
    figure; hold on;
    for ch_ix = 1:size(current_data,1)
        max_vals = squeeze(max(abs(current_data(ch_ix,:,:)),[],3));
        lim = max(std_data(ch_ix,:,:)*std_limit);
        bad = find(max_vals>lim);
        good = setdiff(1:size(current_data,2),bad);
        scatter(max_vals(good),ones([1 numel(good)])*ch_ix,'MarkerEdgeColor','b','Marker','*');
        scatter(max_vals(bad),ones([1 numel(bad)])*ch_ix,...
            'MarkerFaceColor','r','MarkerEdgeColor','r','Marker','o');%,'LineWidth',2);
        scatter(lim,ch_ix,'MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',2,'Marker','d');
    end
    ylim([0 size(current_data,1)+1]);
    ylabel('channel n');
    xlabel('Maximum STD per trial');
    title(strcat('Trial STDs exceeding ',num2str(std_limit)));
end

% Plot results
if(~isempty(plot_struct.plot_chans))
    if isfield(plot_struct,'bad_trial_thresh')
        plot_struct.plot_chans = find(exceeded_channels > plot_struct.bad_trial_thresh);
        if isempty(plot_struct.plot_chans)  % If nothing exceeds threshold, take top 5
            [~,sort_idx] = sort(exceeded_channels,'descend');
            plot_struct.plot_chans = sort_idx(1:5);
        end
        plot_struct.plot_chan_names = plot_struct.plot_chan_names(plot_struct.plot_chans);
        plot_results(current_data, std_data, std_limit, s_rate, ok_epochs, plot_struct, exceeded_channels);
    else
        plot_results(current_data, std_data, std_limit, s_rate, ok_epochs, plot_struct, exceeded_channels);
    end
end

%% Find epochs with data that exceeds std limit function
function [ok_epochs, bad_epochs, exceeded_channels] = find_exceeded_std(eeg_data, std_data, std_limit, s_rate)
% Find data that exceeds threshold

% Get the absolute value of the difference between each sample point and the ensemble mean
%   Ensemble mean will be the robust mean (using a tukey window on the sorted values for each time point)
%   This is less affected by outliers
n_channels = size(eeg_data,1);
n_epochs = size(eeg_data,2);
n_samples = size(eeg_data,3);
ensemble_mean = zeros(n_channels, 1, n_samples);
weight_window = tukeywin(n_epochs,0.40)'; % 60% is square windowed and 20% on each side is tapered
weight_window = weight_window ./ sum(weight_window);
for channel_n = 1:n_channels
  for sample_n = 1:n_samples
    sorted_data = sort(eeg_data(channel_n,:,sample_n));
    ensemble_mean(channel_n,1,sample_n) = weighted_mean(sorted_data,weight_window);
  end
end
% Smooth the ensemble mean by first doing a moving maximum and then a moving average
movmax_window = ones(floor(s_rate*0.050),1); % Window for moving max is about 50ms
smooth_length = (length(movmax_window)*2)/n_samples; % This is a proportion compared to data length. 2*movmax_window.
pad_length = length(movmax_window); % The number of samples to add to each side before smoothing to mitigate edge effects
for channel_n = 1:n_channels
  % Pad the data on each side;
  curr_data = padarray(squeeze(ensemble_mean(channel_n,1,:)),pad_length,'symmetric','both');
  curr_data = smooth(curr_data,smooth_length,'moving');
  curr_data = imdilate(curr_data,movmax_window);
  curr_data = smooth(curr_data,smooth_length,'moving');
  ensemble_mean(channel_n,1,:) = curr_data((pad_length+1):(end-pad_length));
end
abs_eeg_data = abs(eeg_data-repmat(ensemble_mean,[1 size(eeg_data,2) 1]));

exceeded_eeg_data = zeros(size(abs_eeg_data));

n_channels = size(eeg_data,1);
for channel_n = 1:n_channels
  for epoch_n = 1:n_epochs
    % For each channel and trial, does the absolute value exceed std*limit computed across trials within electrode
    exceeded_eeg_data(channel_n,epoch_n,:) = abs_eeg_data(channel_n,epoch_n,:) > (std_data(channel_n,:,:)*std_limit);
  end
end
all_epochs = max(max(exceeded_eeg_data,[],3),[],1);
n_bad_epochs = 0;
n_ok_epochs = 0;
ok_epochs = [];
bad_epochs = [];
for epoch_n = 1:n_epochs
  if(all_epochs(epoch_n)>0)
    n_bad_epochs=n_bad_epochs+1;
    bad_epochs(n_bad_epochs)=epoch_n;
  else
    n_ok_epochs=n_ok_epochs+1;
    ok_epochs(n_ok_epochs)=epoch_n;
  end
end

% Determine channels that exceeded threshold
exceeded_channels = zeros(n_channels,1);
for channel_n = 1:n_channels
  exceeded_channels(channel_n) = sum(max(exceeded_eeg_data(channel_n,:,:),[],3));
end


%% Get standard deviation function
function [eeg_std] = get_std_of_data(eeg_data, bad_channels, s_rate)

n_channels = size(eeg_data,1);
n_samples = size(eeg_data,3);

% Calculate the ensemble standard deviation
eeg_std = std(eeg_data,[],2);

% Smooth the standard deviation by first doing a moving maximum and then a moving average
movmax_window = ones(floor(s_rate*0.050),1); % Window for moving max is about 50ms
smooth_length = (length(movmax_window)*2)+1; % This is a proportion compared to data length. Twice as long as movmax_window.
pad_length = length(movmax_window); % The number of samples to add to each side before smoothing to mitigate edge effects
for channel_n = 1:n_channels
  % Pad the data on each side;
  curr_data = padarray(squeeze(eeg_std(channel_n,1,:)),pad_length,'symmetric','both');
  curr_data = smooth(curr_data,smooth_length,'moving');
  curr_data = imdilate(curr_data,movmax_window);
  curr_data = smooth(curr_data,smooth_length,'moving');
  eeg_std(channel_n,1,:) = curr_data((pad_length+1):(end-pad_length));
end

% Make standard deviation in 'bad_channels' equal a large number
if(~isempty(bad_channels))
  eeg_std(bad_channels,1,:) = 1000;
end

%% Plot results function
function plot_results(eeg_data, std_data, std_limit, s_rate, ok_epochs, plot_struct, exceeded_channels)

plot_chans = plot_struct.plot_chans;

figure;%('Position', [20 50 1200 900]);

% Scale the data by dividing by the largest standard deviation in each channel
max_amp = squeeze(max(std_data*std_limit,[],3));
min_amp = squeeze(min(std_data*-std_limit,[],3));
scale_amp = max(max_amp,min_amp);
eeg_data =  eeg_data ./ repmat(scale_amp,[1 size(eeg_data,2) size(eeg_data,3)]);
std_data =  std_data ./ repmat(scale_amp,[1 size(std_data,2) size(std_data,3)]);

%   Ensemble mean will be the robust mean (using a tukey window on the sorted values for each time point)
%   This is less affected by outliers
n_channels = size(eeg_data,1);
n_epochs = size(eeg_data,2);
n_samples = size(eeg_data,3);
ensemble_mean = zeros(n_channels, 1, n_samples);
weight_window = tukeywin(n_epochs,0.40)'; % 60% is square windowed and 20% on each side is tapered
weight_window = weight_window ./ sum(weight_window);
for channel_n = 1:n_channels
  for sample_n = 1:n_samples
    sorted_data = sort(eeg_data(channel_n,:,sample_n));
    ensemble_mean(channel_n,1,sample_n) = weighted_mean(sorted_data,weight_window);
  end
end
% Smooth the ensemble mean by first doing a moving maximum and then a moving average
movmax_window = ones(floor(s_rate*0.050),1); % Window for moving max is about 50ms
smooth_length = (length(movmax_window)*2)/n_samples; % This is a proportion compared to data length. 2*movmax_window.
pad_length = length(movmax_window); % The number of samples to add to each side before smoothing to mitigate edge effects
for channel_n = 1:n_channels
  % Pad the data on each side;
  curr_data = padarray(squeeze(ensemble_mean(channel_n,1,:)),pad_length,'symmetric','both');
  curr_data = smooth(curr_data,smooth_length,'moving');
  curr_data = imdilate(curr_data,movmax_window);
  curr_data = smooth(curr_data,smooth_length,'moving');
  ensemble_mean(channel_n,1,:) = curr_data((pad_length+1):(end-pad_length));
end
mean_eeg_data = ensemble_mean;

% Before
subplot(1,2,1); hold on;
set(gca, 'Position', [0.07 0.16 0.45 0.83]);
last_max = 0;
for chan_n = 1:numel(plot_chans)
  y_origin(chan_n) = last_max - min(min(min(std_data(plot_chans(chan_n),:,:)*-std_limit))) + 0.25;
  curr_std_data = std_data(plot_chans(chan_n),:,:)*std_limit;
  plot(squeeze(eeg_data(plot_chans(chan_n),:,:))'+y_origin(chan_n));
  plot(squeeze(mean_eeg_data(plot_chans(chan_n),:,:)+curr_std_data)+y_origin(chan_n), 'r','LineWidth', 2);
  plot(squeeze(mean_eeg_data(plot_chans(chan_n),:,:)-curr_std_data)+y_origin(chan_n),'r','LineWidth', 2);
  last_max = y_origin(chan_n) + max(max(max(std_data(plot_chans(chan_n),:,:)*std_limit)));
end
ylim([-0.2 last_max+0.2]);
xlim([0 size(eeg_data,3)]);
set(gca,'YTick', y_origin, 'YTickLabel', plot_struct.plot_chan_names);
set(gca,'XTick', []);
set(gca,'FontSize',5);

% After
subplot(1,2,2); hold on;
set(gca, 'Position', [0.53 0.16 0.45 0.83]);
for chan_n = 1:numel(plot_chans)
  plot(squeeze(eeg_data(plot_chans(chan_n),ok_epochs,:))'+y_origin(chan_n));
end
ylim([-0.2 last_max+0.2]);
xlim([0 size(eeg_data,3)]);
set(gca,'YTick', []);
set(gca,'XTick', []);

% Rejected epochs per channel
if length(exceeded_channels) > 1
  subplot(50,1,50); hold on;
  set(gca, 'Position', [0.02 0.03 0.96 0.12]);
  for b_ix = 1:numel(exceeded_channels)
    b = bar(b_ix,exceeded_channels(b_ix));
    if find(plot_chans==b_ix)
      set(b,'FaceColor','r');
    end
  end
  xlim([0 length(exceeded_channels)]);
  set(gca,'XTick', [1 5:5:(length(exceeded_channels)-1) length(exceeded_channels)]);
end

hold off;
pause;
% close all;


function y = weighted_mean(x,w,dim)
%WMEAN   Weighted Average or mean value.
%   For vectors, WMEAN(X,W) is the weighted mean value of the elements in X
%   using non-negative weights W. For matrices, WMEAN(X,W) is a row vector 
%   containing the weighted mean value of each column.  For N-D arrays, 
%   WMEAN(X,W) is the weighted mean value of the elements along the first 
%   non-singleton dimension of X.
%
%   Each element of X requires a corresponding weight, and hence the size 
%   of W must match that of X.
%
%   WMEAN(X,W,DIM) takes the weighted mean along the dimension DIM of X. 
%
%   Class support for inputs X and W:
%      float: double, single
%
%   Example:
%       x = rand(5,2);
%       w = rand(5,2);
%       wmean(x,w)
%
% Copyright (c) 2009, Adam Auton
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.


if nargin<2
    error('Not enough input arguments.');
end

% Check that dimensions of X match those of W.
if(~isequal(size(x), size(w)))
    error('Inputs x and w must be the same size.');
end

% Check that all of W are non-negative.
if (any(w(:)<0))
    error('All weights, W, must be non-negative.');
end

% Check that there is at least one non-zero weight.
if (all(w(:)==0))
    error('At least one weight must be non-zero.');
end

if nargin==2, 
  % Determine which dimension SUM will use
  dim = find(size(x)~=1, 1 );
  if isempty(dim), dim = 1; end
end

y = sum(w.*x,dim)./sum(w,dim);















































































