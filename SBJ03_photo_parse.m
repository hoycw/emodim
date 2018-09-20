function SBJ03_photo_parse(SBJ, block, rt_window, plot_it, save_it)
%
% SBJ [str] -- uniquely identifies the subject, e.g., 'IR54'
% block [int] -- index of which block of data should be analyzed
% rt_window [two element array] -- earliest and latest response onset following the onset of the word presentation
% plot_it [0/1] -- optional. plot_it = 1 to plot detected events
% save_it [0/1] -- whether to save output

n_hist_bins = 100; % Set this to a higher number for more amplitudes. More bins can mean more errors though.

%% File paths
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

% Import helper functions to Matlab path
helper_function_dir_name = [root_dir 'emodim/scripts/utils/'];
addpath(genpath(helper_function_dir_name));
SBJ_vars_cmd = ['run ' root_dir 'emodim/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
if numel(SBJ_vars.raw_file)>1
    block_suffix = strcat('_',SBJ_vars.block_name{block});
else
    block_suffix = SBJ_vars.block_name{block};   % should just be ''
end

evnt_filename = [SBJ_vars.dirs.preproc SBJ '_evnt_clean',block_suffix,'.mat'];
% log_filename = [SBJ_vars.dirs.raw SBJ '_strooptask_log',block_suffix,'.txt'];
%   create directory if needed
output_dir = [SBJ_vars.dirs.events];
output_filename = [output_dir SBJ '_trial_info_auto',block_suffix,'.mat'];

%% Determine event onset sample points
% Load preprocessing variables to get sampling rate
% hdr = ft_read_header([SBJ_vars.dirs.import SBJ '_1000hz.mat']);
% nrl_srate = SBJ_vars.resamp_freq;
% fprintf('WARNING!!! Assuming 32 kHz neural sampling rate!!!\n');
% nrl_srate = 32000;

% Load input data
fprintf('Loading %s\n',evnt_filename);
load(evnt_filename);
n_channels = hdr.n_channels;
n_samples = hdr.n_samples;
s_rate = hdr.sample_rate;
data_photo = evnt(photod_ix,:);
data_photo_orig = data_photo;
% data_mic = evnt(mic_ix,:);
% data_mic_orig = data_mic;

% Bring data down to zero
data_photo = data_photo - min(data_photo);

% Read photodiode data
fprintf('\tReading photodiode data\n');
min_event_length = 0.8 * s_rate;    %trial must be at least 0.8 sec (actually ~1.5s?)
% Decimate to around 1000Hz (don't worry about aliasing)
if s_rate > 1000
    decimate_v = floor(s_rate/1000);
    data_photo_d = data_photo(1:decimate_v:end);
    min_event_length = floor(min_event_length/decimate_v);
else
    fprintf('================================================\n');
    fprintf('WARNING: evnt sample rate is 1000 or less!!\n');
    data_photo_d = data_photo;
    decimate_v = 1;
end
[~, ~, data_shades] = read_photodiode(data_photo_d, min_event_length, 4);  %4 different shades (con, inc, neu, bsln)
clear data_photo;

% Diff to get edges which correspond to onsets and offsets
data_shades = [diff(data_shades) 0]; % Add a point because diff removes one
word_onsets = find(data_shades>0)'; % 1 to 2,3,4 is word onset. Transpose to make column vector
word_onsets = word_onsets*decimate_v; % convert back to ms
fprintf('\t\tFound %d trials in photodiode channel\n', length(word_onsets));

%% Read in log file
% Open file
fprintf('\tReading log file\n');
log_h = fopen(log_filename, 'r');

% Parse log file
file_contents = textscan(log_h, '%d %d %s %s %s %s %f %f %f', 'HeaderLines', 5, 'Delimiter', '\t', 'MultipleDelimsAsOne', 1);
trial_info.block_n = file_contents{1};
trial_info.trial_n = file_contents{2};
trial_info.word = file_contents{3};
trial_info.color = file_contents{4};
trial_info.trialtype = file_contents{5};
trial_info.blocktype = file_contents{6};
trial_info.response_time = file_contents{7};
trial_info.marker_time = file_contents{8};
trial_info.onset_time = file_contents{9};
fprintf('\t\tFound %d trials in log file\n', length(trial_info.trial_n));

% Remove trials to ignore
trial_info.block_n(ignore_trials) = [];
trial_info.trial_n(ignore_trials) = [];
trial_info.word(ignore_trials) = [];
trial_info.color(ignore_trials) = [];
trial_info.trialtype(ignore_trials) = [];
trial_info.blocktype(ignore_trials) = [];
trial_info.response_time(ignore_trials) = [];
trial_info.marker_time(ignore_trials) = [];
trial_info.onset_time(ignore_trials) = [];
fprintf('\t\tIgnoring %d trials\n', length(ignore_trials));

% If log and photodiode have different n_trials, plot and error out
if(length(trial_info.trial_n) ~= length(word_onsets))
    % Plot photodiode data
    plot_photo = data_photo_orig - min(data_photo_orig);
    plot_photo = plot_photo / (max(plot_photo)-min(plot_photo));
    plot_photo = plot_photo + 0.25;
    plot(plot_photo, 'k'); hold on;
    % Plot word onsets
    for word_n = 1:length(word_onsets)
        plot([word_onsets(word_n) word_onsets(word_n)],[1.30 1.40],'r','LineWidth',2);
        plot([word_onsets(word_n) word_onsets(word_n)],[-0.35 0.35],'r','LineWidth',2);
    end
    error('\nNumber of trials in log is different from number of trials found in event channel\n\n');
end

fprintf('\t\tFound %d word onsets\n', length(word_onsets));

% %% Get response onset times
% fprintf('\tObtaining response times\n');
% perblock_trial_n = 0;
% lastblock_n = trial_info.block_n(1);
% for trial_n = 1:length(trial_info.trial_n)
%     perblock_trial_n = perblock_trial_n+1;
%     if(trial_info.block_n(trial_n)>lastblock_n)
%         lastblock_n = lastblock_n+1;
%         perblock_trial_n = 1;
%     end
%     
%     % Load .wav files based on block and trial numbers
%     searchstr = ['B' num2str(trial_info.block_n(trial_n)) '-T' num2str(trial_info.trial_n(trial_n)) '_*'];
%     wavfile = dir(fullfile(wav_dir,searchstr));
%     try
%         [wavdata, wavsr] = audioread(fullfile(wav_dir,wavfile.name));
%     catch
%         % Try trial numbers resetting each block
%         searchstr = ['B' num2str(trial_info.block_n(trial_n)) '-T' ...
%             num2str(perblock_trial_n) '_*'];
%         wavfile = dir(fullfile(wav_dir,searchstr));
%         try
%             [wavdata, wavsr] = audioread(fullfile(wav_dir,wavfile.name));
%         catch
%             error('\n\nError reading\n\t[%s]\n\t[%s]\n', wav_dir, searchstr);
%         end
%     end
%     
%     % Filter the mic data in human voice range
%     voice_range = [300 2000]; % Hz
%     data_mic = X_eegfilt(wavdata',wavsr,voice_range(1),[]);
%     if(voice_range(2) < (wavsr/2))
%         data_mic = X_eegfilt(data_mic,wavsr,[],voice_range(2));
%     end
%     
%     % Calculate the moving variance
%     stdtime = 0.025; % In seconds
%     data_mic = movingstd(data_mic,ceil(stdtime*wavsr),'central');
%     
%     % Zero all values at or below baseline period (after word onset and before first possible response time
%     baseline_peak = max(data_mic(1:floor(rt_window(1)*wavsr)));
%     data_mic = max(data_mic, baseline_peak) - baseline_peak;
%     
%     % Zero all values below 1/10 of max peak
%     data_mic = max(data_mic, max(data_mic)/10) - max(data_mic)/10;
%     
%     % Dilate and then erode the mic amplitudes to remove noise
%     %  Dilate is kind of like doing the moving-maximum. Eroding then makes it "zero-phase".
%     %   This is also know as the morphological open in image processing.
%     %   These functions require the image processing toolbox
%     filt_window = ones(1,floor(wavsr*0.050)); % Window is about 50ms
%     data_mic = imdilate(data_mic,filt_window);
%     data_mic = imerode(data_mic,filt_window);
%     filt_window = ones(1,floor(wavsr*0.100)); % Window is about 100ms
%     data_mic = imerode(data_mic,filt_window);
%     data_mic = imdilate(data_mic,filt_window);
%     
%     % Get integral of mic data
%     datacum = cumtrapz(data_mic);
%     latest_rt = find(datacum > (max(datacum)/2),1);
%     if ~isempty(latest_rt)
%         data_mic = data_mic(1:latest_rt);
%     else
%         data_mic(:) = 0;
%     end
%     
%     % Get derivative of mic data
%     datadiff = diff(data_mic);
%     filt_window = ones(1,floor(wavsr*0.050)); % Window is about 50ms
%     data_mic = imdilate(datadiff,filt_window);
%     data_mic = imerode(data_mic,filt_window);
%     
%     % Get peak value in rt window
%     [~, peak_idx] = max(data_mic);
%     
%     % Find last time point before peak that has value of zero
%     last_zero_idx = find(data_mic(1:peak_idx)==0, 1, 'last');
%     if(isempty(last_zero_idx))
%         last_zero_idx = NaN;
%         fprintf('\t\tNo response peak found at Trial # %d\n', trial_n);
%     end
%     
%     % Subtract marker times
%     markerpnt = floor(trial_info.marker_time(trial_n)*wavsr);
%     trial_info.response_time(trial_n) = (last_zero_idx-markerpnt)/wavsr;
%     
%     % Check that rt is in bounds
%     if(trial_info.response_time(trial_n) > rt_window(2))
%         trial_info.response_time(trial_n) = NaN;
%     elseif(trial_info.response_time(trial_n) < rt_window(1))
%         trial_info.response_time(trial_n) = NaN;
%     end
%     
%     %plot(linspace(0,length(wavdata)/wavsr,length(wavdata)),wavdata); hold on;
%     %plot([trial_info.response_time(trial_n) trial_info.response_time(trial_n)],[min(wavdata) max(wavdata)],'r','LineWidth',3);
%     %pause;
%     %close all;
% end
% 
% fprintf('\t\tMean response time: %1.2f seconds from word onset\n', nanmean(trial_info.response_time(trial_n))); % Ignore NaNs
% 
% 
%% Put all the information into correct structures
trial_info.word_time = word_onsets;
% trial_info.resp_onset = word_onsets + floor(trial_info.response_time*s_rate);

% Make order of condition types [congruent neutral incongruent]
%trial_info.condition_types = unique(trial_info.condition); % This is more general, but since we know what we want:
trial_types = {'con','neu','inc'};
block_types = {'same','minc','mcon'};
k = 0;
for i = 1:length(trial_types)
    for j = 1:length(block_types)
        k = k+1;
        trial_info.condition_types{k} = [block_types{j} '-' trial_types{i}];
    end
end
for trial_n = 1:length(trial_info.trial_n)
    for type_n = 1:length(trial_info.condition_types)
        if(strcmp([trial_info.blocktype{trial_n} '-' trial_info.trialtype{trial_n}], trial_info.condition_types{type_n}))
            trial_info.condition_n(trial_n) = type_n;
        end
    end
end

%% Convert samples from event channel sample rate to ecog channel sample rate
% s_rate_ratio = nrl_srate / hdr.sample_rate;
trial_info.word_time = trial_info.word_time/hdr.sample_rate;
% trial_info.resp_onset = round(trial_info.resp_onset*s_rate_ratio);
trial_info.photo_sample_rate = hdr.sample_rate;

%% Output reaction times for trial types
fprintf('\tMean reaction times per condition:\n');
trial_info.rt_window = rt_window;
trial_info.ignore_trials = ignore_trials;
trial_info.SBJ = SBJ;
trial_info.rt_type = 'auto';
for type_n = 1:length(trial_info.condition_types)
    mean_rt = nanmean(trial_info.response_time(trial_info.condition_n==type_n));
    fprintf('\t\t%s: %1.3f seconds\n', trial_info.condition_types{type_n}, mean_rt);
end

%% Save results
if save_it
    save(output_filename, 'trial_info');
end

%% Plot results
if(plot_it ~= 0)
    figure;%('Position', [100 100 1200 800]);
    hold on;
    
    % Plot microphone data
    
%     plot_mic = data_mic_orig - mean(data_mic_orig);
%     plot_mic = min(plot_mic,5*std(data_mic_orig));
%     plot_mic = max(plot_mic,-5*std(data_mic_orig));
%     plot_mic = plot_mic / (max(plot_mic)-min(plot_mic));
%     plot(linspace(0, hdr.length_in_seconds, n_samples), plot_mic, 'Color', [0.8 0.8 0.5]); hold on;
%     plot([0 hdr.length_in_seconds],[0 0],'k');
    
    % Plot photodiode data
    plot_photo = data_photo_orig - min(data_photo_orig);
    plot_photo = plot_photo / (max(plot_photo)-min(plot_photo));
    plot_photo = plot_photo + 0.25;
    plot(linspace(0, hdr.length_in_seconds, n_samples), plot_photo, 'Color', [0.5 0.8 0.8]);
    plot([0 hdr.length_in_seconds],[0.25 0.25],'k');
    
    % Plot word onsets
    for word_n = 1:length(trial_info.word_time)
        plot([trial_info.word_time(word_n) trial_info.word_time(word_n)]/s_rate,[1.30 1.40],'b','LineWidth',2);
        plot([trial_info.word_time(word_n) trial_info.word_time(word_n)]/s_rate,[-0.35 0.35],'b','LineWidth',2);
    end
    
%     % Plot resp onsets
%     for resp_n = 1:length(trial_info.resp_onset)
%         plot([trial_info.resp_onset(resp_n)/s_rate_ratio trial_info.resp_onset(resp_n)/s_rate_ratio]/s_rate,[1.35 1.45],'g','LineWidth',2);
%         plot([trial_info.resp_onset(resp_n)/s_rate_ratio trial_info.resp_onset(resp_n)/s_rate_ratio]/s_rate,[-0.30 0.30],'g','LineWidth',2);
%     end
    
end

end

