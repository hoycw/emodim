function SBJ04_photo_parse(SBJ, block, plot_it, save_it)
% INPUTS:
%   SBJ [str] - uniquely identifies the subject, e.g., 'IR54'
%   block [int] - index of which block of data should be analyzed
%   plot_it [0/1] - optional. plot_it = 1 to plot detected events
%   save_it [0/1] - whether to save output

%% File paths
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
elseif exist('/Users/lapate/','dir');root_dir = '/Users/lapate/knight/';ft_dir = '/Users/lapate/knight/hoycw/Apps/fieldtrip/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'emodim/scripts/utils/']);

% Set up SBJ and directory info
SBJ_vars_cmd = ['run ' root_dir 'emodim/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
if numel(SBJ_vars.raw_file)>1
    block_suffix = strcat('_',SBJ_vars.block_name{block});
else
    block_suffix = SBJ_vars.block_name{block};   % should just be ''
end

evnt_filename = [SBJ_vars.dirs.preproc SBJ '_evnt_clean',block_suffix,'.mat'];
output_filename = [SBJ_vars.dirs.events SBJ '_trial_info_auto',block_suffix,'.mat'];

%% Determine event onset sample points
% Load input data
fprintf('Loading %s\n',evnt_filename);
load(evnt_filename);
n_channels = numel(evnt.label);
n_samples = size(evnt.trial{1},2);
s_rate = evnt.fsample;
data_photo = evnt.trial{1};
data_photo_orig = data_photo;

% Bring data down to zero
data_photo = data_photo - min(data_photo);

% Read photodiode data
fprintf('\tReading photodiode data\n');
min_event_length = 0.8 * s_rate;    %trial must be at least 0.8 sec (actually ~1.5s?)
[evnt_on, evnt_off, data_shades] = read_photodiode(data_photo, min_event_length, 2);  %2 different shades (bsln, evnt)
if save_it
    fig_fname = [SBJ_vars.dirs.events SBJ '_photo_segmentation.fig'];
    saveas(gcf,fig_fname);
end
clear data_photo;

% Diff to get edges which correspond to onsets and offsets
data_shades = [diff(data_shades) 0]; % Add a point because diff removes one
video_onsets = find(data_shades>0)'; % 1 to 2 is video onset. Transpose to make column vector
fprintf('\t\tFound %d trials in photodiode channel\n', length(video_onsets));

% Add in first video with no photodiode
if SBJ_vars.restart{block}
    %   `ffmpeg -i 0008.mp4` says duration = 2.1 s; actual photodiode in IR78 said 2.167s
    first_len = diff(video_onsets(1:2));
    if (first_len > 2.200*evnt.fsample) || (first_len < 2.050*evnt.fsample) % if the first video is NOT there
        video_onsets = [video_onsets(1)-2.1*evnt.fsample; video_onsets];
        fprintf('\t\tAdded back missing first video to make total videos found in photodiode: %d\n',length(video_onsets));
    end
end

% Plot photodiode event durations to check consistency
if plot_it
    figure;
    dur = evnt_off{2}-evnt_on{2};
    plot(dur);
    ylabel('Photodiode Durations');
    xlabel('Trial');
    title(['[min, mean, max] = [' num2str(min(dur)) ',' num2str(mean(dur)) ',' num2str(max(dur)) ']']);
    if any(dur<60)
        warning(['WARNING!!! ' num2str(sum(dur<60)) ' photodiode events are less than 60 ms!']);
        disp(dur(dur<60));
    end
    if save_it
        fig_fname = [SBJ_vars.dirs.events SBJ '_photo_durations.png'];
        saveas(gcf,fig_fname);
    end
end

%% Read in log file
% Open file
fprintf('\tReading log file\n');
log_h = fopen([SBJ_vars.dirs.events SBJ '_eventInfo' block_suffix '.txt'], 'r');

% Parse log file
file_contents = textscan(log_h, '%f %d', 'Delimiter', ',', 'MultipleDelimsAsOne', 1);
trial_info.video_id = file_contents{2};
trial_info.log_onset_time = file_contents{1};
fprintf('\t\tFound %d trials in log file\n', length(trial_info.video_id));

% Remove trials to ignore
trial_info.video_id(ignore_trials) = [];
trial_info.log_onset_time(ignore_trials) = [];
trial_info.ignore_trials = ignore_trials;
fprintf('\t\tIgnoring %d trials\n', length(ignore_trials));

% For some SBJs, use log time for 1st video to be consistent with other timing
if ismember(SBJ,{'IR71','IR74', 'IR75', 'IR84'})
    video_onsets(1) = video_onsets(2)-diff(trial_info.log_onset_time(1:2))*evnt.fsample;
end
% For IR78, toss ignore_trials from video_onsets here rather than zeroing
% to baseline in photodiode signal because it's messy and read_photodiode fails
if strcmp(SBJ,'IR78')
    video_onsets(ignore_trials) = [];
    fprintf('IR78: Tossing %i trials from photod video_onsets\n',numel(ignore_trials));
end

% If log and photodiode have different n_trials, plot and error out
if (length(trial_info.video_id) ~= length(video_onsets))
    % Plot photodiode data
    plot_photo = data_photo_orig - min(data_photo_orig);
    plot_photo = plot_photo / (max(plot_photo)-min(plot_photo));
    plot_photo = plot_photo + 0.25;
    figure; hold on;
    plot(plot_photo, 'k');
    % Plot video onsets
    for video_n = 1:length(video_onsets)
        plot([video_onsets(video_n) video_onsets(video_n)],[1.30 1.40],'r','LineWidth',2);
        plot([video_onsets(video_n) video_onsets(video_n)],[-0.35 0.35],'r','LineWidth',2);
    end
    error('\nNumber of trials in log is different from number of trials found in event channel\n\n');
end

% Compare onset differences between photodiode and log times
log_times = trial_info.log_onset_time-trial_info.log_onset_time(1);
video_times = (video_onsets-video_onsets(1))/evnt.fsample;
donsets = log_times-video_times;
dphoto = diff(video_onsets/evnt.fsample);
dlog   = diff(trial_info.log_onset_time);
ddif   = dlog-dphoto;
fprintf('\tMax difference in photodiode - log event onsets = %f\n',max(abs(donsets)));
fprintf('\tMax difference in photodiode - log event durations = %f\n',max(abs(ddif)));

trial_info.video_onsets = video_onsets;

%% Save results
if save_it
    save(output_filename,'-v7.3', 'trial_info');
end

%% Plot results
if(plot_it ~= 0)
    % Plot event onsets
    figure;%('Position', [100 100 1200 800]);
    hold on;
    
    % Plot photodiode data
    plot_photo = data_photo_orig - min(data_photo_orig);
    plot_photo = plot_photo / (max(plot_photo)-min(plot_photo));
    plot_photo = plot_photo + 0.25;
    plot(linspace(0, n_samples/s_rate, n_samples), plot_photo, 'Color', [0.5 0.8 0.8]);
    plot([0 n_samples/s_rate],[0.25 0.25],'k');
    
    % Plot video onsets
    for video_n = 1:length(trial_info.video_onsets)
        plot([trial_info.video_onsets(video_n) trial_info.video_onsets(video_n)]/s_rate,[1.30 1.40],'b','LineWidth',2);
        plot([trial_info.video_onsets(video_n) trial_info.video_onsets(video_n)]/s_rate,[-0.35 0.35],'b','LineWidth',2);
    end
    
    if save_it
        fig_fname = [SBJ_vars.dirs.events SBJ '_events.fig'];
        saveas(gcf,fig_fname);
    end
    
    % Plot difference between log and photodiode
    figure; hold on;
    subplot(2,1,1);
    plot(donsets);
    ylabel('Log-Photo Onsets');
    xlabel('Video');
    title(['max(abs(donsets)) = ' num2str(max(abs(donsets)))]);
    
    subplot(2,1,2);
    plot(ddif);
    ylabel('Log-Photo Durations');
    xlabel('Video');
    title(['max(abs(diff(durations))) = ' num2str(max(abs(ddif)))]);

    if save_it
        fig_fname = [SBJ_vars.dirs.events SBJ '_log_photo_QA.png'];
        saveas(gcf,fig_fname);
    end
end

end

