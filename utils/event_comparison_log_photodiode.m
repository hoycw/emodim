dphoto = diff(video_onsets);

log_times = trial_info.log_onset_time;
n_gone = numel(log_times)-numel(video_onsets);
max_dif = zeros([n_gone 1]);
min_dif = zeros([n_gone 1]);
avg_dif = zeros([n_gone 1]);
for i = 1:n_gone
    dlog = diff(log_times(i:i+length(video_onsets)-1));
    ddif = dlog-dphoto;
    min_dif(i) = min(ddif);
    max_dif(i) = max(ddif);
    avg_dif(i) = mean(ddif);
end

log_times = trial_info.log_onset_time;
log_times = log_times-log_times(end);
video_onsets = (video_onsets_orig-video_onsets_orig(end))/1000;

log_match = zeros([numel(video_onsets) 1]);
for pv = 1:numel(video_onsets)
    diff_vec = log_times-video_onsets(pv);
    [~,log_match(pv)] = min(abs(diff_vec));
end

figure;
plot(log_times,'k');
hold on;
scatter(log_match,log_times(log_match),'r');