video_onsets_orig = trial_info.video_onsets;
dphoto = diff(trial_info.video_onsets);

log_times = trial_info.log_onset_time;
n_gone = numel(log_times)-numel(trial_info.video_onsets);
if n_gone == 0
    n_gone = 1;
end
max_dif = zeros([n_gone 1]);
min_dif = zeros([n_gone 1]);
avg_dif = zeros([n_gone 1]);
for i = 1:n_gone
    dlog = diff(log_times(i:i+length(trial_info.video_onsets)-1));
    ddif = dlog-dphoto*evnt.fsample;
    min_dif(i) = min(ddif);
    max_dif(i) = max(ddif);
    avg_dif(i) = mean(ddif);
end

log_times = trial_info.log_onset_time;
log_times = log_times-log_times(1);
video_times = (video_onsets_orig-video_onsets_orig(1))/evnt.fsample;

log_match = zeros([numel(video_times) 1]);
for pv = 1:numel(video_times)
    diff_vec = log_times-video_times(pv);
    [~,log_match(pv)] = min(abs(diff_vec));
end

figure;
plot(log_times,'k');
hold on;
scatter(log_match,log_times(log_match),'r');