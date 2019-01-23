%% HFA spike detection
hfa_thresh = 8;
hfa_detect = [];
for c = 1:numel(hfa_trl.label)
    c_thresh = std(squeeze(hfa_data.trial{1}(c,:)))*hfa_thresh;
    for t = 1:numel(hfa_trl.trial)
        if sum(hfa_trl.trial{t}(c,:)>c_thresh) > 0.1*hfa_data.fsample
            hfa_detect = [hfa_detect(:); t];clear
        end
    end
end
hfa_detect = unique(hfa_detect);
hfa_dif_epochs = hfa_dif_trl.sampleinfo(hfa_detect,:);

% hfa_dif_epochs = [];
% for c = 1:numel(hfa_dif_trl.label)
%     
% end

%% Variance 
nch = numel(hfa_trl.label);
ntr = numel(hfa_trl.trial);
var_dist = zeros([nch ntr]);
for c = 1:nch
    for t = 1:ntr
        var_dist(c,t) = nanstd(hfa_trl.trial{t}(c,:));
    end
end

% % Plot raw var values
% figure;
% plot(var_dist')

% plot all var dist
n_bins = 1000;
bins = linspace(min(min(var_dist)),max(max(var_dist)),n_bins);
log_bins = linspace(0,max(log(var_dist(:))),100);
var_hist = zeros([nch, numel(bins)-1]);
var_hist_log = zeros([nch, numel(log_bins)-1]);
for c = 1:nch
    [var_hist(c,:),~] = histcounts(var_dist(c,:),bins);
    [var_hist_log(c,:),~] = histcounts(log(var_dist(c,:)),log_bins);
end
figure;
subplot(2,2,1);
plot(bins(2:end),var_hist);
subplot(2,2,3);
plot(log_bins(2:end),var_hist_log);
subplot(2,2,2);hold on;
for c = 1:nch
    scatter(var_hist(c,:),ones(size(var_hist(c,:)))*c);
end
subplot(2,2,4);hold on;
for c = 1:nch
    scatter(var_hist_log(c,:),ones(size(var_hist_log(c,:)))*c);
end

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

