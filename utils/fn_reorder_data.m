function [ordered_data] = fn_reorder_data(data, labels)
%% Re-order an data struct to match order of labels provided
% If labels is empty, sort alphabetically
if isempty(labels)
    probes  = cell(size(data.label));
    lab_num = zeros(size(data.label));
    for l = 1:numel(data.label)
        % Find name of the probe
        probes{l} = data.label{l}(regexp(data.label{l},'\D'));
        % Find number of each channel
        lab_num(l) = str2num(data.label{l}(regexp(data.label{l},'\d')));
    end
    probes_sorted = sort(unique(probes));
    
    % Sort labels by probe then contact number
    labels = {};
    for p_ix = 1:numel(probes_sorted)
        cur_lab = data.label(strcmp(probes,probes_sorted{p_ix}));
        [~,sort_n_idx] = sort(lab_num(strcmp(probes,probes_sorted{p_ix})));
        labels = [labels; cur_lab(sort_n_idx)];
    end
end

ordered_data = data;
% Order them to match labels provided
[matches,order_idx] = ismember(labels,data.label);
if sum(matches)~=numel(labels)
    fprintf('Mismatched electrodes from input labels not found in data:\n');
    disp(label(matches==0));
    error('Mismatch between labels provided and data.label!');
else
    fields = fieldnames(data);
    for f = 1:numel(fields)
        if size(eval(['data.' fields{f}]),1) == numel(data.label)
            eval(['ordered_data.' fields{f} ' = data.' fields{f} '(order_idx(matches==1),:);']);
        elseif size(eval(['data.' fields{f}]),1) > numel(data.label)
            warning(['data field "' fields{f} '" has more elements than data.label, and will not be re-ordered!!!']);
        end 
    end 
    for t = 1:numel(data.trial)
        ordered_data.trial{t} = data.trial{t}(order_idx(matches==1),:);
    end
end

end