function [ordered_elec] = fn_reorder_elec(elec, labels)
%% Re-order an elec struct to match order of labels provided
% If labels is empty, sort alphabetically
if isempty(labels)
    probes  = cell(size(elec.label));
    lab_num = zeros(size(elec.label));
    for l = 1:numel(elec.label)
        % Find name of the probe
        probes{l} = elec.label{l}(regexp(elec.label{l},'\D'));
        % Find number of each channel
        lab_num(l) = str2num(elec.label{l}(regexp(elec.label{l},'\d')));
    end
    probes_sorted = sort(unique(probes));
    
    % Sort labels by probe then contact number
    labels = {};
    for p_ix = 1:numel(probes_sorted)
        cur_lab = elec.label(strcmp(probes,probes_sorted{p_ix}));
        [~,sort_n_idx] = sort(lab_num(strcmp(probes,probes_sorted{p_ix})));
        labels = [labels; cur_lab(sort_n_idx)];
    end
end

ordered_elec = elec;
% Order them to match labels provided
[matches,order_idx] = ismember(labels,elec.label);
if sum(matches)~=numel(labels)
    fprintf('Mismatched electrodes not found in elec:\n');
    disp(labels(matches==0));
    error('Mismatch between labels provided and elec.label!');
else
    fields = fieldnames(elec);
    for f = 1:numel(fields)
        if size(eval(['elec.' fields{f}]),1) == numel(elec.label)
            eval(['ordered_elec.' fields{f} ' = elec.' fields{f} '(order_idx(matches==1),:);']);
        elseif size(eval(['elec.' fields{f}]),1) > numel(elec.label)
            warning(['elec field "' fields{f} '" has more elements than elec.label, and will not be re-ordered!!!']);
        end 
    end 
end

end