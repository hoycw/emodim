function [labels, colors] = fn_roi_label_styles(roi_id)
%% Converts the name of a set of ROIs into labels, plotting colors
% condition_name: [str] ROI label to determine which set of sub-ROIs to load
%   'ROI'- all individual ROIs (at smallest scale of atlas)
%   'gROI'- general ROIs, AKA lobes ('INS', 'LPFC', 'MPFC', 'OFC', 'PAR', 'TMP')
%   any gROI - specific sub-ROIs in that general ROI
% colors from http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=3
%
% % einfo_col is now defunct! % einfo_col = 2 for specific ROIs, 3 for general ROIs

% if length(cond_lab) == 1
switch roi_id
    case 'ROI'
        load('~/PRJ_Stroop/data/full_roi_lists.mat');
        labels = all_rois;
        % Exclude FWM, '', OUT
        labels(strmatch('FWM',labels,'exact')) = [];
        labels(strmatch('TWM',labels,'exact')) = [];
        labels(strmatch('OUT',labels,'exact')) = [];
        labels(strmatch('',labels,'exact')) = [];
        % einfo_col = 2;
    case 'Yeo7'
        labels = {'Vis','SM','DAttn','VAttn','Limb','FP','Def'};
    case 'Main3'
        labels = {'LPFC','MPFC','INS'};
        % einfo_col = 3;
    case 'mgROI'
        labels = {'LPFC','MPFC','INS','OFC'};
        % einfo_col = 3;
    case 'gROI'
        labels = {'LPFC','MPFC','INS','OFC','PAR','TMP','MTL','OCC'};
        % einfo_col = 3;
    case 'mnLPFC'
        labels = {'DLPFC','VLPFC','PM','aMCC','preSMA','SMA'};
        % einfo_col = 2;
    case 'thryROI'
        labels = {'DLPFC','VLPFC','PM','aMCC','preSMA','SMA','daINS','vaINS','FO'};
        % einfo_col = 2;
    case 'PAR'
        labels = {'S1','SPL','IPL','Precuneus'};
    case 'TMP'
        labels = {'STS'};
    case 'LPFC'
        labels = {'FPC','DLPFC','VLPFC','PM','M1'};
        % einfo_col = 2;
    case 'MPFC'
        labels = {'ACC','preSMA','aMCC','SMA','pMCC'};
        % einfo_col = 2;
    case 'INS'
        labels = {'vaINS','daINS','FO','mINS','pINS'};
        % einfo_col = 2;
    case 'OFC'
        labels = {'mOFC','lOFC'};
        % einfo_col = 2;
    case 'MTL'
        labels = {'HPC','AMG'};
    case {'tissue', 'tissueC'}
        labels = {'GM','WM','CSF','OUT'};
        % einfo_col = [];
    case 'all'
        load('~/PRJ_Stroop/data/full_roi_lists.mat');
        labels = all_rois;
        % einfo_col = 2;
    otherwise
        error(strcat('Unknown roi_id: ',roi_id));
end

% Get colors
colors = cell(size(labels));
for roi_ix = 1:numel(labels)
    colors{roi_ix} = fn_roi2color(labels{roi_ix});
end

end
