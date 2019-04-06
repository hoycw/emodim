function [RGB] = fn_atlas2color(atlas_id,rois)
%% Returns the RGB color code to plot a given ROI
% INPUTS:
%   label [cell array of strs] - name (or list of names) of the ROI label
%       right now only Yeo7 and Yeo17
% OUTPUTS:
%   RGB [3 floats] - RGB for this ROI
%

switch atlas_id
    case 'Yeo7'
        labels = {
            'Vis'%'7Networks_1'
            'SM'%'7Networks_2'
            'DAttn'%'7Networks_3'
            'VAttn'%'7Networks_4'
            'Limb'%'7Networks_5'
            'FP'%'7Networks_6'
            'Def'%'7Networks_7'
            'OUT'
            };
        colors = [
            120  18 134;
            70 130 180;
            0 118  14;
            196  58 250;
            220 248 164;
            230 148  34;
            205  62  78;
            0 0 0
            ];
    case 'Yeo17'
        labels = {
            '17Networks_1'
            '17Networks_2'
            '17Networks_3'
            '17Networks_4'
            '17Networks_5'
            '17Networks_6'
            '17Networks_7'
            '17Networks_8'
            '17Networks_9'
            '17Networks_10'
            '17Networks_11'
            '17Networks_12'
            '17Networks_13'
            '17Networks_14'
            '17Networks_15'
            '17Networks_16'
            '17Networks_17'
            'no_label_found'
            };
        colors = [
            120  18 134;
            255   0   0;
            70 130 180;
            42 204 164;
            74 155  60;
            0 118  14;
            196  58 250;
            255 152 213;
            220 248 164;
            122 135  50;
            119 140 176;
            230 148  34;
            135  50  74;
            12  48 255;
            0   0 130;
            255 255   0;
            205  62  78
            0 0 0
            ];
end

RGB = zeros([numel(rois) 3]);
for r = 1:numel(rois)
    RGB(r,:) = colors(strcmp(rois{r},labels),:)./256;
end

end
