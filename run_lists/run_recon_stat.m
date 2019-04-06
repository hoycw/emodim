pipeline_id = 'main_ft';
view_space  = 'pat';
reg_type    = '';
show_labels = 1;
hemi        = 'b';
thresh      = 0.01;
sig_method  = 'prop_boot';

SBJs = {'CP24','IR66','IR68','IR74'};
fn_view_recon_atlas_grp_stat(SBJs, pipeline_id, 'v', show_labels, 'l', 'Dx', 'gROI', 0, thresh, sig_method);
fn_view_recon_atlas_grp_stat(SBJs, pipeline_id, 'v', show_labels, 'r', 'Dx', 'gROI', 0, thresh, sig_method);

for s = 1:numel(SBJs)
    fprintf('\n====================== %s ========================\n',SBJs{s});
%     fn_view_recon_stat(SBJs{s},pipeline_id, view_space, reg_type, show_labels, 'b', 1, thresh, sig_method, 'cnts');
    fn_view_recon_stat(SBJs{s},pipeline_id, view_space, reg_type, show_labels, 'l', 0, thresh, sig_method, 'cnts');
    if ~strcmp(SBJs{s},'IR68')
        fn_view_recon_stat(SBJs{s},pipeline_id, view_space, reg_type, show_labels, 'r', 0, thresh, sig_method, 'cnts');
    end
    
%     fn_view_recon_stat(SBJs{s},pipeline_id, view_space, reg_type, show_labels, 'l', thresh, 'prop_boot', 'mcmp');
%     if ~strcmp(SBJs{s},'IR68')
%         fn_view_recon_stat(SBJs{s},pipeline_id, view_space, reg_type, show_labels, 'r', thresh, 'prop_boot', 'mcmp');
%     end
%     fn_view_recon_stat(SBJs{s},pipeline_id, view_space, reg_type, show_labels, hemi, thresh, 'norminv');
end

%% View angles
% CP24 L
view([-110 -20]);
view([90 0]);
% CP24 R
view([120 -10]);
view([-90 0]);
% IR66 L
view([-130 0]);
% IR66 R
view([130 0]);
% IR68
view([-90 0]);
view([-150 10]);
% IR74 L
view([-50 25]);
% IR74 R
view([-90 0]);