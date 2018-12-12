pipeline_id = 'main_ft';
view_space = 'pat';
reg_type = '';
show_labels =1 ;
hemi = 'b';
thresh = 0.01;

SBJs = {'CP24','IR66','IR68'};%,'IR74'};
for s = 1:numel(SBJs)
    fprintf('\n====================== %s ========================\n',SBJs{s});
    fn_view_recon_stat(SBJs{s},pipeline_id, view_space, reg_type, show_labels, 'l', thresh, 'prop_boot', 'mcmp');
    if ~strcmp(SBJs{s},'IR68')
        fn_view_recon_stat(SBJs{s},pipeline_id, view_space, reg_type, show_labels, 'r', thresh, 'prop_boot', 'mcmp');
    end
%     fn_view_recon_stat(SBJs{s},pipeline_id, view_space, reg_type, show_labels, hemi, thresh, 'norminv');
end