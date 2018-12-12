pipeline_id = 'main_ft';
view_space = 'pat';
reg_type = '';
show_labels =1 ;
hemi = 'b';
thresh = 0.01;

SBJs = {'IR74'};%{'CP24','IR66','IR68','IR74'};
for s = 4%1:numel(SBJs)
    fprintf('\n====================== %s ========================\n',SBJs{s});
    fn_view_recon_stat(SBJs{s},pipeline_id, view_space, reg_type, show_labels, 'l', thresh, 'prop_boot');
    if ~strcmp(SBJs{s},'IR68')
        fn_view_recon_stat(SBJs{s},pipeline_id, view_space, reg_type, show_labels, 'l', thresh, 'prop_boot');
    end
%     fn_view_recon_stat(SBJs{s},pipeline_id, view_space, reg_type, show_labels, hemi, thresh, 'norminv');
end