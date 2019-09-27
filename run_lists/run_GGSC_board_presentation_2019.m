%% Run list for GGSC Advisory Board Meeting Sept. 2019
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'Apps/fieldtrip/'];
elseif exist('/Users/lapate/','dir');root_dir = '/Users/lapate/knight/';ft_dir = '/Users/lapate/knight/hoycw/Apps/fieldtrip/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

SBJs    = {'CP24','IR61','IR66','IR68','IR71','IR74','IR75','IR78','IR84'};
proc_id = 'main_ft';

%%
addpath([root_dir 'emodim/scripts/']);
addpath([root_dir 'emodim/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Plot Group Recons
reg_type    = 'v';
show_labels = 0;
hemi        = 'r';
mirror      = 0;
atlas_id    = 'Dx';
roi_id      = 'gROI';
plot_out    = 0;
save_fig    = 1;
fig_ftype   = 'png';

fn_view_recon_atlas_grp(SBJs, proc_id, reg_type, show_labels,...
                        hemi, atlas_id, roi_id, mirror, plot_out,...
                        'save_fig', save_fig, 'fig_ftype', fig_ftype)
fn_view_recon_atlas_grp(SBJs, proc_id, reg_type, show_labels,...
                        hemi, atlas_id, roi_id, mirror, plot_out,...
                        'save_fig', save_fig, 'fig_ftype', fig_ftype,'view_angle','med')

%% Plot Group STAT Recons
view_space  = 'mni';
reg_type    = 'v';
show_labels = 0;
hemi        = 'l';
atlas_id    = 'Dx';
roi_id      = 'gROI';
mirror      = 1;
plot_out    = 0;
stat_id     = 'HGmLogLowpass';
thresh      = 0.01;
sig_method  = 'prop_boot';
save_fig    = 1;
fig_ftype   = 'png';

fn_view_recon_atlas_grp_stat(SBJs, proc_id, stat_id, thresh, sig_method, reg_type, show_labels,...
                                 hemi, atlas_id, roi_id, mirror, plot_out, 'save_fig', save_fig, 'fig_ftype', fig_ftype)
fn_view_recon_atlas_grp_stat(SBJs, proc_id, stat_id, thresh, sig_method, reg_type, show_labels,...
                                 hemi, atlas_id, roi_id, mirror, plot_out, 'save_fig', save_fig, 'fig_ftype', fig_ftype,...
                                 'view_angle', 'med')


% for s = 1:numel(SBJs)
%     fprintf('\n====================== %s ========================\n',SBJs{s});
% %     fn_view_recon_stat(SBJs{s},pipeline_id, view_space, reg_type, show_labels, 'b', 1, thresh, sig_method, 'cnts');
%     fn_view_recon_stat(SBJs{s},pipeline_id, view_space, reg_type, show_labels, 'l', 0, thresh, sig_method, 'cnts');
%     if ~strcmp(SBJs{s},'IR68')
%         fn_view_recon_stat(SBJs{s},pipeline_id, view_space, reg_type, show_labels, 'r', 0, thresh, sig_method, 'cnts');
%     end
% 
% %     fn_view_recon_stat(SBJs{s},pipeline_id, view_space, reg_type, show_labels, 'l', thresh, 'prop_boot', 'mcmp');
% %     if ~strcmp(SBJs{s},'IR68')
% %         fn_view_recon_stat(SBJs{s},pipeline_id, view_space, reg_type, show_labels, 'r', thresh, 'prop_boot', 'mcmp');
% %     end
% %     fn_view_recon_stat(SBJs{s},pipeline_id, view_space, reg_type, show_labels, hemi, thresh, 'norminv');
% end
