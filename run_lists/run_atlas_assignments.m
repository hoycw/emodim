% SBJs = {'CP24','IR61','IR66','IR68','IR71','IR74','IR75','IR78','IR84'};
% Jacob anat development SBJs = {'IR68','IR74','IR78','IR84'};
SBJs = {'IR61','IR71'}; % remove -ROI

%% Directories
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'hoycw/Apps/fieldtrip/'];
elseif exist('/Users/lapate/','dir');root_dir = '/Users/lapate/knight/';ft_dir = '/Users/lapate/knight/hoycw/Apps/fieldtrip/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%%
addpath([root_dir 'emodim/scripts/']);
addpath([root_dir 'emodim/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Compare to atlases
for s = 2:numel(SBJs)
    fn_compile_elec_struct(SBJs{s},'main_ft','pat','',1);
    fn_compile_elec_struct(SBJs{s},'main_ft','mni','v',1);
    fn_save_elec_atlas(SBJs{s},'main_ft','pat','','DK',1);
    fn_save_elec_atlas(SBJs{s},'main_ft','pat','','Dx',1);
%     fn_save_elec_tissue(SBJs{s},'main_ft','pat','','Dx');
%     fn_save_elec_tissue(SBJs{s},'main_ft','pat','','DK');
%     fn_save_elec_atlas(SBJs{s},'main_ft','mni','v','Yeo7');
%     fn_save_elec_atlas(SBJs{s},'main_ft','mni','v','Yeo17');
end

%% MNI Check
% SBJ = 'IR21';
% 
% % Compare patient and MNI in ortho
% fn_view_recon(SBJ,'main_ft','ortho','pat','',1,'b');
% fn_view_recon(SBJ,'main_ft','ortho','mni','v',1,'b');
% 
% % Check atlas assignments
% fn_view_recon_atlas(SBJ,pipeline_id,'pat','',1,'b','DK','gROI');
% fn_view_recon_atlas(SBJ,pipeline_id,'pat','',1,'b','Dx','gROI');
% fn_view_recon_atlas(SBJ,pipeline_id,'mni','v',1,'b','Yeo7','Yeo7');
% 
