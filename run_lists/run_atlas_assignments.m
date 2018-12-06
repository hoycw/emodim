SBJs = {'IR68'};

%% Create basic elec structs
fn_compile_elec_struct(SBJs{1},'main_ft','pat','');
fn_compile_elec_struct(SBJs{1},'main_ft','mni','v');

%% Compare to atlases
for s = 1:numel(SBJs)
    fn_save_elec_atlas(SBJs{s},'main_ft','pat','','DK');
    fn_save_elec_atlas(SBJs{s},'main_ft','pat','','Dx');
    fn_save_elec_tissue(SBJs{s},'main_ft','pat','','Dx');
    fn_save_elec_tissue(SBJs{s},'main_ft','pat','','DK');
    fn_save_elec_atlas(SBJs{s},'main_ft','mni','v','Yeo7');
    fn_save_elec_atlas(SBJs{s},'main_ft','mni','v','Yeo17');
end

%% MNI Check
SBJ = 'IR21';

% Compare patient and MNI in ortho
fn_view_recon(SBJ,'main_ft','ortho','pat','',1,'b');
fn_view_recon(SBJ,'main_ft','ortho','mni','v',1,'b');

% Check atlas assignments
fn_view_recon_atlas(SBJ,pipeline_id,'pat','',1,'b','DK','gROI');
fn_view_recon_atlas(SBJ,pipeline_id,'pat','',1,'b','Dx','gROI');
fn_view_recon_atlas(SBJ,pipeline_id,'mni','v',1,'b','Yeo7','Yeo7');

