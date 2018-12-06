function fn_save_elec_atlas(SBJ, pipeline_id, view_space, reg_type, atlas_name)
%% Plot a reconstruction with electrodes
% INPUTS:
%   SBJ [str] - subject ID to plot
%   pipeline_id [str] - name of analysis pipeline, used to pick elec file
%   view_space [str] - {'pat', 'mni'}
%   reg_type [str] - {'v', 's'} choose volume-based or surface-based registration
%   atlas_name [str] - {'DK','Dx'} are the only ones implemented so far

[root_dir, ~] = fn_get_root_dir();
SBJ_vars_cmd = ['run ' root_dir 'emodim/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

view_angle = [-90 0];
if strcmp(reg_type,'v') || strcmp(reg_type,'s')
    reg_suffix = ['_' reg_type];
else
    reg_suffix = '';
end

%% Load elec struct
load([SBJ_vars.dirs.recon,SBJ,'_elec_',pipeline_id,'_',view_space,reg_suffix,'.mat']);

%% Load Atlas
fprintf('Using atlas: %s\n',atlas_name);
if strcmp(atlas_name,'DK')                  
    atlas      = ft_read_atlas(SBJ_vars.recon.fs_DK); % Desikan-Killiany (+volumetric)
    atlas.coordsys = 'acpc';
elseif strcmp(atlas_name,'Dx')
    atlas      = ft_read_atlas(SBJ_vars.recon.fs_Dx); % Destrieux (+volumetric)
    atlas.coordsys = 'acpc';
elseif strcmp(atlas_name,'Yeo7')
    atlas = fn_read_atlas(atlas_name);
    atlas.coordsys = 'mni';
elseif strcmp(atlas_name,'Yeo17')
    atlas = fn_read_atlas(atlas_name);
    atlas.coordsys = 'mni';
else
    error(['atlas_name unknown: ' atlas_name]);
end
atlas.name = atlas_name;

%% Match elecs to atlas ROIs
elec = fn_atlas_lookup(elec,atlas,'min_qry_rng',1,'max_qry_rng',5);

%% Save elec strcut with atlas labels
out_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',pipeline_id,'_',view_space,reg_suffix,'_',atlas_name,'.mat'];
fprintf('Saving %s\n',out_fname);
fprintf('==================================================================\n');
save(out_fname,'-v7.3','elec');

end
