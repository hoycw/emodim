function fn_save_elec_tissue(SBJ, pipeline_id, view_space, reg_type, atlas_name)
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
else
    error(['atlas_name unknown: ' atlas_name]);
end
atlas.name = atlas_name;

%% Match elecs to atlas ROIs
elec = fn_atlas_lookup(elec,atlas,'min_qry_rng',5,'max_qry_rng',5);

% Check that all atlas_prob add to 1
for e = 1:numel(elec.label)
    if elec.atlas_prob(e)+sum(elec.atlas_prob2{e})<0.99999  % sometimes it's 0.99999999999999988898
        error(['Electrode ' elec.label{e} ' has atlas_prob = '...
            num2str(elec.atlas_prob(e)+sum(elec.atlas_prob2{e}))]);
    end
end

%% Convert atlas labels and probabilities to GM probability
% usedqueryrange search sizes: 1 = 1; 3 = 7; 5 = 33
elec.tissue_labels = {'GM','WM','CSF','OUT'};
elec.tissue_prob = zeros([numel(elec.label) numel(elec.tissue_labels)]);

% Assign atlas labels to tissue type
elec.tissue       = fn_atlas2roi_labels(elec.atlas_label,atlas_name,'tissue');
for e = 1:numel(elec.label)
    % Compute Probability of Tissue Types {GM, WM, CSF, OUT}
    elec.tissue_prob(e,strcmp(elec.tissue{e},elec.tissue_labels)) = ...
        elec.tissue_prob(e,strcmp(elec.tissue{e},elec.tissue_labels)) + elec.atlas_prob(e);
    
    % Check for secondary matches and add to total
    if ~isempty(elec.atlas_label2{e})
        elec.tissue2{e} = fn_atlas2roi_labels(elec.atlas_label2{e},atlas_name,'tissue');
        for roi = 1:numel(elec.tissue2{e})
            elec.tissue_prob(e,strcmp(elec.tissue2{e}{roi},elec.tissue_labels)) = ...
                elec.tissue_prob(e,strcmp(elec.tissue2{e}{roi},elec.tissue_labels)) + elec.atlas_prob2{e}(roi);
        end
    end
end

%% Save elec strcut with atlas labels
out_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',pipeline_id,'_',view_space,reg_suffix,'_',atlas_name,'_tis.mat'];
fprintf('Saving %s\n',out_fname);
fprintf('==================================================================\n');
save(out_fname,'-v7.3','elec');

end
