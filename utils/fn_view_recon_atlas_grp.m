function fn_view_recon_atlas_grp(SBJs, proc_id, reg_type, show_labels,...
                                 hemi, atlas_id, roi_id, mirror, plot_out, varargin)
%% Plot a reconstruction with electrodes
% INPUTS:
%   SBJs [cell array str] - subject IDs to plot
%   proc_id [str] - name of analysis pipeline, used to pick elec file
%   plot_type [str] - {'ortho', '3d'} choose 3 slice orthogonal plot or 3D surface rendering
%   reg_type [str] - {'v', 's'} choose volume-based or surface-based registration
%   show_labels [0/1] - plot the electrode labels
%   hemi [str] - {'l', 'r', 'b'} hemisphere to plot
%   atlas_id [str] - {'DK','Dx','Yeo7','Yeo17'}
%   roi_id [str] - ROI grouping by which to color the atlas ROIs
%       'gROI','mgROI','main3' - general ROIs (lobes or broad regions)
%       'ROI','thryROI','LPFC','MPFC','OFC','INS' - specific ROIs (within these larger regions)
%       'Yeo7','Yeo17' - colored by Yeo networks
%       'tissue','tissueC' - colored by tisseu compartment, e.g., GM vs WM vs OUT
%   mirror [0/1] - plot the other hemi reflected onto single side
%   plot_out [0/1] - include electrodes that don't have an atlas label or in hemi?

[root_dir, ~] = fn_get_root_dir();

%% Handle variables
% Handle variable inputs
if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'view_angle')
            view_angle = varargin{v+1};
        elseif strcmp(varargin{v},'mesh_alpha') && varargin{v+1}>0 && varargin{v+1}<=1
            mesh_alpha = varargin{v+1};
        elseif strcmp(varargin{v},'mesh_type') && ischar(varargin{v+1})
            mesh_type = varargin{v+1};
        elseif strcmp(varargin{v},'save_fig')
            save_fig = varargin{v+1};
        elseif strcmp(varargin{v},'fig_ftype') && ischar(varargin{v+1})
            fig_ftype = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

%% Define default options
if ~exist('view_angle','var')
    view_angle = fn_get_view_angle(hemi,roi_id);
    view_str = 'def';
end
% Adjust view angle if custom
if ischar(view_angle)
    view_str = view_angle;
    switch view_angle
        case 'med'
            view_angle = fn_get_view_angle(hemi,'MPFC');
        case 'lat'
            view_angle = fn_get_view_angle(hemi,'LPFC');
        case 'ven'
            view_angle = fn_get_view_angle(hemi,'OFC');
    end
end
if show_labels; lab_arg = 'label'; else; lab_arg = 'off'; end
reg_suffix = ['_' reg_type];    % MNI space
if ~exist('mesh_type','var'); mesh_type = 'pial'; end
if ~exist('save_fig','var'); save_fig = 0; end
if ~exist('fig_ftype','var'); fig_ftype = 'fig'; end
% Assume SEEG
if ~exist('mesh_alpha','var'); mesh_alpha = 0.3; end

% ROI info
if any(strcmp(roi_id,{'mgROI','gROI','main3','lat','deep','gPFC'}))
    roi_field = 'gROI';
else
    roi_field = 'ROI';
end

%% Load elec struct
elec     = cell([numel(SBJs) 1]);
good_sbj = true(size(SBJs));
all_roi_colors = [];
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    SBJ_vars_cmd = ['run ' root_dir 'emodim/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Load MNI elecs
    elec_mni_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_mni',reg_suffix,'.mat'];
    tmp = load(elec_mni_fname); elec{sbj_ix} = tmp.elec;
    
    % Load manually adjusted gROI labels --> !!! fix this hack once I have elec_final!
    if ~strcmp(roi_id,'gROI'); error('only doing gROI now!'); end
    elec_gROI_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat_Dx_gROI.mat'];
    tmp = load(elec_gROI_fname); elec_gROI = tmp.elec;
    if numel(elec{sbj_ix}.label)~=numel(elec_gROI.label) || ~all(strcmp(elec{sbj_ix}.label,elec_gROI.label))
        error('mismatched elec labels between mni and pat_gROI!');
    end
    elec{sbj_ix}.(roi_field) = elec_gROI.roi;
    elec{sbj_ix}.color = fn_roi2color(elec_gROI.roi);
    
    % Append SBJ name to labels
    for e_ix = 1:numel(elec{sbj_ix}.label)
        elec{sbj_ix}.label{e_ix} = [SBJs{sbj_ix} '_' elec{sbj_ix}.label{e_ix}];
    end
    
%     % Match elecs to colors
%     elec{sbj_ix}.color = fn_roi2color(elec{sbj_ix}.(roi_field));
    
    if ~plot_out
        % Remove electrodes that aren't in atlas ROIs & hemisphere
        if mirror
            good_elecs = fn_select_elec_lab_match(elec{sbj_ix}, 'b', atlas_id, roi_id);
        else
            good_elecs = fn_select_elec_lab_match(elec{sbj_ix}, hemi, atlas_id, roi_id);
        end
    else
        % Remove electrodes that aren't in hemisphere
        if mirror
            good_elecs = fn_select_elec_lab_match(elec{sbj_ix}, 'b', [], []);
        else
            good_elecs = fn_select_elec_lab_match(elec{sbj_ix}, hemi, [], []);
        end
    end
    
    % Mirror hemispheres
    if mirror
        elec{sbj_ix}.chanpos(~strcmp(elec{sbj_ix}.hemi,hemi),1) = ...
                        -elec{sbj_ix}.chanpos(~strcmp(elec{sbj_ix}.hemi,hemi),1);
        hemi_str = [hemi 'b'];
    else
        hemi_str = hemi;
    end
    
    % fn_select_elec messes up if you try to toss all elecs
    if isempty(good_elecs)
        elec{sbj_ix} = {};
        good_sbj(sbj_ix) = false;
    else
        cfgs = [];
        cfgs.channel = good_elecs;
        elec{sbj_ix} = fn_select_elec(cfgs, elec{sbj_ix});
        all_roi_colors = [all_roi_colors; elec{sbj_ix}.color];
    end
    clear SBJ SBJ_vars SBJ_vars_cmd
end

%% Combine elec structs
elec = ft_appendsens([],elec{good_sbj});
elec.color = all_roi_colors;    % appendsens strips that field

%% Load brain recon
mesh = fn_load_recon_mesh([],'mni',reg_type,mesh_type,hemi);

%% 3D Surface + Grids (3d, pat/mni, vol/srf, 0/1)
out_dir = [root_dir 'emodim/results/GRP_recons/' atlas_id '_' roi_id '/'];
if ~exist(out_dir,'dir')
    [~,~] = mkdir(out_dir);
end
plot_name = ['GRP_gROI_' hemi_str '_' view_str];
f = figure('Name',plot_name);

% Plot 3D mesh
ft_plot_mesh(mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);

% Plot electrodes on top
for e = 1:numel(elec.label)
    cfgs = []; cfgs.channel = elec.label(e);
    elec_tmp = fn_select_elec(cfgs,elec);
    ft_plot_sens(elec_tmp, 'elecshape', 'sphere',...
                 'facecolor', elec_tmp.color, 'label', lab_arg);
end

view(view_angle); material dull; lighting gouraud;
l = camlight;
fprintf(['To reset the position of the camera light after rotating the figure,\n' ...
    'make sure none of the figure adjustment tools (e.g., zoom, rotate) are active\n' ...
    '(i.e., uncheck them within the figure), and then hit ''l'' on the keyboard\n'])
set(f, 'windowkeypressfcn',   @cb_keyboard);

if save_fig
    saveas(f, [out_dir plot_name '.' fig_ftype]);
end

end
