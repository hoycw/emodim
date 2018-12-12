function fn_view_recon_atlas_gROItmp(SBJ, pipeline_id, view_space, reg_type, show_labels, hemi, atlas_name, roi_style, plot_out)%, view_angle)
%% Plot a reconstruction with electrodes
% INPUTS:
%   SBJ [str] - subject ID to plot
%   pipeline_id [str] - name of analysis pipeline, used to pick elec file
%   plot_type [str] - {'ortho', '3d'} choose 3 slice orthogonal plot or 3D surface rendering
%   view_space [str] - {'pat', 'mni'}
%   reg_type [str] - {'v', 's'} choose volume-based or surface-based registration
%   show_labels [0/1] - plot the electrode labels
%   hemi [str] - {'l', 'r', 'b'} hemisphere to plot
%   plot_out [0/1] - exclude electrodes that don't match atlas or aren't in hemisphere
%   atlas_name [str] - {'DK','Dx','Yeo7','Yeo17'}
%   roi_style [str] - ROI grouping by which to color the atlas ROIs
%       'gROI','mgROI','main3' - general ROIs (lobes or broad regions)
%       'ROI','thryROI','LPFC','MPFC','OFC','INS' - specific ROIs (within these larger regions)
%       'Yeo7','Yeo17' - colored by Yeo networks
%       'tissue','tissueC' - colored by tisseu compartment, e.g., GM vs WM vs OUT

[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
SBJ_vars_cmd = ['run ' root_dir 'emodim/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

view_angle = [-90 0];
if strcmp(reg_type,'v') || strcmp(reg_type,'s')
    reg_suffix = ['_' reg_type];
else
    reg_suffix = '';
end
if strcmp(roi_style,'tissue') || strcmp(roi_style,'tissueC')
    tis_suffix = '_tis';
else
    tis_suffix = '';
end

%% Load elec struct
try
    elec_atlas_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',pipeline_id,'_',view_space,reg_suffix,'_',atlas_name,tis_suffix,'_',roi_style,'.mat'];
    load(elec_atlas_fname);
    
catch
    answer = input(['Could not load requested file: ' elec_atlas_fname ...
        '\nDo you want to run the atlas matching now? "y" or "n"\n'],'s');
    if strcmp(answer,'y')
        fn_save_elec_atlas(SBJ,pipeline_id,view_space,reg_type,atlas_name);
    else
        error('not running atlas assignment, exiting...');
    end
end

%% Remove electrodes that aren't in atlas ROIs
if ~plot_out
    atlas_out_elecs = elec.label(strcmp(elec.atlas_label,'no_label_found'));
    if ~strcmp(hemi,'b')
        hemi_out_elecs = elec.label(~strcmp(elec.hemi,hemi));
    end
    cfgs = []; cfgs.channel = [{'all'} fn_ch_lab_negate(atlas_out_elecs) fn_ch_lab_negate(hemi_out_elecs)];
    elec = fn_select_elec(cfgs, elec);
end

%% Load brain recon
if strcmp(view_space,'pat')
    if strcmp(hemi,'r') || strcmp(hemi,'l')
        mesh = ft_read_headshape([SBJ_vars.dirs.recon 'Surfaces/' SBJ '_cortex_' hemi 'h.mat']);
    elseif strcmp(hemi,'b')
        mesh = ft_read_headshape({[SBJ_vars.dirs.recon 'Surfaces/' SBJ '_cortex_rh.mat'],...
                                    [SBJ_vars.dirs.recon 'Surfaces/' SBJ '_cortex_lh.mat']});
    else
        error(['Unknown hemisphere selected: ' hemi]);
    end
    mesh.coordsys = 'acpc';
elseif strcmp(view_space,'mni')
    if strcmp(reg_type,'v')
        if strcmp(hemi,'r')
            load([ft_dir 'template/anatomy/surface_pial_right.mat']);
        elseif strcmp(hemi,'l')
            load([ft_dir 'template/anatomy/surface_pial_left.mat']);
        elseif strcmp(hemi,'b')
            load([ft_dir 'template/anatomy/surface_pial_both.mat']);
        else
            error(['Unknown hemisphere option: ' hemi]);
        end
%         mesh.coordsys = 'mni';
    elseif strcmp(reg_type,'s')
        if strcmp(hemi,'r') || strcmp(hemi,'l')
            mesh = ft_read_headshape([root_dir 'emodim/data/atlases/freesurfer/fsaverage/' hemi 'h.pial']);
        elseif strcmp(hemi,'b')
            error('hemisphere "b" not yet implemented for reg_type: "srf"!');
            mesh = ft_read_headshape([ft_dir 'subjects/fsaverage/surf/' hemi 'h.pial']);
        else
            error(['Unknown hemisphere option: ' hemi]);
        end
        mesh.coordsys = 'fsaverage';
    else
        error(['Unknown registration type (reg_type): ' reg_type]);
    end
else
    error(['Unknown view_space: ' view_space]);
end

%% Load Atlas
fprintf('Using atlas: %s\n',atlas_name);
% if strcmp(atlas_name,'DK')                  
%     atlas      = ft_read_atlas(SBJ_vars.recon.fs_DK); % Desikan-Killiany (+volumetric)
%     atlas.coordsys = 'acpc';
% elseif strcmp(atlas_name,'Dx')
%     atlas      = ft_read_atlas(SBJ_vars.recon.fs_Dx); % Destrieux (+volumetric)
%     atlas.coordsys = 'acpc';
% elseif strcmp(atlas_name,'Yeo7')
%     atlas = fn_read_atlas(atlas_name);
%     atlas.coordsys = 'mni';
% elseif strcmp(atlas_name,'Yeo17')
%     atlas = fn_read_atlas(atlas_name);
%     atlas.coordsys = 'mni';
% else
%     error(['atlas_name unknown: ' atlas_name]);
% end
% atlas.name = atlas_name;
% % elec.elecpos_fs   = elec.elecpos;

%% Match elecs to atlas ROIs
if any(strcmp(atlas_name,{'DK','Dx','Yeo7'}))
%     elec.roi       = fn_atlas2roi_labels(elec.atlas_label,atlas_name,roi_style);
    if strcmp(roi_style,'tissueC')
        elec.roi_color = fn_tissue2color(elec);
    elseif strcmp(atlas_name,'Yeo7')
        elec.roi_color = fn_atlas2color(atlas_name,elec.roi);
    else
        elec.roi_color = fn_roi2color(elec.roi);
    end
elseif any(strcmp(atlas_name,{'Yeo17'}))
    elec.roi       = elec.atlas_label;
    elec.roi_color = fn_atlas2color(atlas_name,elec.roi);
end

%% 3D Surface + Grids (3d, pat/mni, vol/srf, 0/1)
h = figure;

% Plot 3D mesh
mesh_alpha = 0.7;
if any(strcmp(SBJ_vars.ch_lab.probe_type,'seeg'))
    mesh_alpha = 0.2;
end
ft_plot_mesh(mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);

% Plot electrodes on top
cfgs = [];
for e = 1:numel(elec.label)
    cfgs.channel = elec.label{e};
    elec_tmp = fn_select_elec(cfgs, elec);
    if show_labels
        ft_plot_sens(elec_tmp, 'elecshape', 'sphere', 'facecolor', elec_tmp.roi_color, 'label', 'label');
    else
        ft_plot_sens(elec_tmp, 'elecshape', 'sphere', 'facecolor', elec_tmp.roi_color);
    end
end

view(view_angle); material dull; lighting gouraud;
l = camlight;
fprintf(['To reset the position of the camera light after rotating the figure,\n' ...
    'make sure none of the figure adjustment tools (e.g., zoom, rotate) are active\n' ...
    '(i.e., uncheck them within the figure), and then hit ''l'' on the keyboard\n'])
set(h, 'windowkeypressfcn',   @cb_keyboard);

%% Plot SEEG data in 3D
% % Create volumetric mask of ROIs from fs parcellation/segmentation
% atlas = ft_read_atlas([SBJ_dir 'freesurfer/mri/aparc+aseg.mgz']);
% atlas.coordsys = 'acpc';
% cfg = [];
% cfg.inputcoord = 'acpc';
% cfg.atlas = atlas;
% cfg.roi = {'Right-Hippocampus', 'Right-Amygdala'};
% mask_rha = ft_volumelookup(cfg, atlas);

%% EXTRA CRAP:
% % Plot HFA from bipolar channels via clouds around electrode positions
% cfg = [];
% cfg.funparameter = 'powspctrm';
% cfg.funcolorlim = [-.5 .5];
% cfg.method = 'cloud';
% cfg.slice = '3d';
% cfg.nslices = 2;
% cfg.facealpha = .25;
% ft_sourceplot(cfg, freq_sel2, mesh_rha);
% view([120 40]); lighting gouraud; camlight;
% 
% % 2D slice version:
% cfg.slice = '2d';
% ft_sourceplot(cfg, freq_sel2, mesh_rha);
% 
% %% View grid activity on cortical mesh
% cfg = [];
% cfg.funparameter = 'powspctrm';
% cfg.funcolorlim = [-.5 .5];
% cfg.method = 'surface';
% cfg.interpmethod = 'sphere_weighteddistance';
% cfg.sphereradius = 8;
% cfg.camlight = 'no';
% ft_sourceplot(cfg, freq_sel, pial_lh);
% 
% %% Prepare and plot 2D layout
% % Make layout
% cfg = [];
% cfg.headshape = pial_lh;
% cfg.projection = 'orthographic';
% cfg.channel = {'LPG*', 'LTG*'};
% cfg.viewpoint = 'left';
% cfg.mask = 'convex';
% cfg.boxchannel = {'LTG30', 'LTG31'};
% lay = ft_prepare_layout(cfg, freq);
% % Plot interactive
% cfg = [];
% cfg.layout = lay;
% cfg.showoutline = 'yes';
% ft_multiplotTFR(cfg, freq_blc);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_keyboard(h, eventdata)

if isempty(eventdata)
  % determine the key that corresponds to the uicontrol element that was activated
  key = get(h, 'userdata');
else
  % determine the key that was pressed on the keyboard
  key = parseKeyboardEvent(eventdata);
end
% get focus back to figure
if ~strcmp(get(h, 'type'), 'figure')
  set(h, 'enable', 'off');
  drawnow;
  set(h, 'enable', 'on');
end

if strcmp(key, 'l') % reset the light position
  delete(findall(h,'Type','light')) % shut out the lights
  camlight; lighting gouraud; % add a new light from the current camera position
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = parseKeyboardEvent(eventdata)

key = eventdata.Key;
% handle possible numpad events (different for Windows and UNIX systems)
% NOTE: shift+numpad number does not work on UNIX, since the shift
% modifier is always sent for numpad events
if isunix()
  shiftInd = match_str(eventdata.Modifier, 'shift');
  if ~isnan(str2double(eventdata.Character)) && ~isempty(shiftInd)
    % now we now it was a numpad keystroke (numeric character sent AND
    % shift modifier present)
    key = eventdata.Character;
    eventdata.Modifier(shiftInd) = []; % strip the shift modifier
  end
elseif ispc()
  if strfind(eventdata.Key, 'numpad')
    key = eventdata.Character;d
  end
end

if ~isempty(eventdata.Modifier)
  key = [eventdata.Modifier{1} '+' key];
end
