function fn_view_recon(SBJ, pipeline_id, plot_type, view_space, reg_type, show_labels, hemi)%, view_angle)
%% Plot a reconstruction with electrodes
% INPUTS:
%   SBJ [str] - subject ID to plot
%   pipeline_id [str] - name of analysis pipeline, used to pick elec file
%   plot_type [str] - {'ortho', '3d'} choose 3 slice orthogonal plot or 3D surface rendering
%   view_space [str] - {'pat', 'mni'}
%   reg_type [str] - {'v', 's'} choose volume-based or surface-based registration
%   show_labels [0/1] - plot the electrode labels
%   hemi [str] - {'l', 'r', 'b'} hemisphere to plot

[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

view_angle = [-90 0];
if strcmp(reg_type,'v') || strcmp(reg_type,'s')
    reg_suffix = ['_' reg_type];
else
    reg_suffix = '';
end

%% Load elec struct
load([SBJ_vars.dirs.recon,SBJ,'_elec_',pipeline_id,'_',view_space,reg_suffix,'.mat']);

%% Load brain recon
if strcmp(view_space,'pat')
    if strcmp(plot_type,'3d')
        if strcmp(hemi,'r') || strcmp(hemi,'l')
            mesh = ft_read_headshape([SBJ_vars.dirs.recon 'Surfaces/' SBJ '_cortex_' hemi 'h.mat']);
        elseif strcmp(hemi,'b')
            mesh = ft_read_headshape({[SBJ_vars.dirs.recon 'Surfaces/' SBJ '_cortex_rh.mat'],...
                [SBJ_vars.dirs.recon 'Surfaces/' SBJ '_cortex_lh.mat']});
        else
            error(['Unknown hemisphere selected: ' hemi]);
        end
        mesh.coordsys = 'acpc';
    elseif strcmp(plot_type,'ortho')
        mri = ft_read_mri(SBJ_vars.recon.fs_T1);
        mri.coordsys = 'acpc';
    else
        error(['Unknown plot_type: ' plot_type]);
    end
elseif strcmp(view_space,'mni')
    if strcmp(plot_type,'3d')
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
                mesh = ft_read_headshape([root_dir 'PRJ_Stroop/data/atlases/freesurfer/fsaverage/' hemi 'h.pial']);
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
    elseif strcmp(plot_type,'ortho')
        if strcmp(reg_type,'v')
            mri = ft_read_mri([ft_dir 'template/anatomy/single_subj_T1_1mm.nii']);
            mri.coordsys = 'mni';
        elseif strcmp(reg_type,'s')
            error('ortho plot with surface based registration doesnt make sense!');
        end
    else
        error(['Unknown plot_type: ' plot_type]);
    end
else
    error(['Unknown view_space: ' view_space]);
end

%% Orthoplot (pat/mni, v only, 0/1 labels)
if strcmp(plot_type,'ortho')
    % ft_electrodeplacement only plots elec.elecpos, so swap in chanpos
    elec.elecpos = elec.chanpos;
    if isfield(elec,'tra')
        elec = rmfield(elec, 'tra');    % error if elec.tra shows the difference between original elecpos and new chanpos post-reref
    end
    cfg = [];
    cfg.elec = elec;
    ft_electrodeplacement(cfg, mri);
%     ft_plot_ortho(mri.anatomy, 'transform', mri.transform);%, 'style', 'intersect');
%     if show_labels
%         ft_plot_sens(elec, 'label', 'on', 'fontcolor', 'w');
%     else
%         ft_plot_sens(elec, 'label', 'off');
%     end
end

%% 3D Surface + Grids (3d, pat/mni, v/s, 0/1)
if strcmp(plot_type,'3d')
    h = figure;
    
    % Plot 3D mesh
    mesh_alpha = 0.8;
    if any(strcmp(SBJ_vars.ch_lab.probe_type,'seeg'))
        mesh_alpha = 0.4;
    end
    ft_plot_mesh(mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);
    
    % Plot electrodes on top
    if show_labels
        ft_plot_sens(elec, 'elecshape', 'sphere', 'label', 'label');
    else
        ft_plot_sens(elec, 'elecshape', 'sphere');
    end
    
    view(view_angle); material dull; lighting gouraud;
    l = camlight;
    fprintf(['To reset the position of the camera light after rotating the figure,\n' ...
        'make sure none of the figure adjustment tools (e.g., zoom, rotate) are active\n' ...
        '(i.e., uncheck them within the figure), and then hit ''l'' on the keyboard\n'])
    set(h, 'windowkeypressfcn',   @cb_keyboard);
end

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
