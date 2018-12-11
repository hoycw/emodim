function fn_view_recon_atlas_grp_sphr(SBJs, pipeline_id, reg_type, show_labels, hemi, atlas_id, roi_id, plot_out)%, view_angle)
%% Plot a reconstruction with electrodes
% INPUTS:
%   SBJ [str] - subject ID to plot
%   pipeline_id [str] - name of analysis pipeline, used to pick elec file
%   plot_type [str] - {'ortho', '3d'} choose 3 slice orthogonal plot or 3D surface rendering
%   reg_type [str] - {'v', 's'} choose volume-based or surface-based registration
%   show_labels [0/1] - plot the electrode labels
%   hemi [str] - {'l', 'r', 'b'} hemisphere to plot
%   altas_id
%   roi_id
%   plot_out [0/1] - include electrodes that don't have an atlas label?

[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
view_space = 'mni';

view_angle = [-90 0];
if strcmp(reg_type,'v') || strcmp(reg_type,'s')
    reg_suffix = ['_' reg_type];
else
    reg_suffix = '';
end
if strcmp(roi_id,'tissue') || strcmp(roi_id,'tissueC')
    tis_suffix = '_tis';
else
    tis_suffix = '';
end

% ROI info
[roi_list, ~] = fn_roi_label_styles(roi_id);
fprintf('Using atlas: %s\n',atlas_id);

%% Load elec struct
elec = cell([numel(SBJs) 1]);
good_sbj = true(size(SBJs));
all_roi_labels = {};
all_roi_colors = [];
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    SBJ_vars_cmd = ['run ' root_dir 'emodim/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    try
        if strcmp(atlas_id,'Dx') || strcmp(atlas_id,'DK')   % Cover atlases defined on SBJ surfaces
            elec_atlas_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',pipeline_id,...
                '_',view_space,reg_suffix,'.mat'];
            tmp = load(elec_atlas_fname); elec{sbj_ix} = tmp.elec;
            
            elec_atlas_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',pipeline_id,...
                '_pat_',atlas_id,tis_suffix,'.mat'];
            tmp = load(elec_atlas_fname);
            elec{sbj_ix}.atlas_label = tmp.elec.atlas_label;
            elec{sbj_ix}.atlas_name = tmp.elec.atlas_name;
            elec{sbj_ix}.hemi = tmp.elec.hemi;%!!! fix this! why IR32 is different?
        else
            elec_atlas_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',pipeline_id,...
                '_',view_space,reg_suffix,'_',atlas_id,tis_suffix,'.mat'];
            tmp = load(elec_atlas_fname); elec_sbj{sbj_ix} = tmp.elec;
        end
    catch
%         answer = input(['Could not load requested file: ' elec_atlas_fname ...
%             '\nDo you want to run the atlas matching now? "y" or "n"\n'],'s');
%         if strcmp(answer,'y')
%             fn_save_elec_atlas(SBJ,pipeline_id,view_space,reg_type,atlas_id);
%         else
            error([elec_atlas_fname 'doesnt exist, exiting...']);
%         end
    end
    
    % Append SBJ name to labels
    for e_ix = 1:numel(elec{sbj_ix}.label)
        elec{sbj_ix}.label{e_ix} = [SBJs{sbj_ix} '_' elec{sbj_ix}.label{e_ix}];
    end
    
    % Match elecs to atlas ROIs
    if any(strcmp(atlas_id,{'DK','Dx','Yeo7'}))
        elec{sbj_ix}.roi       = fn_atlas2roi_labels(elec{sbj_ix}.atlas_label,atlas_id,roi_id);
        if strcmp(roi_id,'tissueC')
            elec{sbj_ix}.roi_color = fn_tissue2color(elec{sbj_ix});
        elseif strcmp(atlas_id,'Yeo7')
            elec{sbj_ix}.roi_color = fn_atlas2color(atlas_id,elec{sbj_ix}.roi);
        else
            elec{sbj_ix}.roi_color = fn_roi2color(elec{sbj_ix}.roi);
        end
    elseif any(strcmp(atlas_id,{'Yeo17'}))
        elec{sbj_ix}.roi       = elec{sbj_ix}.atlas_label;
        elec{sbj_ix}.roi_color = fn_atlas2color(atlas_id,elec{sbj_ix}.roi);
    end
    
    % Remove electrodes that aren't in atlas ROIs
    if ~plot_out
        atlas_in_elecs = {};
        for roi_ix = 1:numel(roi_list)
            atlas_in_elecs = [atlas_in_elecs; elec{sbj_ix}.label(strcmp(elec{sbj_ix}.roi,roi_list(roi_ix)))];
        end
    else
        atlas_in_elecs = elec{sbj_ix}.label;
    end
    if ~strcmp(hemi,'b')
        hemi_in_elecs = elec{sbj_ix}.label(strcmp(elec{sbj_ix}.hemi,hemi));
    else
        hemi_in_elecs = elec{sbj_ix}.label;
    end
    % fn_select_elec messes up if you try to toss all elecs
    good_elecs = intersect(atlas_in_elecs, hemi_in_elecs);
    if numel(intersect(elec{sbj_ix}.label,good_elecs))==0
        elec{sbj_ix} = {};
        good_sbj(sbj_ix) = false;
    else
        cfgs = [];
        cfgs.channel = good_elecs;
        elec{sbj_ix} = fn_select_elec(cfgs, elec{sbj_ix});
        all_roi_labels = [all_roi_labels; elec{sbj_ix}.roi];
        all_roi_colors = [all_roi_colors; elec{sbj_ix}.roi_color];
    end
    clear SBJ SBJ_vars SBJ_vars_cmd
end

%% Combine elec structs
elec = ft_appendsens([],elec{good_sbj});
elec.roi       = all_roi_labels;    % appendsens strips that field
elec.roi_color = all_roi_colors;    % appendsens strips that field

%% Load brain recon
if strcmp(view_space,'pat')
    error('This is a group plot, not a patient plot!');
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
            mesh = ft_read_headshape([root_dir 'emodim/data/atlases/fsaverage/' hemi 'h.pial']);
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

%% 3D Surface + Grids (3d, pat/mni, vol/srf, 0/1)
h = figure;

% Plot 3D mesh
mesh_alpha = 0.2;
ft_plot_mesh(mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);

% Plot electrodes on top
for e = 1:numel(elec.label)
    cfgs = []; cfgs.channel = elec.label(e);
    elec_tmp = fn_select_elec(cfgs,elec);
    if show_labels
        ft_plot_sens(elec_tmp, 'elecshape', 'sphere',...
            'facecolor', elec_tmp.roi_color, 'label', 'label');
    else
        ft_plot_sens(elec_tmp, 'elecshape', 'sphere',...
            'facecolor', elec_tmp.roi_color);
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
