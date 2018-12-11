function fn_view_recon_stat(SBJ, pipeline_id, stat_id, an_id, view_space, reg_type, show_labels, hemi)
%% Plot a reconstruction with electrodes colored according to statistics
%   FUTURE 1: this is a static brain, need to adapt to a movie!
%   FUTURE 2: add option for stat_var to be a cell with 2nd stat for edge
% INPUTS:
%   SBJ [str] - subject ID to plot
%   pipeline_id [str] - name of analysis pipeline, used to pick elec file
%   stat_id [str] - ID of the stats
%       'actv': red for active, blue for deactive, yellow for both
%       NOPE: 'CI': inc vs. con via ft statistics (not run for all patients!)
%       'RT': correlation with RT (red for significant)
%       'CNI': ANOVA of congruence (red for sig)
%       'pcon': ANOVA of proportion congruence (red for sig)
%   an_id [str] - analysis ID for preprocessing, filtering, etc.
%   view_space [str] - {'pat', 'mni'}
%   reg_type [str] - {'v', 's'} choose volume-based or surface-based registration
%   show_labels [0/1] - plot the electrode labels
%   hemi [str] - {'l', 'r', 'b'} hemisphere to plot

[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
SBJ_vars_cmd = ['run ' root_dir 'emodim/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

ns_color = [0 0 0];
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
    if strcmp(reg_type,'vol')
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
    elseif strcmp(reg_type,'srf')
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

%% Load Stats
% Determine options: {'actv','CI','RT','CNI','pcon'}
if strcmp(stat_id,'actv')
    load([SBJ_vars.dirs.proc SBJ '_actv_ROI_' an_id '.mat']);
    elec = fn_reorder_elec(elec,hfa.label);
    grp_lab = {};
%     grp_colors = {'k','r','b'};
    elec_colors = cell([numel(hfa.label) 1]);
    for ch_ix = 1:numel(hfa.label)
        if any(strcmp(actv_ch,hfa.label{ch_ix}))
            actv_epochs = actv_ch_epochs{strcmp(actv_ch,hfa.label{ch_ix})};
            epoch_signs = zeros([size(actv_epochs,1) 1]);
            for ep_ix = 1:size(actv_epochs,1)
                sig_chunk_ix = [find(hfa.time==actv_epochs(ep_ix,1))...
                    find(hfa.time==actv_epochs(ep_ix,2))];
                % Find sign of (de)activation
                if 0<=squeeze(mean(mean(hfa.powspctrm(:,ch_ix,1,sig_chunk_ix(1):sig_chunk_ix(2)),1),4))
                    epoch_signs(ep_ix) = 1;
                else
                    epoch_signs(ep_ix) = -1;
                end
            end
            if any(epoch_signs==1) && any(epoch_signs==-1)
                elec_colors{ch_ix,1} = [0.75 0 0.75];
            elseif any(epoch_signs==1)
                elec_colors{ch_ix,1} = [1 0 0];
            elseif any(epoch_signs==-1)
                elec_colors{ch_ix,1} = [0 0 1];
            end
        else
            elec_colors{ch_ix,1} = ns_color;
        end
    end
elseif strcmp(stat_id,'CSE')
    load([SBJ_vars.dirs.proc,SBJ,'_',stat_id,'_ROI_',an_id,'.mat']);
    elec = fn_reorder_elec(elec,stat.label);
    grp_lab = {};
	[~, cond_colors, ~] = fn_condition_label_styles(stat_id);
    elec_colors = cell([numel(stat.label) 1]);
    for ch_ix = 1:numel(stat.label)
        if any(stat.mask(ch_ix,1,:))
            sig_epoch = find(stat.mask(ch_ix,1,:));
            if any(diff(sig_epoch)>1)
                error([SBJ ' has multiple sig epochs in channel ' stat.label{ch_ix}]);
            end
            % Find sign of (de)activation
            cI_mean = squeeze(mean(mean(hfa{1}.powspctrm(:,ch_ix,1,sig_epoch),1),4));
            iI_mean = squeeze(mean(mean(hfa{2}.powspctrm(:,ch_ix,1,sig_epoch),1),4));
            if cI_mean>=iI_mean
                elec_colors{ch_ix,1} = cond_colors{1};
            else
                elec_colors{ch_ix,1} = cond_colors{2};
            end
        else
            elec_colors{ch_ix,1} = ns_color;
        end
    end    
else    % ANOVA
    eval(['run ' root_dir 'emodim/scripts/stat_vars/' stat_id '_vars.m']);
    [grp_lab, grp_colors, ~] = fn_group_label_styles(model_lab);
    [rt_lab, rt_color, ~]    = fn_group_label_styles('RT');
    % % Load RTs
    % load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
    
    f_name = [SBJ_vars.dirs.proc SBJ '_ANOVA_ROI_' stat_id '_' an_id '.mat'];
    load(f_name,'stat','w2');
    elec = fn_reorder_elec(elec,stat.label);
    
    % FDR correct pvalues for ANOVA
%     win_lim = {}; win_center = {};
    elec_colors = cell([numel(w2.label) numel(grp_lab)+1]);
    for ch_ix = 1:numel(stat.label)
        pvals = squeeze(w2.pval(:,ch_ix,:));
        [~, ~, ~, qvals] = fdr_bh(pvals);%,0.05,'pdep','yes');
        
        % Consolidate to binary sig/non-sig
        for grp_ix = 1:numel(grp_lab)
            if any(qvals(grp_ix,:)<0.05,2)
                elec_colors{ch_ix,grp_ix} = grp_colors{grp_ix};  % sig
            else
                elec_colors{ch_ix,grp_ix} = ns_color;  % non-sig
            end
        end
        if any(stat.mask(ch_ix,1,:))
            elec_colors{ch_ix,numel(grp_lab)+1} = rt_color{:};  % sig
        else
            elec_colors{ch_ix,numel(grp_lab)+1} = ns_color;  % non-sig
        end
    end
    
%     % Get Sliding Window Parameters
%     win_lim    = fn_sliding_window_lim(stat.time,win_len,win_step);
%     win_center = round(mean(win_lim,2));
    
%     % Convert % explained variance to 0-100 scale
%     w2.trial = w2.trial*100;
    
    
%     % Trim data to plotting epoch
%     %   NOTE: stat should be on stat_lim(1):stat_lim(2)+0.001 time axis
%     %   w2 should fit within that since it's averaging into a smaller window
%     cfg_trim = [];
%     cfg_trim.latency = plt_vars.plt_lim_SR;
%     stat = ft_selectdata(cfg_trim,stat);
end

%% 3D Surface + Grids (3d, pat/mni, vol/srf, 0/1)
f = {};
for grp_ix = 1:numel(grp_lab)+1
    if strcmp(stat_id,'actv')
        plot_name = [SBJ '_actv_HFA_' an_id];
    elseif strcmp(stat_id,'CSE')
        plot_name = [SBJ '_' stat_id '_' an_id];
    else
        if grp_ix<=numel(grp_lab)
            plot_name = [SBJ '_ANOVA_' grp_lab{grp_ix} '_' stat_id '_' an_id];
        else
            plot_name = [SBJ '_ANOVA_' rt_lab{1} '_' stat_id '_' an_id];
        end
    end
    f{grp_ix} = figure('Name',plot_name);
        
    % Plot 3D mesh
    mesh_alpha = 0.8;
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
            lab_arg = 'label';
        else
            lab_arg = 'off';
        end
        ft_plot_sens(elec_tmp, 'elecshape', 'sphere', 'facecolor', elec_colors{e,grp_ix}, 'label', lab_arg);
    end
    
    view(view_angle); material dull; lighting gouraud;
    l = camlight;
    fprintf(['To reset the position of the camera light after rotating the figure,\n' ...
        'make sure none of the figure adjustment tools (e.g., zoom, rotate) are active\n' ...
        '(i.e., uncheck them within the figure), and then hit ''l'' on the keyboard\n'])
    set(f{grp_ix}, 'windowkeypressfcn',   @cb_keyboard);
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
