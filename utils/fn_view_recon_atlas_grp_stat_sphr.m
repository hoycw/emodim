function fn_view_recon_atlas_grp_stat_sphr(SBJs, pipeline_id, ...%stat_id, an_id, ...
                        reg_type, show_labels, hemi, atlas_id, roi_id, plot_out)%, view_angle)
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

%% Process stat_id
eval(['run ' root_dir 'emodim/scripts/stat_vars/' stat_id '_vars.m']);
if strcmp(stat_id,'actv') || strcmp(stat_id,'CSE')
    cond_lab = stat_id;
elseif strcmp(stat_id,'corrRT_CNI_pcon_WL200_WS50')
    % Get condition info
    [grp_lab, ~, ~] = fn_group_label_styles(model_lab);
    % if rt_correlation
    [rt_lab, ~, ~]     = fn_group_label_styles('RT');
    % end
    cond_lab = [grp_lab rt_lab];
else
    error(['Unknown stat_id: ' stat_id]);
end

% Prep report
out_dir = [root_dir 'emodim/results/HFA/'];
sig_report_fname = [out_dir 'GRP_' stat_id '_' an_id '_' atlas_id '_' roi_id '_sig_report.txt'];
if exist(sig_report_fname)
    system(['mv ' sig_report_fname ' ' sig_report_fname(1:end-4) '_bck.txt']);
end
sig_report = fopen(sig_report_fname,'a');

%% Load Atlas
% ROI info
[roi_list, ~] = fn_roi_label_styles(roi_id);

fprintf('Using atlas: %s\n',atlas_id);
% % if strcmp(atlas_id,'DK')                  
% %     atlas      = ft_read_atlas(SBJ_vars.recon.fs_DK); % Desikan-Killiany (+volumetric)
% %     atlas.coordsys = 'acpc';
% % elseif strcmp(atlas_id,'Dx')
% %     atlas      = ft_read_atlas(SBJ_vars.recon.fs_Dx); % Destrieux (+volumetric)
% %     atlas.coordsys = 'acpc';
% if strcmp(atlas_id,'Yeo7')
%     atlas = fn_read_atlas(atlas_id);
%     atlas.coordsys = 'mni';
% elseif strcmp(atlas_id,'Yeo17')
%     atlas = fn_read_atlas(atlas_id);
%     atlas.coordsys = 'mni';
% % else
% %     error(['atlas_name unknown: ' atlas_id]);
% end
% atlas.name = atlas_id;
% % elec.elecpos_fs   = elec.elecpos;

%% Load elec struct
elec_sbj = cell([numel(SBJs) numel(cond_lab)]);
good_sbj = true([numel(SBJs) numel(cond_lab)]);
all_roi_labels = cell([numel(cond_lab) 1]);
all_roi_colors = cell([numel(cond_lab) 1]);
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    SBJ_vars_cmd = ['run ' root_dir 'emodim/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    try
        if strcmp(atlas_id,'Dx') || strcmp(atlas_id,'DK')   % Cover atlases defined on SBJ surfaces
            elec_atlas_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',pipeline_id,...
                '_',view_space,reg_suffix,'.mat'];
            tmp = load(elec_atlas_fname); elec_sbj{sbj_ix,1} = tmp.elec;
            
            elec_atlas_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',pipeline_id,...
                '_pat_',atlas_id,tis_suffix,'.mat'];
            tmp = load(elec_atlas_fname);
            elec_sbj{sbj_ix,1}.atlas_label = tmp.elec.atlas_label;
            elec_sbj{sbj_ix,1}.atlas_name = tmp.elec.atlas_name;
        else
            elec_atlas_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',pipeline_id,...
                '_',view_space,reg_suffix,'_',atlas_id,tis_suffix,'.mat'];
            tmp = load(elec_atlas_fname); elec_sbj{sbj_ix,1} = tmp.elec;
        end
        if strcmp(atlas_id,'Dx')
            elec_atlas_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',pipeline_id,...
                '_pat_',atlas_id,tis_suffix,'.mat'];
        else
            elec_atlas_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',pipeline_id,...
                '_',view_space,reg_suffix,'_',atlas_id,tis_suffix,'.mat'];
        end
        tmp = load(elec_atlas_fname); elec_sbj{sbj_ix,1} = tmp.elec;
    catch
%         answer = input(['Could not load requested file: ' elec_atlas_fname ...
%             '\nDo you want to run the atlas matching now? "y" or "n"\n'],'s');
%         if strcmp(answer,'y')
%             fn_save_elec_atlas(SBJ,pipeline_id,view_space,reg_type,atlas_id);
%         else
        error([elec_atlas_fname 'doesnt exist, exiting...']);
%         end
    end
    
    % Update labels with SBJ name
    for e_ix = 1:numel(elec_sbj{sbj_ix,1}.label)
        elec_sbj{sbj_ix,1}.label{e_ix} = [SBJs{sbj_ix} '_' elec_sbj{sbj_ix,1}.label{e_ix}];
    end
    
    % Match elecs to atlas ROIs
    if any(strcmp(atlas_id,{'DK','Dx','Yeo7'}))
        elec_sbj{sbj_ix,1}.roi       = fn_atlas2roi_labels(elec_sbj{sbj_ix,1}.atlas_label,atlas_id,roi_id);
        if strcmp(roi_id,'tissueC')
            elec_sbj{sbj_ix,1}.roi_color = fn_tissue2color(elec_sbj{sbj_ix,1});
        elseif strcmp(atlas_id,'Yeo7')
            elec_sbj{sbj_ix,1}.roi_color = fn_atlas2color(atlas_id,elec_sbj{sbj_ix,1}.roi);
        else
            elec_sbj{sbj_ix,1}.roi_color = fn_roi2color(elec_sbj{sbj_ix,1}.roi);
        end
    elseif any(strcmp(atlas_id,{'Yeo17'}))
        elec_sbj{sbj_ix,1}.roi       = elec_sbj{sbj_ix,1}.atlas_label;
        elec_sbj{sbj_ix,1}.roi_color = fn_atlas2color(atlas_id,elec_sbj{sbj_ix,1}.roi);
    end
    
    % Copy for other conditions
    for cond_ix = 2:numel(cond_lab)
        elec_sbj{sbj_ix,cond_ix} = elec_sbj{sbj_ix,1};
    end
    
    % Load Stats
    % Determine options: {'actv','CI','RT','CNI','pcon'}
    sig_ch = cell(size(cond_lab));
    if strcmp(stat_id,'actv')
        load([SBJ_vars.dirs.proc SBJ '_actv_ROI_' an_id '.mat'],'actv_ch');
        sig_ch{1} = actv_ch;
        clear actv_ch
    elseif strcmp(stat_id,'CSE')
        load([SBJ_vars.dirs.proc,SBJ,'_',stat_id,'_ROI_',an_id,'.mat'],'stat');
        for ch_ix = 1:numel(stat.label)
            if any(stat.mask(ch_ix,1,:))
                sig_ch{1} = [sig_ch{1} stat.label(ch_ix)];
            end
        end
        clear stat
    else    % ANOVA
        eval(['run ' root_dir 'emodim/scripts/stat_vars/' stat_id '_vars.m']);        
        f_name = [SBJ_vars.dirs.proc SBJ '_ANOVA_ROI_' stat_id '_' an_id '.mat'];
        load(f_name,'stat','w2');
        
        % FDR correct pvalues for ANOVA
        for ch_ix = 1:numel(stat.label)
            pvals = squeeze(w2.pval(:,ch_ix,:));
            [~, ~, ~, qvals] = fdr_bh(pvals);%,0.05,'pdep','yes');
            
            % Consolidate to binary sig/non-sig
            for cond_ix = 1:numel(cond_lab)
                if strcmp(cond_lab{cond_ix},'RT') && any(stat.mask(ch_ix,1,:))
                    sig_ch{cond_ix} = [sig_ch{cond_ix} {[SBJs{sbj_ix} '_' stat.label{ch_ix}]}];
                elseif any(strcmp(cond_lab{cond_ix},{'CNI','pcon'})) && any(qvals(cond_ix,:)<0.05,2)
                    sig_ch{cond_ix} = [sig_ch{cond_ix} {[SBJs{sbj_ix} '_' w2.label{ch_ix}]}];
                end
            end
        end
        clear stat w2
    end
    
    % Select elecs in atlas ROIs and hemisphere to plot
    if ~plot_out
        atlas_in_elecs = {};
        for roi_ix = 1:numel(roi_list)
            atlas_in_elecs = [atlas_in_elecs;...
                elec_sbj{sbj_ix,1}.label(strcmp(elec_sbj{sbj_ix,1}.roi,roi_list(roi_ix)))];
        end
    else
        atlas_in_elecs = elec_sbj{sbj_ix,1}.label;
    end
    if ~strcmp(hemi,'b')
        hemi_in_elecs = elec_sbj{sbj_ix,1}.label(strcmp(elec_sbj{sbj_ix,1}.hemi,hemi));
    else
        hemi_in_elecs = elec_sbj{sbj_ix,1}.label;
    end
    
    % Select sig elecs
    fprintf(sig_report,'===============================================================\n');
    for cond_ix = 1:numel(cond_lab)
        fprintf(sig_report,'\t%s - %s  = %i / %i (%.02f) sig elecs:\n',SBJs{sbj_ix},cond_lab{cond_ix},...
            numel(sig_ch{cond_ix}),numel(elec_sbj{sbj_ix,cond_ix}.label),...
            100*numel(sig_ch{cond_ix})/numel(elec_sbj{sbj_ix,cond_ix}.label));
        for sig_ix = 1:numel(sig_ch{cond_ix})
            e_ix = find(strcmp(elec_sbj{sbj_ix,cond_ix}.label,sig_ch{cond_ix}{sig_ix}));
            fprintf(sig_report,'%s - %s (%s)\n',sig_ch{cond_ix}{sig_ix},elec_sbj{sbj_ix,cond_ix}.atlas_label{e_ix},...
                elec_sbj{sbj_ix,cond_ix}.hemi{e_ix});
        end
        
        % Select sig elecs && elecs matching atlas
        % fn_select_elec messes up if you try to toss all elecs
        good_elecs = intersect(intersect(atlas_in_elecs, hemi_in_elecs), sig_ch{cond_ix});
        if numel(intersect(elec_sbj{sbj_ix,cond_ix}.label,good_elecs))==0
            elec_sbj{sbj_ix,cond_ix} = {};
            good_sbj(sbj_ix,cond_ix) = false;
            warning('WARNING!!! All sig_ch are out of atlas and/or hemisphere!');
            fprintf(sig_report,'\t!!! 0 sig_ch remain\n');
        else
            cfgs = [];
            cfgs.channel = good_elecs;
            elec_sbj{sbj_ix,cond_ix} = fn_select_elec(cfgs, elec_sbj{sbj_ix,cond_ix});
            all_roi_labels{cond_ix} = [all_roi_labels{cond_ix}; elec_sbj{sbj_ix,cond_ix}.roi];
            all_roi_colors{cond_ix} = [all_roi_colors{cond_ix}; elec_sbj{sbj_ix,cond_ix}.roi_color];
            fprintf(sig_report,'\t%i sig_ch remain\n',numel(elec_sbj{sbj_ix,cond_ix}.label));
        end
        fprintf(sig_report,'---------------------------------------------------------------\n');
    end
    fprintf(sig_report,'===============================================================\n');
    
    clear SBJ SBJ_vars SBJ_vars_cmd
end

% Combine elec structs
elec = cell([numel(cond_lab) 1]);
for cond_ix = 1:numel(cond_lab)
    elec{cond_ix} = ft_appendsens([],elec_sbj{good_sbj(:,cond_ix),cond_ix});
    elec{cond_ix}.roi       = all_roi_labels{cond_ix};    % appendsens strips that field
    elec{cond_ix}.roi_color = all_roi_colors{cond_ix};    % appendsens strips that field
end

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
f = cell(size(cond_lab));
for cond_ix = 1%:numel(cond_lab)
    if strcmp(stat_id,'actv') || strcmp(stat_id,'CSE')
        plot_name = ['GRP_' stat_id '_' an_id];
    else
        plot_name = ['GRP_ANOVA_' cond_lab{cond_ix} '_' stat_id '_' an_id];
    end
    f{cond_ix} = figure('Name',plot_name);
    
    % Plot 3D mesh
    mesh_alpha = 0.2;
    ft_plot_mesh(mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);
    
    % Plot electrodes on top
    for e = 1:numel(elec{cond_ix}.label)
        cfgs = []; cfgs.channel = elec{cond_ix}.label(e);
        elec_tmp = fn_select_elec(cfgs,elec{cond_ix});
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
    set(f{cond_ix}, 'windowkeypressfcn',   @cb_keyboard);
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
