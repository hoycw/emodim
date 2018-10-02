function elec = viewrecon(ptid)

%% VIEWRECON loads the reconstruction for a given patient. Electrodes can be
% viewed in patient space or normalized space, and as an orthoplot or in 3D.
% If the recon is unfinished or in progress, viewrecon should tell you how
% far along the recon is
%
% Use as: loadrecon(patient_ID) where patient_ID is a string
% Example: loadrecon('IR57')

% by Sandon Griffin (2017)

if any(strfind(ptid(1:2), 'OS')); hosp_code = 'Oslo';
elseif any(strfind(ptid(1:2), 'IR')); hosp_code = 'Irvine';
elseif any(strfind(ptid(1:2), 'CH')); hosp_code = 'Childrens';
elseif any(strfind(ptid(1:2), 'AM')); hosp_code = 'Albany';
elseif any(strfind(ptid(1:2), 'WV')); hosp_code = 'Wadsworth';
elseif any(strfind(ptid(1:2), 'GP')) || any(strfind(ptid(1:2), 'EC')); hosp_code = 'UCSF';
elseif any(strfind(ptid(1:2), 'SF')); hosp_code = 'SFC';
elseif any(strfind(ptid(1:2), 'CH')); hosp_code = 'Childrens';
elseif any(strfind(ptid(1:2), 'JH')); hosp_code = 'Hopkins';
elseif strfind(ptid(1:2), 'P' == 1); hosp_code = 'Germany';
elseif any(strfind(ptid(1:2), 'CP')); hosp_code = 'CPMC';
elseif any(strfind(ptid(1:2), 'NY')); hosp_code = 'New_York';
end
subdir = ['/home/knight/ecog/DATA_FOLDER/' hosp_code '/' ptid '/'];
ogrecondir = dir([subdir '3D_Images/Recon_*']);
if ~isempty(ogrecondir)
  recon_date = ogrecondir.name(end-7:end); % 3 letter month code + year
elseif ~isempty(dir([subdir '3D_Images/DICOM_Anon/CT*'])) && ~isempty(dir([subdir '3D_Images/DICOM_Anon/MR*']));
  error('No recon has been started yet for this subject, but we do have CT and MR scans so one should be started soon');
else
  error('No CT and MR scans were found');
end
if exist([subdir '3D_Images/Recon_' recon_date '/FT_Pipeline/']);
  recondir = [subdir '3D_Images/Recon_' recon_date '/FT_Pipeline/'];
else
  error('This appears to be an older case that has not been converted from the old pipeline to the fieldtrip pipeline, please try viewing the recon in BioImage Suite');
end
if ~exist([recondir '/Electrodes/' ptid '_elec_acpc.mat']) && ~exist([recondir '/Electrodes/' ptid '_elec_acpc_f.mat'])
  error('The recon has been started, but is not yet finished')
end

method = input('What do you want to see? Enter "orthoplot" or "3d"\n', 's');
if strcmp(method, 'orthoplot')
  methodortho = input('Do you want to see the electrodes in patient space or normalized space?\nEnter "patient" or "normalized"\n', 's');
  if strcmp(methodortho, 'patient')
    ct_acpc_f = ft_read_mri([recondir 'Scans/' ptid '_CT_acpc_f.nii']);
    ct_acpc_f.coordsys = 'acpc';
    % fsmri_acpc = ft_read_mri(['Scans/' patient_id '_fsMR_pre_acpc.nii']);
    if exist([recondir 'Scans/' ptid '_fsMR_post_acpc.nii']) && exist([recondir 'Scans/' ptid '_fsMR_pre_acpc.nii']);
      scantype = input('There is a post-op and pre-op MR for this subject\nEnter "post" to toggle between the post-op MR and the CT\nEnter "pre" to toggle between the pre-op MR and the post-op MR', 's');
      if strcmp(scantype, 'post')
        fsmri_post_acpc = ft_read_mri([recondir 'Scans/' ptid '_fsMR_post_acpc_f.nii']);
        fsmri_post_acpc.coordsys = 'acpc';
       
        if exist([recondir '/Electrodes/' ptid '_elec_acpc.mat']);
          load([recondir '/Electrodes/' ptid '_elec_acpc.mat']);
          cfg = [];
          cfg.elec = elec_acpc;
          ft_electrodeplacement(cfg, fsmri_post_acpc, ct_acpc_f);
        elseif exist([recondir '/Electrodes/' ptid '_elec_acpc_f.mat']);
          load([recondir 'Electrodes/' ptid '_elec_acpc_f.mat']);
          cfg = [];
          cfg.elec = elec_acpc_f;
          ft_electrodeplacement(cfg, fsmri_post_acpc, ct_acpc_f);
        end 
%         cfg = [];
%         cfg.elec = elec_acpc;
%         ft_electrodeplacement(cfg, fsmri_post_acpc, ct_acpc_f);
      elseif strcmp(scantype, 'pre')
        if exist([recondir '/Electrodes/' ptid '_elec_acpc_f.mat'])
          load([recondir '/Electrodes/' ptid '_elec_acpc_f.mat']);
        else
          error('The post-op MRI has not been snapped to the pre-op MRI yet, so electrodes can only be viewed in the post-op MRI')
        end
        fsmri_pre_acpc = ft_read_mri([recondir 'Scans/' ptid '_fsMR_pre_acpc.nii']);
        fsmri_pre_acpc.coordsys = 'acpc';
        
        fsmri_post_acpc_f = ft_read_mri([recondir 'Scans/' ptid '_fsMR_post_acpc_f.nii']);
        fsmri_post_acpc_f.coordsys = 'acpc';

        cfg = [];
        cfg.elec = elec_acpc_f;
        ft_electrodeplacement(cfg, ct_acpc_f, fsmri_post_acpc_f, fsmri_pre_acpc);
      else
        error('You must enter "post" or "pre"')
      end
    elseif exist([recondir 'Scans/' ptid '_fsMR_pre_acpc.nii']) == 2
      fsmri_pre_acpc = ft_read_mri([recondir 'Scans/' ptid '_fsMR_pre_acpc.nii']);
      fsmri_pre_acpc.coordsys = 'acpc';
      
      cfg = [];
      if exist([recondir '/Electrodes/' ptid '_elec_acpc_f.mat'])
        load([recondir '/Electrodes/' ptid '_elec_acpc_f.mat'])
        cfg.elec = elec_acpc_f;
      elseif exist([recondir '/Electrodes/' ptid '_elec_acpc.mat'])
        load([recondir '/Electrodes/' ptid '_elec_acpc.mat'])
        cfg.elec = elec_acpc;
      end
      
      ft_electrodeplacement(cfg, ct_acpc_f, fsmri_pre_acpc);
    elseif exist([recondir 'Scans/' ptid '_fsMR_post_acpc.nii']) == 2
      fsmri_post_acpc = ft_read_mri([recondir 'Scans/' ptid '_fsMR_post_acpc.nii']);
      fsmri_post_acpc.coordsys = 'acpc';
      
      cfg = [];
      if exist([recondir '/Electrodes/' ptid '_elec_acpc_f.mat'])
        load([recondir '/Electrodes/' ptid '_elec_acpc_f.mat'])
        cfg.elec = elec_acpc_f;
      elseif exist([recondir '/Electrodes/' ptid '_elec_acpc.mat'])
        load([recondir '/Electrodes/' ptid '_elec_acpc.mat'])
        cfg.elec = elec_acpc;
      end
      ft_electrodeplacement(cfg, ct_acpc_f, fsmri_post_acpc);
    end
  elseif strcmp(methodortho, 'normalized')
    mri_template = ft_read_mri('/home/knight/smg211/fieldtrip/template/anatomy/single_subj_T1_1mm.nii');
    mri_template.coordsys = 'mni';
    cfg = [];
    
    if exist([recondir 'Electrodes/' ptid '_elec_mni_v.mat'])
      load([recondir 'Electrodes/' ptid '_elec_mni_v.mat'])
      cfg.elec = elec_mni_v;
    elseif exist([recondir 'Electrodes/' ptid '_elec_mni_frv.mat'])
      load([recondir 'Electrodes/' ptid '_elec_mni_frv.mat'])
      cfg.elec = elec_mni_frv;
    else
      error('The recon has been started but the electrodes have not been volume-based normalized yet')
    end
    
    ft_electrodeplacement(cfg, mri_template);
  else
    error('you must enter "patient" or "normalized"')
  end
elseif strcmp(method, '3d')
  warning('To adjust specific aspects of the 3d figure, such as transparency or adding electrode labels, please explicitly use ft_plot_mesh (to plot the surface) and ft_plot_sens (to plot the electrodes) and refer to the documentation for those functions.')
  h = figure;
  method3d = input('Do you want to see the electrodes in patient space or normalized space?\nEnter "patient" or "normalized"\n', 's');
  if strcmp(method3d, 'patient')
    load([recondir 'Surfaces/' ptid '_cortex_lh.mat']);
    load([recondir 'Surfaces/' ptid '_cortex_rh.mat']);
    %view([-55 10]); 
    view([0 90]); 
    
%     h = camlight;
    methodlabels = input('Do you want to see the electrodes with our without labels?\nEnter "labels" or "no labels"\n', 's');
    if strcmp(methodlabels, 'no labels')
      if exist([recondir 'Electrodes/' ptid '_elec_acpc_fr.mat'])
        ft_plot_mesh(cortex_rh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none');
        ft_plot_mesh(cortex_lh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none');
        load([recondir 'Electrodes/' ptid '_elec_acpc_fr.mat'])
        ft_plot_sens(elec_acpc_fr, 'elecshape', 'sphere');
      elseif exist([recondir 'Electrodes/' ptid '_elec_acpc_r.mat'])
        ft_plot_mesh(cortex_rh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none');
        ft_plot_mesh(cortex_lh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none');
        load([recondir 'Electrodes/' ptid '_elec_acpc_r.mat'])
        ft_plot_sens(elec_acpc_r, 'elecshape', 'sphere');
      elseif exist([recondir 'Electrodes/' ptid '_elec_acpc.mat'])
        load([recondir 'Electrodes/' ptid '_elec_acpc.mat'])
        ft_plot_mesh(cortex_rh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', 0.3);
        ft_plot_mesh(cortex_lh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', 0.3);
        ft_plot_sens(elec_acpc, 'elecshape', 'sphere');
        warning('Assuming that you are trying to view depth electrodes in a 3d plot, the surface mesh has been made transparent to allow viewing electrodes beneath the surface. For a more accurate visualization of the location of depth electrodes, use orthoplot. If this case contains surface electrodes, they have not been projected to the surface yet, so you may not be able to see them properly in this figure.');
      else
        error('The recon has been started but the electrodes have not been projected to the surface yet or there are no surface electrodes')
      end
    elseif strcmp(methodlabels, 'labels')
      if exist([recondir 'Electrodes/' ptid '_elec_acpc_fr.mat'])
        ft_plot_mesh(cortex_rh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none');
        ft_plot_mesh(cortex_lh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none');
        load([recondir 'Electrodes/' ptid '_elec_acpc_fr.mat'])
        ft_plot_sens(elec_acpc_fr, 'elecshape', 'sphere','label', 'label');
      elseif exist([recondir 'Electrodes/' ptid '_elec_acpc_r.mat'])
        ft_plot_mesh(cortex_rh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none');
        ft_plot_mesh(cortex_lh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none');
        load([recondir 'Electrodes/' ptid '_elec_acpc_r.mat'])
        ft_plot_sens(elec_acpc_r, 'elecshape', 'sphere','label', 'label');
      elseif exist([recondir 'Electrodes/' ptid '_elec_acpc.mat'])
        load([recondir 'Electrodes/' ptid '_elec_acpc.mat'])
        ft_plot_mesh(cortex_rh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', 0.3);
        ft_plot_mesh(cortex_lh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', 0.3);
        ft_plot_sens(elec_acpc, 'elecshape', 'sphere','label', 'label');
        warning('Assuming that you are trying to view depth electrodes in a 3d plot, the surface mesh has been made transparent to allow viewing electrodes beneath the surface. For a more accurate visualization of the location of depth electrodes, use orthoplot. If this case contains surface electrodes, they have not been projected to the surface yet, so you may not be able to see them properly in this figure.');
      else
        error('The recon has been started but the electrodes have not been projected to the surface yet or there are no surface electrodes')
      end
    end
  elseif strcmp(method3d, 'normalized')
    method3dnorm = input('Do you want to see the volume normalized-based or surface-based normalized electrodes?\nEnter "volume" or "surface"\n', 's');
    if strcmp(method3dnorm, 'volume')
      surface_template_l = load('/home/knight/smg211/fieldtrip/template/anatomy/surface_pial_left.mat');
      surface_template_r = load('/home/knight/smg211/fieldtrip/template/anatomy/surface_pial_right.mat');
      ft_plot_mesh(surface_template_l.mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none');
      ft_plot_mesh(surface_template_r.mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none');
      view([0 90]); lighting gouraud; camlight;
      
      methodlabels = input('Do you want to see the electrodes with our without labels?\nEnter "labels" or "no labels"\n', 's');
      if strcmp(methodlabels, 'no labels')
        if exist([recondir 'Electrodes/' ptid '_elec_mni_v.mat'])
          load([recondir 'Electrodes/' ptid '_elec_mni_v.mat'])
          ft_plot_sens(elec_mni_v, 'elecshape', 'sphere');
        elseif exist([recondir 'Electrodes/' ptid '_elec_mni_frv.mat'])
          load([recondir 'Electrodes/' ptid '_elec_mni_frv.mat'])
          ft_plot_sens(elec_mni_frv, 'elecshape', 'sphere');
        else
          error('The recon has been started but the electrodes have not been volume-based normalized yet')
        end
        
      elseif strcmp(methodlabels, 'labels')
        if exist([recondir 'Electrodes/' ptid '_elec_mni_v.mat'])
          load([recondir 'Electrodes/' ptid '_elec_mni_v.mat'])
          ft_plot_sens(elec_mni_v, 'elecshape', 'sphere','label', 'label');
        elseif exist([recondir 'Electrodes/' ptid '_elec_mni_frv.mat'])
          load([recondir 'Electrodes/' ptid '_elec_mni_frv.mat'])
          ft_plot_sens(elec_mni_frv, 'elecshape', 'sphere','label', 'label');
        else
          error('The recon has been started but the electrodes have not been volume-based normalized yet')
        end
      end
      
      
    elseif strcmp(method3dnorm, 'surface')
      if ~exist([recondir 'Electrodes/' ptid '_elec_mni_s.mat']) && ~exist([recondir 'Electrodes/' ptid '_elec_fsavg_frs.mat'])
        error('The recon has been started but the electrodes have either not been surface-based normalized or there are only depth, in which case surface-based normalization is not possible')
      end
      
      fshome = '/usr/local/freesurfer_x86_64-5.3.0';
      fs_surface_template_l = ft_read_headshape([fshome '/subjects/fsaverage/surf/lh.pial']);
      fs_surface_template_l.coordsys = 'tal';
      fs_surface_template_r = ft_read_headshape([fshome '/subjects/fsaverage/surf/rh.pial']);
      fs_surface_template_r.coordsys = 'tal';
      ft_plot_mesh(fs_surface_template_l, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none');
      ft_plot_mesh(fs_surface_template_r, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none');
      view([0 90]); lighting gouraud; camlight;
      
      methodlabels = input('Do you want to see the electrodes with our without labels?\nEnter "labels" or "no labels"\n', 's');
      if strcmp(methodlabels, 'no labels')
        if exist([recondir 'Electrodes/' ptid '_elec_mni_s.mat'])
          load([recondir 'Electrodes/' ptid '_elec_mni_s.mat'])
          ft_plot_sens(elec_mni_s, 'elecshape', 'sphere');
        elseif exist([recondir 'Electrodes/' ptid '_elec_fsavg_frs.mat'])
          load([recondir 'Electrodes/' ptid '_elec_fsavg_frs.mat'])
          ft_plot_sens(elec_fsavg_frs, 'elecshape', 'sphere');
        else
          error('The recon has been started but the electrodes have either not been surface-based normalized or there are only depth, in which case surface-based normalization is not possible')
        end
      elseif strcmp(methodlabels, 'labels')
        if exist([recondir 'Electrodes/' ptid '_elec_mni_s.mat'])
          load([recondir 'Electrodes/' ptid '_elec_mni_s.mat'])
          ft_plot_sens(elec_mni_s, 'elecshape', 'sphere','label', 'label');
        elseif exist([recondir 'Electrodes/' ptid '_elec_fsavg_frs.mat'])
          load([recondir 'Electrodes/' ptid '_elec_fsavg_frs.mat'])
          ft_plot_sens(elec_fsavg_frs, 'elecshape', 'sphere','label', 'label');
        else
          error('The recon has been started but the electrodes have either not been surface-based normalized or there are only depth, in which case surface-based normalization is not possible')
        end
      end
    else
      error('you must enter "volume" or "surface"')
    end
  else
    error('you must enter "patient" or "normalized"')
  end
  
  l = camlight; lighting gouraud; material dull;
  fprintf(['To reset the position of the camera light after rotating the figure,\n' ...
    'make sure none of the figure adjustment tools (e.g., zoom, rotate) are active\n' ...
    '(i.e., uncheck them within the figure), and then hit ''l'' on the keyboard\n'])
  set(h, 'windowkeypressfcn',   @cb_keyboard);
else
  error('you must enter "orthoplot" or "3d"')
end

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
