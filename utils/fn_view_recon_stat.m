function fn_view_recon_stat(SBJ, pipeline_id, view_space, reg_type, show_labels, hemi, plot_out, thresh, sig_method, coloring)%stat_id, an_id, 
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
%   plot_out [0/1] - exclude electrodes that don't match atlas or aren't in hemisphere
%   thresh [float] - alpha value (0.05)
%   sig_method [str] - 'prop_boot' or 'norminv'
%   coloring [str] - type of coloring to use
%       'cnts': a continuous inferno colormap for plotting R2 individual
%           (one brain per model)
%       'mcmp': model comparison, using sig across 3 models as colormap

[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
addpath([root_dir 'emodim/scripts/Colormaps/']);
SBJ_vars_cmd = ['run ' root_dir 'emodim/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

stat_id  = {'OM','RM','RT'};
stat_lab = {'OrigMean', 'ReginaMean', 'ReginaTimeAnn'};

ns_color = [0 0 0];
view_angle = [-90 0];
if strcmp(reg_type,'v') || strcmp(reg_type,'s')
    reg_suffix = ['_' reg_type];
else
    reg_suffix = '';
end

%% Load elec struct
load([SBJ_vars.dirs.recon,SBJ,'_elec_',pipeline_id,'_',view_space,reg_suffix,'_Dx_gROI.mat']);

%% Remove electrodes that aren't in atlas ROIs
if ~plot_out
%     atlas_out_elecs = elec.label(strcmp(elec.atlas_label,'no_label_found'));
    if ~strcmp(hemi,'b')
        hemi_out_elecs = elec.label(~strcmp(elec.hemi,hemi));
    else
        hemi_out_elecs = {};
    end
    cfgs = []; cfgs.channel = [{'all'} fn_ch_lab_negate(hemi_out_elecs)];%fn_ch_lab_negate(atlas_out_elecs) 
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

%% Load and process stats
% Load stats
load([root_dir 'emodim/data/predCorr.mat']);

% Check elec match
if numel(elec.label)~=numel(predCorr.(SBJ).label) || ~isempty(setdiff(elec.label,predCorr.(SBJ).label))
    warning(['WARNING!!! Mismatch in electrodes in stat (n=' ...
        num2str(numel(predCorr.(SBJ).label)) ') and elec (n=' num2str(numel(elec.label)) ')!']);
    cfgs = []; cfgs.channel = predCorr.(SBJ).label;
    elec = fn_select_elec(cfgs,elec);
%     error('Mismatch in electrodes in stat and elec!');
end

% Get data and threshold
cors = cell(size(stat_lab));
boot = cell(size(stat_lab));
for st = 1:numel(stat_id)
    cors{st}  = predCorr.(SBJ).(stat_lab{st}).Correlations;
    boot{st}  = predCorr.(SBJ).(stat_lab{st}).Bootstrap;
end

% Find significant values
sig  = zeros([numel(elec.label) numel(stat_lab)]);
stat = zeros([numel(elec.label) numel(stat_lab)]);
% cut = zeros([numel(elec.label) numel(stat_lab) 2]); % Lower, Upper
max_r = -1;
min_r = 1;
% figure('Name',[SBJ '_' sig_method]);
for st = 1:numel(stat_lab)
    % Get limits of colormap
    if max(cors{st})>max_r
        max_r = max(cors{st});
    end
    if min(cors{st})<min_r
        min_r = min(cors{st});
    end
    % Check significance
    if thresh
        fprintf('%s Significant Electrodes:\n\t',stat_lab{st});
        for e = 1:numel(elec.label)
            if strcmp(sig_method,'prop_boot')
                % Proportion of bootstrap values above zero
                stat(e,st) = 1-sum(boot{st}(:,e)>0)/size(boot{st},1);
                if stat(e,st) <= thresh
                    sig(e,st) = 1;
                    fprintf('%s\t',elec.label{e});
                end
            elseif strcmp(sig_method,'norminv')
                % STDs of bootstrap above 0
                sd = std(boot{st}(:,e),0,1);
                stat(e,st) = -norminv(thresh)*sd;
                if cors{st}(e) >= stat(e,st)
                    sig(e,st) = 1;
                    fprintf('%s\t',elec.label{e});
                end
            end
            % Erroneous (thought it was permutation distribution)
%             boot_sort = sort(boot{st}(:,e));
%             boot_cut_ix = int32(round(size(boot{st},1)*thresh/2));
%             cut(e,st,:) = [boot_sort(boot_cut_ix) boot_sort(end-boot_cut_ix)];
%             if cors{st}(e) <= cut(e,st,1) || cors{st}(e) >= cut(e,st,2)
%                 sig(e,st) = 1;
%             end
        end
        fprintf('\nTotal n_sig = %i / %i\n',sum(sig(:,st)),numel(cors{st}));
    end
    % Plot the correlations vs. threshold
%     subplot(numel(stat_lab),1,st);
%     if strcmp(sig_method,'norminv')
%         scatter(1:numel(elec.label),cors{st},'k'); hold on;
%         scatter(1:numel(elec.label),stat(:,st),'r');
%         legend('r','SD thresh');
%     elseif strcmp(sig_method,'prop_boot')
%         plot(stat(:,st));hold on; line(xlim,[thresh thresh],'Color','r');
%         ylabel('pval');
%     end
%     xlabel('elec#');
% %     plot(cors{st},'k');hold on;plot(squeeze(cut(:,st,:)),'r');
% %     legend('Corr','Sig Thresh');
%     title([stat_lab{st} ':' num2str(sum(sig(:,st))) 'sig']);
%     figure;
%     for e = 1:numel(elec.label);
%         histogram(boot{st}(:,e),30);hold on;
%         line([cors{st}(e) cors{st}(e)],ylim,'Color','r');
%         title([stat_lab{st} ':' elec.label{e} ',sig=' num2str(sig(e,st)) ',' num2str(stat(e,st))]);
%         pause;hold off;
%     end
end
out_fname = [SBJ_vars.dirs.proc SBJ '_model_fit_corr_stats_' sig_method '_' num2str(thresh) '.mat'];
save(out_fname,'-v7.3','sig','stat','cors','boot');

% Map to colors
if strcmp(coloring,'cnts')
    cmap = colormap(inferno());% viridis()
    % cmap = colormap('cool');
    cmap_idx = linspace(min_r,max_r,size(cmap,1));
    
    elec_colors = cell(size(stat_lab));
    for st = 1:numel(stat_lab)
        elec_colors{st} = zeros([numel(elec.label) 3]);
        for e = 1:numel(elec.label)
            r_match = find(cmap_idx<=cors{st}(e));
            elec_colors{st}(e,:) = cmap(r_match(end),:);
        end
    end
elseif strcmp(coloring,'mcmp')
    elec_colors = sig./1.5;
end

%% 3D Surface + Grids (3d, pat/mni, vol/srf, 0/1)
if strcmp(coloring,'cnts')
    f = {};
    for st = 1:numel(stat_lab)
        plot_name = [SBJ '_HFA_' stat_id{st} '_sigCorr'];
        f{st} = figure('Name',plot_name);
        
        % Plot 3D mesh
        mesh_alpha = 0.8;
        if any(strcmp(SBJ_vars.ch_lab.probe_type,'seeg'))
            mesh_alpha = 0.2;
        end
        ft_plot_mesh(mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);
        
        % Plot electrodes on top
        cfgs = [];
        for e = 1:numel(elec.label)
            if thresh==0 || sig(e,st)
                cfgs.channel = elec.label{e};
                elec_tmp = fn_select_elec(cfgs, elec);
                if show_labels
                    lab_arg = 'label';
                else
                    lab_arg = 'off';
                end
                ft_plot_sens(elec_tmp, 'elecshape', 'sphere', 'facecolor', elec_colors{st}(e,:), 'label', lab_arg);
            end
        end
        colorbar;
        colormap(cmap);
        caxis([min_r max_r]);
        
        view(view_angle); material dull; lighting gouraud;
        l = camlight;
        fprintf(['To reset the position of the camera light after rotating the figure,\n' ...
            'make sure none of the figure adjustment tools (e.g., zoom, rotate) are active\n' ...
            '(i.e., uncheck them within the figure), and then hit ''l'' on the keyboard\n'])
        set(f{st}, 'windowkeypressfcn',   @cb_keyboard);
    end
elseif strcmp(coloring,'mcmp')
    plot_name = [SBJ '_HFA_' stat_id{:} '_sigCorr'];
    f = figure('Name',plot_name);
    
    % Plot 3D mesh
    mesh_alpha = 0.8;
    if any(strcmp(SBJ_vars.ch_lab.probe_type,'seeg'))
        mesh_alpha = 0.2;
    end
    ft_plot_mesh(mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);
    
    % Plot electrodes on top
    cfgs = [];
    for e = 1:numel(elec.label)
        if thresh==0 || any(sig(e,:))
            cfgs.channel = elec.label{e};
            elec_tmp = fn_select_elec(cfgs, elec);
            if show_labels
                lab_arg = 'label';
            else
                lab_arg = 'off';
            end
            ft_plot_sens(elec_tmp, 'elecshape', 'sphere', 'facecolor', elec_colors(e,:), 'label', lab_arg);
        end
    end
    
    view(view_angle); material dull; lighting gouraud;
    l = camlight;
    fprintf(['To reset the position of the camera light after rotating the figure,\n' ...
        'make sure none of the figure adjustment tools (e.g., zoom, rotate) are active\n' ...
        '(i.e., uncheck them within the figure), and then hit ''l'' on the keyboard\n'])
    set(f, 'windowkeypressfcn',   @cb_keyboard);
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
