function fn_view_recon_atlas_grp_stat(SBJs, proc_id, stat_id, thresh, sig_method, reg_type, show_labels,...
                                 hemi, atlas_id, roi_id, mirror, plot_out, varargin)
%% Plot a reconstruction with electrodes
% INPUTS:
%   SBJs [cell array str] - subject IDs to plot
%   proc_id [str] - name of analysis pipeline, used to pick elec file
%   stat_id [str] - which set of results to load
%   thresh [float] - alpha value (0.05)
%   sig_method [str] - 'prop_boot' or 'norminv'
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
addpath([root_dir 'emodim/scripts/Colormaps/']);

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

%% Load and process stats
load([root_dir 'emodim/data/predCorr_' stat_id '.mat']);
if strcmp(stat_id,'EpocEx'); predCorr = predCorrEpocEx; end
stat_lab = setdiff(fieldnames(predCorr.(SBJs{1})),{'label'});

%% Load elec struct
cors = cell([numel(SBJs) numel(stat_lab)]);
boot = cell([numel(SBJs) numel(stat_lab)]);
elec_sbj = cell([numel(SBJs) numel(stat_lab)]);
good_sbj = true([numel(SBJs) numel(stat_lab)]);
all_roi_labels = cell([numel(stat_lab) 1]);
all_roi_colors = cell([numel(stat_lab) 1]);
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    SBJ_vars_cmd = ['run ' root_dir 'emodim/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Load MNI elecs
    elec_mni_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_mni',reg_suffix,'.mat'];
    tmp = load(elec_mni_fname); elec_sbj{sbj_ix,1} = tmp.elec;
    
    % Load manually adjusted gROI labels --> !!! fix this hack once I have elec_final!
    if ~strcmp(roi_id,'gROI'); error('only doing gROI now!'); end
    elec_gROI_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat_Dx_gROI.mat'];
    tmp = load(elec_gROI_fname); elec_gROI = tmp.elec;
%     if strcmp(SBJ,'IR66')
% %         error('undo this hack!');
%         elec_gROI.hemi = [elec_gROI.hemi(1:6); {'l';'l'}; elec_gROI.hemi(7:51);...
%                             {'r';'r'}; elec_gROI.hemi(52:end)];
%         elec_gROI.roi  = [elec_gROI.roi(1:6); {'AMG';'AMG'}; elec_gROI.roi(7:51);...
%                             {'AMG';'AMG'}; elec_gROI.roi(52:end)];
    if numel(elec_sbj{sbj_ix,1}.label)~=numel(elec_gROI.label) ||...
           ~all(strcmp(elec_sbj{sbj_ix,1}.label,elec_gROI.label))
        error('pat_gROI and mni elec labels mismatch!');
    end
    elec_sbj{sbj_ix,1}.(roi_field) = elec_gROI.roi;
    if any(~strcmp(elec_sbj{sbj_ix,1}.hemi,elec_gROI.hemi))
        warning([SBJ ' ' num2str(sum(~strcmp(elec_sbj{sbj_ix,1}.hemi,elec_gROI.hemi))) ...
            ' different hemi on mni vs. gROI elecs!']);
    end
    elec_sbj{sbj_ix,1}.hemi = elec_gROI.hemi;
    elec_sbj{sbj_ix,1}.color = fn_roi2color(elec_gROI.roi);
    
    % Check elec match
    if numel(elec_sbj{sbj_ix,1}.label)~=numel(predCorr.(SBJ).label) || ...
            ~isempty(setdiff(elec_sbj{sbj_ix,1}.label,predCorr.(SBJ).label))
        warning(['WARNING!!! Mismatch in electrodes in stat (n=' ...
            num2str(numel(predCorr.(SBJ).label)) ') and elec (n=' num2str(numel(elec_sbj{sbj_ix,1}.label)) ')!']);
        cfgs = []; cfgs.channel = predCorr.(SBJ).label;
        elec_sbj{sbj_ix,1} = fn_select_elec(cfgs,elec_sbj{sbj_ix,1});
        elec_sbj{sbj_ix,1} = fn_reorder_elec(elec_sbj{sbj_ix,1},predCorr.(SBJ).label);
        %     error('Mismatch in electrodes in stat and elec!');
    end
    
    % Append SBJ name to labels
    for e_ix = 1:numel(elec_sbj{sbj_ix,1}.label)
        elec_sbj{sbj_ix,1}.label{e_ix} = [SBJs{sbj_ix} '_' elec_sbj{sbj_ix,1}.label{e_ix}];
    end
    
    if ~plot_out
        % Remove electrodes that aren't in atlas ROIs & hemisphere
        if mirror
            good_elecs = fn_select_elec_lab_match(elec_sbj{sbj_ix,1}, 'b', atlas_id, roi_id);
        else
            good_elecs = fn_select_elec_lab_match(elec_sbj{sbj_ix,1}, hemi, atlas_id, roi_id);
        end
    else
        % Remove electrodes that aren't in hemisphere
        if mirror
            good_elecs = fn_select_elec_lab_match(elec_sbj{sbj_ix,1}, 'b', [], []);
        else
            good_elecs = fn_select_elec_lab_match(elec_sbj{sbj_ix,1}, hemi, [], []);
        end
    end
    
    % Mirror hemispheres
    if mirror
        elec_sbj{sbj_ix,1}.chanpos(~strcmp(elec_sbj{sbj_ix,1}.hemi,hemi),1) = ...
                        -elec_sbj{sbj_ix,1}.chanpos(~strcmp(elec_sbj{sbj_ix,1}.hemi,hemi),1);
        hemi_str = [hemi 'b'];
    else
        hemi_str = hemi;
    end
    
    % Copy for other conditions
    for st_ix = 2:numel(stat_lab)
        elec_sbj{sbj_ix,st_ix} = elec_sbj{sbj_ix,1};
    end
    
    % Select sig elecs
    for st_ix = 1:numel(stat_lab)
        % Get data
        cors{sbj_ix,st_ix}  = predCorr.(SBJ).(stat_lab{st_ix}).Correlations;
        boot{sbj_ix,st_ix}  = predCorr.(SBJ).(stat_lab{st_ix}).Bootstrap;
        
        % Threshold data
        sig  = false(size(elec_sbj{sbj_ix,st_ix}.label));
        stat = zeros(size(elec_sbj{sbj_ix,st_ix}.label));
        fprintf('%s Significant Electrodes:\n\t',stat_lab{st_ix});
        for e = 1:numel(elec_sbj{sbj_ix,st_ix}.label)
            if strcmp(sig_method,'prop_boot')
                % Proportion of bootstrap values above zero
                stat(e) = 1-sum(boot{sbj_ix,st_ix}(:,e)>0)/size(boot{sbj_ix,st_ix},1);
                if stat(e) <= thresh
                    sig(e) = true;
                    fprintf('%s\t',elec_sbj{sbj_ix,st_ix}.label{e});
                end
            elseif strcmp(sig_method,'norminv')
                error('dont use this until sure its appropriate, see email Slides ready to edit on 12/11/18');
                % STDs of bootstrap above 0
                sd = std(boot{sbj_ix,st_ix}(:,e),0,1);
                stat(e) = -norminv(thresh)*sd;
                if cors{sbj_ix,st_ix}(e) >= stat(e)
                    sig(e) = true;
                    fprintf('%s\t',elec_sbj{sbj_ix,st_ix}.label{e});
                end
            end
        end
        fprintf('\nTotal n_sig = %i / %i\n',sum(sig),numel(cors{sbj_ix,st_ix}));
        
        % Select sig elecs && elecs matching atlas
        % fn_select_elec messes up if you try to toss all elecs
        plot_elecs = intersect(good_elecs, elec_sbj{sbj_ix,st_ix}.label(sig));
        fprintf('\nTotal plotting: %i\n',numel(plot_elecs));
        if isempty(plot_elecs)
            elec_sbj{sbj_ix,st_ix} = {};
            good_sbj(sbj_ix,st_ix) = false;
            warning('WARNING!!! All sig_ch are out of atlas and/or hemisphere!');
        else
            cfgs = [];
            cfgs.channel = plot_elecs;
            elec_sbj{sbj_ix,st_ix} = fn_select_elec(cfgs, elec_sbj{sbj_ix,st_ix});
            all_roi_labels{st_ix} = [all_roi_labels{st_ix}; elec_sbj{sbj_ix,st_ix}.(roi_field)];
            all_roi_colors{st_ix} = [all_roi_colors{st_ix}; elec_sbj{sbj_ix,st_ix}.color];
        end
    end
    clear SBJ SBJ_vars SBJ_vars_cmd
end

%% Combine elec structs
elec = cell([numel(stat_lab) 1]);
for st_ix = 1:numel(stat_lab)
    elec{st_ix} = ft_appendsens([],elec_sbj{good_sbj(:,st_ix),st_ix});
    elec{st_ix}.roi   = all_roi_labels{st_ix};    % appendsens strips that field
    elec{st_ix}.color = all_roi_colors{st_ix};    % appendsens strips that field
end

%% Load brain recon
mesh = fn_load_recon_mesh([],'mni',reg_type,mesh_type,hemi);

%% 3D Surface + Grids (3d, pat/mni, vol/srf, 0/1)
out_dir = [root_dir 'emodim/results/HFA/GRP_recons/' stat_id '/' sig_method '_' ...
    num2str(thresh) '/' atlas_id '_' roi_id '/'];
if ~exist(out_dir,'dir')
    [~,~] = mkdir(out_dir);
end
% f = cell(size(stat_lab));
% for stat_ix = 1:numel(stat_lab)
%     plot_name = ['GRP_' stat_lab{stat_ix} '_sig_' atlas_id '_' roi_id];
%     f{stat_ix} = figure('Name',plot_name);
plot_name = ['GRP_' stat_id '_allModels_' sig_method '_' num2str(thresh) '_' ...
                atlas_id '_' roi_id '_' hemi_str '_' view_str];
f = figure('Name',plot_name);

% Plot 3D mesh
ft_plot_mesh(mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);

% Plot electrodes on top
plotted_elecs = {};
for stat_ix = 1:numel(stat_lab)
    for e = 1:numel(elec{stat_ix}.label)
        if ~any(strcmp(plotted_elecs,elec{stat_ix}.label(e)))
            plotted_elecs = [plotted_elecs; elec{stat_ix}.label(e)];
            cfgs = []; cfgs.channel = elec{stat_ix}.label(e);
            elec_tmp = fn_select_elec(cfgs,elec{stat_ix});
            ft_plot_sens(elec_tmp, 'elecshape', 'sphere', 'facecolor', elec_tmp.color, 'label', lab_arg);
        end
    end
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
