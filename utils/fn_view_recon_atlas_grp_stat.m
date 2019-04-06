function fn_view_recon_atlas_grp_stat(SBJs, pipeline_id, reg_type, show_labels, hemi,...
                                 atlas_id, roi_id, plot_out, thresh, sig_method, varargin)%, coloring)
%% Plot a reconstruction with electrodes
% INPUTS:
%   SBJs [cell array str] - subject IDs to plot
%   pipeline_id [str] - name of analysis pipeline, used to pick elec file
%   stat_id [str] - ID of the stats
%   an_id [str] - analysis ID for preprocessing, filtering, etc.
%   reg_type [str] - {'v', 's'} choose volume-based or surface-based registration
%   show_labels [0/1] - plot the electrode labels
%   hemi [str] - {'l', 'r', 'b'} hemisphere to plot
%   atlas_id [str] - {'DK','Dx','Yeo7','Yeo17'}
%   roi_id [str] - ROI grouping by which to color the atlas ROIs
%       'gROI','mgROI','main3' - general ROIs (lobes or broad regions)
%       'ROI','thryROI','LPFC','MPFC','OFC','INS' - specific ROIs (within these larger regions)
%       'Yeo7','Yeo17' - colored by Yeo networks
%       'tissue','tissueC' - colored by tisseu compartment, e.g., GM vs WM vs OUT
%   plot_out [0/1] - include electrodes that don't have an atlas label or in hemi?
%   thresh [float] - alpha value (0.05)
%   sig_method [str] - 'prop_boot' or 'norminv'

% NOT USING THIS (YET), only gROI coloring
%   coloring [str] - type of coloring to use
%       'cnts': a continuous inferno colormap for plotting R2 individual
%           (one brain per model)
%       'mcmp': model comparison, using sig across 3 models as colormap

[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
addpath([root_dir 'emodim/scripts/Colormaps/']);

%% Handle variables
stat_id  = {'OM','RM','RT'};
stat_lab = {'OrigMean', 'ReginaMean', 'ReginaTimeAnn'};

% Handle variable inputs
if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'view_angle')
            view_angle = varargin{v+1};
        elseif strcmp(varargin{v},'mesh_alpha') && varargin{v+1}>0 && varargin{v+1}<=1
            mesh_alpha = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

% Define default options
% view_space = 'mni';
if ~exist('view_angle','var')
    if strcmp(hemi,'l')
        view_angle = [-90 0];
    elseif any(strcmp(hemi,{'r','b'}))
        view_angle = [90 0];
    else
        error(['unknown hemi: ' hemi]);
    end
end
if ~exist('mesh_alpha','var')
    % assume SEEG
    mesh_alpha = 0.3;
end
if show_labels
    lab_arg = 'label';
else
    lab_arg = 'off';
end
if strcmp(reg_type,'v') || strcmp(reg_type,'s')
    reg_suffix = ['_' reg_type];    % MNI space
else
    reg_suffix = '';                % Patient space
end

%% Load and process stats
load([root_dir 'emodim/data/predCorrEpocEx.mat']);

%% Load Elec
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
    
    % Load elec struct
    mni_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',pipeline_id,'_mni',reg_suffix,'.mat'];
    tmp = load(mni_fname); elec_sbj{sbj_ix,1} = tmp.elec;
    roi_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',pipeline_id,'_pat_' atlas_id '_' roi_id '.mat'];
    tmp = load(roi_fname);
    if strcmp(SBJ,'IR66')
        error('undo this hack!');
        tmp.elec.hemi = [tmp.elec.hemi(1:6); {'l';'l'}; tmp.elec.hemi(7:51); {'r';'r'}; tmp.elec.hemi(52:end)];
        tmp.elec.roi  = [tmp.elec.roi(1:6); {'AMG';'AMG'}; tmp.elec.roi(7:51); {'AMG';'AMG'}; tmp.elec.roi(52:end)];
    elseif ~all(strcmp(elec_sbj{sbj_ix,1}.label,tmp.elec.label))
        error('pat_roi and mni elec labels mismatch!');
    end
    elec_sbj{sbj_ix,1}.hemi = tmp.elec.hemi;
    elec_sbj{sbj_ix,1}.roi  = tmp.elec.roi;
%     try
%         elec_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',pipeline_id,'_mni',reg_suffix,'_',atlas_id,'_full.mat'];
%         if exist([elec_fname(1:end-4) '_' roi_id '.mat'],'file')
%             elec_fname = [elec_fname(1:end-4) '_' roi_id '.mat'];
%         end
%         tmp = load(elec_fname); elec_sbj{sbj_ix} = tmp.elec;
%     catch
%         error([elec_fname 'doesnt exist, exiting...']);
%     end
    
    % Check elec match
    if numel(elec_sbj{sbj_ix,1}.label)~=numel(predCorrEpocEx.(SBJ).label) || ...
            ~isempty(setdiff(elec_sbj{sbj_ix,1}.label,predCorrEpocEx.(SBJ).label))
        warning(['WARNING!!! Mismatch in electrodes in stat (n=' ...
            num2str(numel(predCorrEpocEx.(SBJ).label)) ') and elec (n=' num2str(numel(elec_sbj{sbj_ix,1}.label)) ')!']);
        cfgs = []; cfgs.channel = predCorrEpocEx.(SBJ).label;
        elec_sbj{sbj_ix,1} = fn_select_elec(cfgs,elec_sbj{sbj_ix,1});
        %     error('Mismatch in electrodes in stat and elec!');
    end
    
    % Append SBJ name to labels
    for e_ix = 1:numel(elec_sbj{sbj_ix,1}.label)
        elec_sbj{sbj_ix,1}.label{e_ix} = [SBJs{sbj_ix} '_' elec_sbj{sbj_ix,1}.label{e_ix}];
    end
    
    % Match elecs to atlas ROIs
    if any(strcmp(atlas_id,{'DK','Dx','Yeo7'}))
%         if ~isfield(elec_sbj{sbj_ix,1},'man_adj')
%             elec_sbj{sbj_ix,1}.roi       = fn_atlas2roi_labels(elec_sbj{sbj_ix,1}.atlas_lab,atlas_id,roi_id);
%         end
        if strcmp(roi_id,'tissueC')
            elec_sbj{sbj_ix,1}.roi_color = fn_tissue2color(elec_sbj{sbj_ix,1});
        elseif strcmp(atlas_id,'Yeo7')
            elec_sbj{sbj_ix,1}.roi_color = fn_atlas2color(atlas_id,elec_sbj{sbj_ix,1}.roi);
        else
            elec_sbj{sbj_ix,1}.roi_color = fn_roi2color(elec_sbj{sbj_ix,1}.roi);
        end
    elseif any(strcmp(atlas_id,{'Yeo17'}))
%         if ~isfield(elec_sbj{sbj_ix,1},'man_adj')
%             elec_sbj{sbj_ix,1}.roi       = elec_sbj{sbj_ix,1}.atlas_lab;
%         end
        elec_sbj{sbj_ix,1}.roi_color = fn_atlas2color(atlas_id,elec_sbj{sbj_ix,1}.roi);
    end
    
    % Copy for other conditions
    for stat_ix = 2:numel(stat_lab)
        elec_sbj{sbj_ix,stat_ix} = elec_sbj{sbj_ix,1};
    end
    
    % Remove hemi and/or atlas elecs
    if ~plot_out
        % Remove electrodes that aren't in atlas ROIs & hemisphere
        good_elecs = fn_select_elec_lab_match(elec_sbj{sbj_ix}, hemi, atlas_id, roi_id);
    else
        % Remove electrodes that aren't in hemisphere
        good_elecs = fn_select_elec_lab_match(elec_sbj{sbj_ix}, hemi, [], []);
    end
    
    % Select sig elecs
    for stat_ix = 1:numel(stat_lab)
        % Get data
        cors{sbj_ix,stat_ix}  = predCorrEpocEx.(SBJ).(stat_lab{stat_ix}).Correlations;
        boot{sbj_ix,stat_ix}  = predCorrEpocEx.(SBJ).(stat_lab{stat_ix}).Bootstrap;
        
        % Threshold data
        sig  = zeros(size(elec_sbj{sbj_ix,stat_ix}.label));
        stat = zeros(size(elec_sbj{sbj_ix,stat_ix}.label));
        fprintf('%s Significant Electrodes:\n\t',stat_lab{stat_ix});
        for e = 1:numel(elec_sbj{sbj_ix,stat_ix}.label)
            if strcmp(sig_method,'prop_boot')
                % Proportion of bootstrap values above zero
                stat(e) = 1-sum(boot{sbj_ix,stat_ix}(:,e)>0)/size(boot{sbj_ix,stat_ix},1);
                if stat(e) <= thresh
                    sig(e) = 1;
                    fprintf('%s\t',elec_sbj{sbj_ix,stat_ix}.label{e});
                end
            elseif strcmp(sig_method,'norminv')
                error('dont use this until sure its appropriate, see email Slides ready to edit on 12/11/18');
                % STDs of bootstrap above 0
                sd = std(boot{sbj_ix,stat_ix}(:,e),0,1);
                stat(e) = -norminv(thresh)*sd;
                if cors{sbj_ix,stat_ix}(e) >= stat(e)
                    sig(e) = 1;
                    fprintf('%s\t',elec_sbj{sbj_ix,stat_ix}.label{e});
                end
            end
        end
        fprintf('\nTotal n_sig = %i / %i\n',sum(sig),numel(cors{sbj_ix,stat_ix}));
        
        % Select sig elecs && elecs matching atlas
        % fn_select_elec messes up if you try to toss all elecs
        plot_elecs = intersect(good_elecs, elec_sbj{sbj_ix,stat_ix}.label(logical(sig)));
        fprintf('\nTotal plotting: %i\n',numel(plot_elecs));
        if numel(intersect(elec_sbj{sbj_ix,stat_ix}.label,plot_elecs))==0
            elec_sbj{sbj_ix,stat_ix} = {};
            good_sbj(sbj_ix,stat_ix) = false;
            warning('WARNING!!! All sig_ch are out of atlas and/or hemisphere!');
        else
            cfgs = [];
            cfgs.channel = plot_elecs;
            elec_sbj{sbj_ix,stat_ix} = fn_select_elec(cfgs, elec_sbj{sbj_ix,stat_ix});
            all_roi_labels{stat_ix} = [all_roi_labels{stat_ix}; elec_sbj{sbj_ix,stat_ix}.roi];
            all_roi_colors{stat_ix} = [all_roi_colors{stat_ix}; elec_sbj{sbj_ix,stat_ix}.roi_color];
        end
    end
    clear SBJ SBJ_vars SBJ_vars_cmd
end

%% Combine elec structs
elec = cell([numel(stat_lab) 1]);
for stat_ix = 1:numel(stat_lab)
    elec{stat_ix} = ft_appendsens([],elec_sbj{good_sbj(:,stat_ix),stat_ix});
    elec{stat_ix}.roi       = all_roi_labels{stat_ix};    % appendsens strips that field
    elec{stat_ix}.roi_color = all_roi_colors{stat_ix};    % appendsens strips that field
end

%% Load brain recon
mesh = fn_load_recon_mesh([],'mni',reg_type,hemi);

%% 3D Surface + Grids (3d, pat/mni, vol/srf, 0/1)
% f = cell(size(stat_lab));
% for stat_ix = 1:numel(stat_lab)
%     plot_name = ['GRP_' stat_lab{stat_ix} '_sig_' atlas_id '_' roi_id];
%     f{stat_ix} = figure('Name',plot_name);
plot_name = ['GRP_allModels_sig_' atlas_id '_' roi_id];
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
            ft_plot_sens(elec_tmp, 'elecshape', 'sphere', 'facecolor', elec_tmp.roi_color, 'label', lab_arg);
        end
    end
end

view(view_angle); material dull; lighting gouraud;
l = camlight;
fprintf(['To reset the position of the camera light after rotating the figure,\n' ...
    'make sure none of the figure adjustment tools (e.g., zoom, rotate) are active\n' ...
    '(i.e., uncheck them within the figure), and then hit ''l'' on the keyboard\n'])
set(f, 'windowkeypressfcn',   @cb_keyboard);
% set(f{stat_ix}, 'windowkeypressfcn',   @cb_keyboard);
% end


