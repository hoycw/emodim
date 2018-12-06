function fn_compile_elec_struct(SBJ,pipeline_id,view_space,reg_type)
%% Compile ROI info (biggest changes are from single electrodes into bipolar pairs)
%   Goes from smallest to largest (ELEC1-ELEC2, ELEC2-ELEC3, etc.)
%   Pairs are drawn from imported data labels, then preprocessing logic is applied
%
% INPUTS:
%   SBJ [str] - name of subject
%   pipeline_id [str] - name of analysis pipeline
%   view_space [str] - {'pat','mni'} select patient native or mni group space
%   reg_type [str] - {'v', 's'} choose volume-based or surface-based registration

% Set up paths
[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
addpath([root_dir 'emodim/scripts/']);
addpath([root_dir 'emodim/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Load variables
eval(['run ' root_dir 'emodim/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run ' root_dir 'emodim/scripts/proc_vars/' pipeline_id '_proc_vars.m']);

%% Load Elec struct
if strcmp(view_space,'pat')
    if ~isempty(reg_type)
        warning(['view_space = "' view_space '" does not require reg_type; ignoring reg_type = ' reg_type]);
    end
    elec_ext = view_space;
    elec_name = SBJ_vars.recon.elec_pat;
elseif strcmp(view_space,'mni')
    elec_ext = [view_space '_' reg_type];
    if strcmp(reg_type,'v') % volume based
        elec_name = SBJ_vars.recon.elec_mni_v;
    elseif strcmp(reg_type,'s') %surface based
        elec_name = SBJ_vars.recon.elec_mni_s;
        if isempty(elec_name)
            warning(['WARNING!!! elec_mni_s does not exist for ' SBJ ', nothing being compiled.']);
            return
        end
    else
        error(['Unknown reg_type: ' reg_type]);
    end
else
    error(['Unknown view_space: ' view_space]);
end
slash = strfind(elec_name,'/'); elec_suffix = elec_name(slash(end)+numel(SBJ)+2:end-4);
tmp = load(elec_name);
elec_var_name = fieldnames(tmp);
if ~strcmp(elec_var_name,elec_suffix)
    warning(['\t!!!! ' SBJ ' elec names in variable and file names do not match! file=' elec_suffix '; var=' elec_var_name{1}]);
end
eval(['elec = tmp.' elec_var_name{1} ';']); clear tmp;

%% Channel Label Corrections
% Strip Pre/Suffix if Necessary
for ch_ix = 1:numel(elec.label)
    if isfield(SBJ_vars.ch_lab,'prefix')
        elec.label{ch_ix} = strrep(elec.label{ch_ix},SBJ_vars.ch_lab.prefix,'');
    end
    if isfield(SBJ_vars.ch_lab,'suffix')
        elec.label{ch_ix} = strrep(elec.label{ch_ix},SBJ_vars.ch_lab.suffix,'');
    end
end
% Fix any mislabeled channels
if isfield(SBJ_vars.ch_lab,'mislabel')
    for ch_ix = 1:numel(SBJ_vars.ch_lab.mislabel)
        % Future edit: search for the bad label across data, eeg, evnt
        elec.label(strcmp(elec.label,SBJ_vars.ch_lab.mislabel{ch_ix}(1))) = SBJ_vars.ch_lab.mislabel{ch_ix}(2);
    end
end

%% Load data
if numel(SBJ_vars.raw_file)>1
    block_suffix = strcat('_',SBJ_vars.block_name{1});
else
    block_suffix = SBJ_vars.block_name{1};   % should just be ''
end
import_filename = [SBJ_vars.dirs.import SBJ '_',num2str(proc_vars.resample_freq),'hz',block_suffix,'.mat'];
load(import_filename);

% % Original (single electrode) labels
% import  = load([SBJ_vars.dirs.import SBJ '_1000hz.mat']);
% raw_lab = import.data.label;

%% Select imported channels
cfg = []; cfg.channel = data.label;
elec = fn_select_elec(cfg,elec);

% Order them to match data.label
elec = fn_reorder_elec(elec, data.label);

%% Apply montage per probe
left_out_ch = {};
elec_labels = {};
elec_types  = {};
danger_name = false([1 numel(SBJ_vars.ch_lab.probes)]);
name_holder = cell([2 numel(SBJ_vars.ch_lab.probes)]);
elec_reref  = cell([1 numel(SBJ_vars.ch_lab.probes)]);
for d = 1:numel(SBJ_vars.ch_lab.probes)
    cfg = [];
    cfg.channel = ft_channelselection(strcat(SBJ_vars.ch_lab.probes{d},'*'), data.label);
    probe_elec  = fn_select_elec(cfg,elec);
    %     probe_data = ft_selectdata(cfg,data);   % Grab data from this probe to plot in PSD comparison
    %     probe_data.elec = fn_elec_ch_select(elec,cfg.channel);
    
    % Check if the names of these elecs will cause problems
    eeg1010_match = strfind(probe_elec.label,'AF');
    if ~isempty([eeg1010_match{:}])
        danger_name(d)   = true;
        name_holder{1,d} = probe_elec.label;
        name_holder{2,d} = fn_generate_random_strings(numel(probe_elec.label),'',10);
        probe_elec.label = name_holder{2,d};
    end
    
    % Create referencing scheme
    if strcmp(SBJ_vars.ch_lab.ref_type{d},'BP')
        cfg.montage.labelold = cfg.channel;
        [cfg.montage.labelnew, cfg.montage.tra, left_out_ch{d}] = fn_create_ref_scheme_bipolar(cfg.channel);
        cfg.updatesens = 'yes';
        elec_reref{d} = ft_apply_montage(probe_elec, cfg.montage);%, 'feedback', 'none', 'keepunused', 'no', 'balancename', bname);
        %     data_reref{d} = ft_preprocessing(cfg, probe_data);
    else
        elec_reref{d} = probe_elec;
    end
    if d==1
        elec_labels = elec_reref{d}.label;
        elec_types  = repmat(SBJ_vars.ch_lab.probe_type(d),size(elec_reref{d}.label));
    else
        elec_labels = cat(find(size(elec_labels)>1), elec_labels, elec_reref{d}.label);
        elec_types  = cat(find(size(elec_types)>1), elec_types, repmat(SBJ_vars.ch_lab.probe_type(d),size(elec_reref{d}.label)));
    end
end

%% Recombine
cfg = [];
elec = ft_appendsens(cfg,elec_reref{:});
elec.type = 'ieeg';
for e = 1:numel(elec.chantype)
    elec.chantype{e} = elec_types{strcmp(elec.label{e},elec_labels)};
end

% Re-label any problematic channel labels
if any(danger_name)
    for d_ix = find(danger_name)
        for s_ix = 1:numel(name_holder{2,d_ix})
            elec.label{strcmp(elec.label,name_holder{2,2}{s_ix})} = name_holder{1,d_ix}{s_ix};
        end
    end
end

%% Save data
% Check if elec.cfg.previosu got ridiculously large, and keep only first
var_stats = whos('elec');
if var_stats.bytes>1000000
    elec.cfg = rmfield(elec.cfg,'previous');
end
output_filename = strcat(SBJ_vars.dirs.recon,SBJ,'_elec_',pipeline_id,'_',elec_ext,'.mat');
fprintf('============== Saving %s ==============\n',output_filename);
save(output_filename, '-v7.3', 'elec');

end
