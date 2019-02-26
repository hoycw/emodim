%% Run elec sorting
SBJs = {'CP24','IR61','IR66','IR68','IR71','IR74','IR75','IR78','IR84'};
  
%% Directories
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'hoycw/Apps/fieldtrip/'];
elseif exist('/Users/lapate/','dir');root_dir = '/Users/lapate/knight/';ft_dir = '/Users/lapate/knight/hoycw/Apps/fieldtrip/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%%
addpath([root_dir 'emodim/scripts/']);
addpath([root_dir 'emodim/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Compare to atlases
for s = 1:numel(SBJs)
    SBJ = SBJs{s};
    SBJ_vars_cmd = ['run ' root_dir 'emodim/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    % basic pat
    load(strcat(SBJ_vars.dirs.recon,SBJ,'_elec_main_ft_pat.mat'));
    elec = fn_reorder_elec(elec,{});
    save(strcat(SBJ_vars.dirs.recon,SBJ,'_elec_main_ft_pat.mat'),'-v7.3','elec');
    % basic mni
    load(strcat(SBJ_vars.dirs.recon,SBJ,'_elec_main_ft_mni_v.mat'));
    elec = fn_reorder_elec(elec,{});
    save(strcat(SBJ_vars.dirs.recon,SBJ,'_elec_main_ft_mni_v.mat'),'-v7.3','elec');
    % pat DK
    load([SBJ_vars.dirs.recon,SBJ,'_elec_main_ft_pat_DK.mat']);
    elec = fn_reorder_elec(elec,{});
    save([SBJ_vars.dirs.recon,SBJ,'_elec_main_ft_pat_DK.mat'],'-v7.3','elec');
    % pt Dx
    load([SBJ_vars.dirs.recon,SBJ,'_elec_main_ft_pat_Dx.mat']);
    elec = fn_reorder_elec(elec,{});
    save([SBJ_vars.dirs.recon,SBJ,'_elec_main_ft_pat_Dx.mat'],'-v7.3','elec');
    
    clear SBJ SBJ_vars SBJ_vars_cmd
end
