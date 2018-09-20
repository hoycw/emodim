function [root_dir, ft_dir] = fn_get_root_dir()
%% Check with user/OS/root directory, return relevant locations
if exist('/home/knight/','dir')
    root_dir = '/home/knight/';
    ft_dir   = ['/home/knight/hoycw/Apps/fieldtrip/'];
elseif exist('/Volumes/hoycw_clust/','dir')
    root_dir = '/Volumes/hoycw_clust/';
    ft_dir   = '/Users/colinhoy/Code/Apps/fieldtrip/';
else
    error('root directory not found. where are you running this?');
end

end
