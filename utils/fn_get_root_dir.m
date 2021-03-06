function [root_dir, app_dir] = fn_get_root_dir()
%% Check with user/OS/root directory, return relevant locations
if exist('/home/knight/','dir')
    root_dir = '/home/knight/';
    app_dir   = [root_dir 'hoycw/Apps/'];
elseif exist('/Volumes/hoycw_clust/','dir')
    root_dir = '/Volumes/hoycw_clust/';
    app_dir   = '/Users/colinhoy/Code/Apps/';
elseif exist('/Users/lapate/','dir')
    root_dir = '/Users/lapate/knight/';
    app_dir = '/Users/lapate/knight/hoycw/Apps/';
else
    error('root directory not found. where are you running this?');
end

end
