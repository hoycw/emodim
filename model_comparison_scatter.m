SBJs = {'CP24','IR68','IR66'};
sig_method = 'prop_boot';
thresh = 0.01;
stat_id  = {'OM','RM','RT'};
stat_lab = {'OrigMean', 'ReginaMean', 'ReginaTimeAnn'};

[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];

%% Plot
for s = 1:numel(SBJs)
    SBJ = SBJs{s};
    SBJ_vars_cmd = ['run ' root_dir 'emodim/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    sig_fname = [SBJ_vars.dirs.proc SBJ '_model_fit_corr_stats_' sig_method '_' num2str(thresh) '.mat'];
    load(sig_fname,'sig','cors');
    fprintf('%s - %i sig\n',SBJs{s},sum(any(sig,2)));

    % Get correlation bounds
    min_r = 1; max_r = -1;
    for st = 1:numel(stat_lab)
        if min(cors{st})<min_r
            min_r = min(cors{st});
        end
        if max(cors{st})>max_r
            max_r = max(cors{st});
        end
    end
        
    st_pairs = [1 2; 1 3; 2 3];
    figure('Name',SBJ,'Position', [55 285 1046 286]);
    for stp = 1:size(st_pairs,1)
        subplot(1,size(st_pairs,1),stp);
        hold on;
        % To plot bigger for sig in any of 3 models
        scatter(cors{st_pairs(stp,1)}(any(sig,2)),...
                cors{st_pairs(stp,2)}(any(sig,2)),...
                40,sig(any(sig,2),:)./1.5,'*');
        scatter(cors{st_pairs(stp,1)}(~any(sig,2)),...
                cors{st_pairs(stp,2)}(~any(sig,2)),...
                15,sig(~any(sig,2),:)./1.5,'.');
            % To plot bigger only for significant in current 2 models
%         scatter(cors{st_pairs(stp,1)}(any(sig(:,st_pairs(stp,:)),2)),...
%                 cors{st_pairs(stp,2)}(any(sig(:,st_pairs(stp,:)),2)),...
%                 40,sig(any(sig(:,st_pairs(stp,:)),2),:)./1.5,'*');
%         scatter(cors{st_pairs(stp,1)}(~any(sig(:,st_pairs(stp,:)),2)),...
%                 cors{st_pairs(stp,2)}(~any(sig(:,st_pairs(stp,:)),2)),...
%                 15,sig(~any(sig(:,st_pairs(stp,:)),2),:)./1.5,'.');
        xlim([min_r max_r]);
        ylim([min_r max_r]);
        line(xlim,ylim,'Color','k','LineStyle','--');
        xlabel(stat_lab{st_pairs(stp,1)});
        ylabel(stat_lab{st_pairs(stp,2)});
        axis square
    end
    fig_fname = [root_dir 'emodim/results/' SBJ '_model_comparison_scatter.png'];
    saveas(gcf,fig_fname);
end