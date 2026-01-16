function [ f ] = func_plot_blp_diff_revised(table_all_corr, table_all_rsn_corr, chanlocs, vtype, fband, pstate, flag_sign)

% Function to show and contrast band-limited power (BLP) indices between
% the HC and SZ groups. Upper panels show channel-wise topology, lower
% panel shows statistics on the level of resting-state networks. Showing
% channel-wise differences is optional (flag_sign). Use vtype='mean' or
% vtype='var' for mean or variance over time, respectively. Use fband from
% {'delta','theta','alpha','beta','gamma'} for the corresponding frequency
% range, respectively. Use 'mixd','frac' or 'osci' for the power spectrum
% type.
%
% Note that this function uses the topoplot() function of EEGLAB (Delorme
% & Makeig, 2004) and therefore EEGLAB (or topoplot) needs to be added to
% the path.
%
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 02/21/2025


%% collect channel-wise significant differences (exploratory)

vname_mixd = [vtype '_' fband '_mixd'];
vname_osci = [vtype '_' fband '_osci'];

table_blp_all = table_all_corr(contains(table_all_corr.meas,vname_mixd) | contains(table_all_corr.meas,vname_osci),:);
table_sign_blp = table_blp_all(table_blp_all.h_FDR==1,:);

% mixed (raw) spectra
if ~isempty(table_sign_blp)
    ind_mixd = (contains(table_sign_blp.meas,vname_mixd)==1); % select mixd BLP
    p_mixd = table_sign_blp(ind_mixd,:);
else
    p_mixd = [];
end

% oscillatory spectra
if ~isempty(table_sign_blp)
    ind_osci = (contains(table_sign_blp.meas,vname_osci)==1); % select oscillatory BLP
    p_osci = table_sign_blp(ind_osci,:);
else
    p_osci = [];
end

%% topoplots
f = figure('color','w','units','normalized','outerposition',[0 0.05 1 0.9]);
irasa_vars = {'mixd','osci'};

%% obtain channel-wise results

for iv = 1:length(irasa_vars)
    switch iv
        case 1
            stype = 'mixd';
            stype_title = 'Mixed';
            p_blp = p_mixd;
            vname = vname_mixd;
        case 2
            stype = 'osci';
            stype_title = 'Oscillatory';
            p_blp = p_osci;
            vname = vname_osci;
    end

    table_blp = table_blp_all(contains(table_blp_all.meas,vname),:);
    table_blp = sortrows(table_blp,'chID','ascend');

    if strcmp(vtype,'mean')
        cblabel = ['\mu(\' fband ' -BLP^{' pstate '}_{\it{' stype '}})']; % label for colorbar
    elseif strcmp(vtype,'var')
        cblabel = ['\sigma(\' fband ' -BLP^{' pstate '}_{\it{' stype '}})']; % label for colorbar
    else
        error('Set vtype to mean or var')
    end

    % get group expected values
    v_hc = table_blp.E_hc;
    v_sz = table_blp.E_sz;
    v_full = [v_hc; v_sz];

    % denote significant channel-wise differences
    if ~isempty(p_blp)
        ch_sign = p_blp.chID;
    else
        ch_sign = [];
    end

    %% plotting
    switch iv
        case 1
            ax1 = subplot(2,4,2*(iv-1)+5);

            % HC topoplot
            if flag_sign % denote channel-wise between-group differences with asterisk
                topoplot(v_hc, chanlocs, 'emarker', {'.','k',12,2}, 'emarker2', {ch_sign,'*','k',12,2.5}, 'whitebk', 'on');
            else
                topoplot(v_hc, chanlocs, 'emarker', {'.','k',12,2}, 'whitebk', 'on');
            end

            clim([min(v_full) max(v_full)])
            cb1 = colorbar;
            % cb1.Label.String = cblabel;
            % cb1.Label.Rotation = 270;
            cb1.FontSize = 14;
            % cb1.Label.FontSize = 16;
            % cb1.Label.VerticalAlignment = 'bottom';
            if iv == 1
                % ylabel([stype_title ' power, \mu(\' fband '-BLP^{' pstate '}_{\it{' stype '}})' ],'FontSize',16,'visible','on')
            end
            title({'Healthy',[stype_title ' Power']},'FontSize',20,'visible','on')

            ax2 = subplot(2,4,2*(iv-1)+6);

            % SZ topoplot
            if flag_sign % denote channel-wise between-group differences with asterisk
                topoplot(v_sz, chanlocs, 'emarker', {'.','k',12,2}, 'emarker2', {ch_sign,'*','k',12,2.5}, 'whitebk', 'on');
            else
                topoplot(v_sz, chanlocs, 'emarker', {'.','k',12,2}, 'whitebk', 'on');
            end

            clim([min(v_full) max(v_full)])
            cb1 = colorbar;
            % cb1.Label.String = cblabel;
            % cb1.Label.Rotation = 270;
            cb1.FontSize = 14;
            % cb1.Label.FontSize = 16;
            % cb1.Label.VerticalAlignment = 'bottom';

            title({'Schizophrenia',[stype_title ' Power']},'FontSize',20,'visible','on')

            % cchart = summer;
            % c_hc = cchart(1,:);
            % c_sz = cchart(end,:);
            c_hc = [0, 0.45, 1];
            c_sz = [1, 0.45, 0];
        case 2
            ax3 = subplot(2,4,2*(iv-1)+5);

            % HC topoplot
            if flag_sign % denote channel-wise between-group differences with asterisk
                topoplot(v_hc, chanlocs, 'emarker', {'.','k',12,2}, 'emarker2', {ch_sign,'*','k',12,2.5}, 'whitebk', 'on');
            else
                topoplot(v_hc, chanlocs, 'emarker', {'.','k',12,2}, 'whitebk', 'on');
            end

            clim([min(v_full) max(v_full)])
            cb1 = colorbar;
            % cb1.Label.String = cblabel;
            % cb1.Label.Rotation = 270;
            cb1.FontSize = 14;
            % cb1.Label.FontSize = 16;
            % cb1.Label.VerticalAlignment = 'bottom';
            if iv == 2
                % ylabel([stype_title ' power, \mu(\' fband '-BLP^{' pstate '}_{\it{' stype '}})' ],'FontSize',16,'visible','on')
            end
            title({'Healthy',[stype_title ' Power']},'FontSize',20,'visible','on')

            ax4 = subplot(2,4,2*(iv-1)+6);

            % SZ topoplot
            if flag_sign % denote channel-wise between-group differences with asterisk
                topoplot(v_sz, chanlocs, 'emarker', {'.','k',12,2}, 'emarker2', {ch_sign,'*','k',12,2.5}, 'whitebk', 'on');
            else
                topoplot(v_sz, chanlocs, 'emarker', {'.','k',12,2}, 'whitebk', 'on');
            end

            clim([min(v_full) max(v_full)])
            cb1 = colorbar;
            % cb1.Label.String = cblabel;
            % cb1.Label.Rotation = 270;
            cb1.FontSize = 14;
            % cb1.Label.FontSize = 16;
            % cb1.Label.VerticalAlignment = 'bottom';

            title({'Schizophrenia',[stype_title ' Power']},'FontSize',20,'visible','on')

            % cchart = winter;
            % c_hc = cchart(1,:);
            % c_sz = cchart(end,:);
            c_hc = [0, 0.45, 1];
            c_sz = [1, 0.45, 0];
    end

    %% RSN-wise results

    rsn_blp = table_all_rsn_corr(contains(table_all_rsn_corr.meas,vname),:);
    rsn_blp = sortrows(rsn_blp,'rsnID','ascend');

    if strcmp(vtype,'mean')
        str_title = [stype_title ' RSN-wise \mu(\' fband '-BLP^{' pstate '}_{\it{' stype '}})' ];% indicating frequency range in title
        str_ylabel = ['\mu(\' fband ' -BLP^{' pstate '}_{\it{' stype '}})']; % label for colorbar
    elseif strcmp(vtype,'var')
        str_title = [stype_title ' RSN-wise \sigma(\' fband '-BLP^{' pstate '}_{\it{' stype '}})' ];% indicating frequency range in title
        str_ylabel = ['\sigma(\' fband ' -BLP^{' pstate '}_{\it{' stype '}})']; % label for colorbar
    else
        error('Set vtype to mean or var')
    end


    rsn_list = rsn_blp.rsn;
    nrsn = length(rsn_list);

    % create box plots
    subplot(2,4,[2*(iv-1)+1 2*(iv-1)+2])
    for rsn = 1:nrsn
        v_hc = rsn_blp.v_hc{rsn};
        N_hc = length(v_hc);
        v_sz = rsn_blp.v_sz{rsn};
        N_sz = length(v_sz);

        % plot values
        plot((rsn-0.2)*ones(N_hc,1),v_hc,'o','color',0.75*c_hc,'MarkerSize',8,'LineWidth',1)
        hold on
        plot((rsn+0.2)*ones(N_sz,1),v_sz,'o','color',0.75*c_sz,'MarkerSize',8,'LineWidth',1)

        % draw boxes
        pm_hc = fill([rsn-0.05 rsn-0.35 rsn-0.35 rsn-0.05],[quantile(v_hc,0.25), quantile(v_hc,0.25), quantile(v_hc,0.75) quantile(v_hc,0.75)],c_hc,'LineWidth',1.5,'FaceAlpha',0.4);
        pm_sz = fill([rsn+0.05 rsn+0.35 rsn+0.35 rsn+0.05],[quantile(v_sz,0.25), quantile(v_sz,0.25), quantile(v_sz,0.75) quantile(v_sz,0.75)],c_sz,'LineWidth',1.5,'FaceAlpha',0.4);

        % draw percentiles (5, 25, 75, 95)
        plot([rsn-0.2 rsn-0.2],[quantile(v_hc,0.75) quantile(v_hc,0.95)],'k-','LineWidth',1.5)
        plot([rsn-0.3 rsn-0.1],[quantile(v_hc,0.95) quantile(v_hc,0.95)],'k-','LineWidth',1.5)
        plot([rsn-0.2 rsn-0.2],[quantile(v_hc,0.25) quantile(v_hc,0.05)],'k-','LineWidth',1.5)
        plot([rsn-0.3 rsn-0.1],[quantile(v_hc,0.05) quantile(v_hc,0.05)],'k-','LineWidth',1.5)

        plot([rsn+0.2 rsn+0.2],[quantile(v_sz,0.75) quantile(v_sz,0.95)],'k-','LineWidth',1.5)
        plot([rsn+0.1 rsn+0.3],[quantile(v_sz,0.95) quantile(v_sz,0.95)],'k-','LineWidth',1.5)
        plot([rsn+0.2 rsn+0.2],[quantile(v_sz,0.25) quantile(v_sz,0.05)],'k-','LineWidth',1.5)
        plot([rsn+0.1 rsn+0.3],[quantile(v_sz,0.05) quantile(v_sz,0.05)],'k-','LineWidth',1.5)

        % plot means and medians
        plot([rsn-0.35 rsn-0.05],[mean(v_hc) mean(v_hc)],'-','color',0.5*c_hc, 'LineWidth',3)
        plot([rsn-0.35 rsn-0.05],[median(v_hc) median(v_hc)],':','color',0.5*c_hc, 'LineWidth',3)
        plot([rsn+0.05 rsn+0.35],[mean(v_sz) mean(v_sz)],'-','color',0.5*c_sz, 'LineWidth',3)
        plot([rsn+0.05 rsn+0.35],[median(v_sz) median(v_sz)],':','color',0.5*c_sz, 'LineWidth',3)
    end

    ymin = min(get(gca,'ylim'));
    ymax = max(get(gca,'ylim'));
    yrange = ymax - ymin;
    ylim([ymin, ymin + yrange*1.2])

    % plot significant differences
    for i = 1:nrsn
        p = rsn_blp.p(i);
        h_FDR = rsn_blp.h_FDR(i);
        if h_FDR == 1
            plot([i-0.3 i+0.3],[ymin+1.05*yrange ymin+1.05*yrange],'k-','LineWidth',2.5)

            if (p < 0.05) && (p > 0.001)
                plot(i,ymin+1.1*yrange,'k*','MarkerSize',10,'LineWidth',1.5)
            elseif (p < 0.001) && (p > 0.0001)
                plot([i-0.1 i+0.1],[ymin+1.1*yrange ymin+1.1*yrange],'k*','MarkerSize',10,'LineWidth',1.5)
            else
                plot([i-0.2 i i+0.2],[ymin+1.1*yrange ymin+1.1*yrange ymin+1.1*yrange],'k*','MarkerSize',10,'LineWidth',1.5)
            end
        elseif p < 0.05
            plot([i-0.3 i+0.3],[ymin+1.05*yrange ymin+1.05*yrange],':','color',[0.3,0.3,0.3],'LineWidth',1.5)
        end
    end

    % separate Global results from RSN results
    plot([nrsn-0.5 nrsn-0.5],[ymin, ymin + yrange*1.2],'k--','LineWidth',1.5)
    plot([nrsn+0.5 nrsn+0.5],[ymin, ymin + yrange*1.2],'k-','LineWidth',1)

    % set panel parameters
    xlim([0.5 nrsn+2])
    set(gca,'LineWidth',1,'box','on','XTick',(1:nrsn),'XTickLabel',rsn_list,'FontSize',14);
    legend([pm_hc,pm_sz],{'HC','SZ'},'FontSize',14,'location','east');
    xlabel('Resting-state networks','FontSize',16)
    ylabel(str_ylabel, 'FontSize', 16, 'visible', 'on')
    title(str_title,'FontSize',20,'visible','on')

end

% set colormaps
for iv = 1:length(irasa_vars)
    switch iv
        case 1
            colormap(ax1,jet)
            colormap(ax2,jet)
        case 2
            colormap(ax3,jet)
            colormap(ax4,jet)
    end
end


end