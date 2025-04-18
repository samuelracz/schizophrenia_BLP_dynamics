function [ f ] = func_plot_EC_EO_BLP_response(struct_hc_stat, struct_sz_stat, chanlocs, vtype, fband, flag_sign)


%% topoplots
f = figure('color','w','units','normalized','outerposition',[0.025 0.05 0.95 0.8]);
irasa_vars = {'mixd','frac','osci'};
pstate = 'dS';

%% obtain channel-wise results

for iv = 1:length(irasa_vars)
    % HC
    switch iv
        case 1
            stype = 'mixd';
            stype_title = 'Mixed';
        case 2
            stype = 'frac';
            stype_title = 'Fractal';
        case 3
            stype = 'osci';
            stype_title = 'Oscillatory';
    end
    
    table_blp_hc = sortrows(struct_hc_stat.(irasa_vars{iv}).(fband).(vtype),'chID','ascend');

    if strcmp(vtype,'mean')
        str_title = [stype_title ' power, \Delta\mu(\' fband '-BLP^{' pstate '}_{\it{' stype '}})' ];% indicating frequency range in title
        cblabel = ['\Delta\mu(\' fband ' -BLP^{' pstate '}_{\it{' stype '}})']; % label for colorbar
    elseif strcmp(vtype,'var')
        str_title = [stype_title ' power, \Delta\sigma(\' fband '-BLP^{' pstate '}_{' stype '})' ];% indicating frequency range in title
        cblabel = ['\Delta\sigma(\' fband ' -BLP^{' pstate '}_{\it{' stype '}})']; % label for colorbar
    else
        error('Set vtype to mean or var')
    end

    % get group expected values
    v_hc_ec = table_blp_hc.E_ec;
    v_hc_eo = table_blp_hc.E_eo;
    v_diff = v_hc_ec - v_hc_eo;

    % denote significant channel-wise differences
    p_blp_hc = table_blp_hc(table_blp_hc.h_FDR==1,:);
    if ~isempty(p_blp_hc)
        ch_sign = p_blp_hc.chID;
    else
        ch_sign = [];
    end

    %% plotting
    switch iv
        case 1
            ax1 = subplot(2,3,iv);
        case 2
            ax2 = subplot(2,3,iv);
        case 3
            ax3 = subplot(2,3,iv);
    end

    % HC topoplot
    if flag_sign % denote channel-wise between-group differences with asterisk
        topoplot(v_diff, chanlocs, 'emarker', {'.','k',20,2}, 'emarker2', {ch_sign,'*','k',18,2.5}, 'whitebk', 'on');
    else
        topoplot(v_diff, chanlocs, 'emarker', {'.','k',20,2}, 'whitebk', 'on');
    end

    clim([min(v_diff) max(v_diff)])
    cb = colorbar;
    cb.Label.String = cblabel;
    cb.Label.Rotation = 270;
    cb.FontSize = 16;
    cb.Label.FontSize = 28;
    cb.Label.VerticalAlignment = 'bottom';
    if iv == 1
        ylabel('Healthy EC-EO','FontSize',24,'visible','on')
    end
    title(str_title,'FontSize',28,'visible','on')

    %% SZ
    switch iv
        case 1
            stype = 'mixd';
            stype_title = 'Mixed';
        case 2
            stype = 'frac';
            stype_title = 'Fractal';
        case 3
            stype = 'osci';
            stype_title = 'Oscillatory';
    end
    
    table_blp_sz = sortrows(struct_sz_stat.(irasa_vars{iv}).(fband).(vtype),'chID','ascend');

    if strcmp(vtype,'mean')
        str_title = [stype_title ' power, \Delta\mu(\' fband '-BLP^{' pstate '}_{\it{' stype '}})' ];% indicating frequency range in title
        cblabel = ['\Delta\mu(\' fband ' -BLP^{' pstate '}_{\it{' stype '}})']; % label for colorbar
    elseif strcmp(vtype,'var')
        str_title = [stype_title ' power, \Delta\sigma(\' fband '-BLP^{' pstate '}_{' stype '})' ];% indicating frequency range in title
        cblabel = ['\Delta\sigma(\' fband ' -BLP^{' pstate '}_{\it{' stype '}})']; % label for colorbar
    else
        error('Set vtype to mean or var')
    end

    % get group expected values
    v_sz_ec = table_blp_sz.E_ec;
    v_sz_eo = table_blp_sz.E_eo;
    v_diff = v_sz_ec - v_sz_eo;

    % denote significant channel-wise differences
    p_blp_sz = table_blp_sz(table_blp_sz.h_FDR==1,:);
    if ~isempty(p_blp_sz)
        ch_sign = p_blp_sz.chID;
    else
        ch_sign = [];
    end

    %% plotting
    switch iv
        case 1
            ax4 = subplot(2,3,iv+3);
        case 2
            ax5 = subplot(2,3,iv+3);
        case 3
            ax6 = subplot(2,3,iv+3);
    end

    % HC topoplot
    if flag_sign % denote channel-wise between-group differences with asterisk
        topoplot(v_diff, chanlocs, 'emarker', {'.','k',20,2}, 'emarker2', {ch_sign,'*','k',18,2.5}, 'whitebk', 'on');
    else
        topoplot(v_diff, chanlocs, 'emarker', {'.','k',20,2}, 'whitebk', 'on');
    end

    clim([min(v_diff) max(v_diff)])
    cb = colorbar;
    cb.Label.String = cblabel;
    cb.Label.Rotation = 270;
    cb.FontSize = 16;
    cb.Label.FontSize = 28;
    cb.Label.VerticalAlignment = 'bottom';
    if iv == 1
        ylabel('Schizophrenia EC-EO','FontSize',24,'visible','on')
    end
    title(str_title,'FontSize',28,'visible','on')
end

% set colormaps
for iv = 1:length(irasa_vars)
    switch iv
        case 1
            colormap(ax1,summer)
            colormap(ax4,summer)
        case 2
            colormap(ax2,autumn)
            colormap(ax5,autumn)
        case 3
            colormap(ax3,winter)
            colormap(ax6,winter)
    end
end

end