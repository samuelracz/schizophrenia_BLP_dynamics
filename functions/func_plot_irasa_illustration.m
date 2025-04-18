function [fig1, output_hc, output_sz] = func_plot_irasa_illustration(struct_hc, struct_sz, ch_lab, flag_fig)

% Function to illustrate the IRASA method on grand average spectra in the
% HC and SZ groups. Channel selection is via channel labels.
%
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 04/18/2025

%% variables and constants

spec_hc = struct_hc.spec; %............HC spectra
spec_sz = struct_sz.spec; %............SZ spectra

chlab = struct_hc.chlab; %.............channel labels
if strcmp(ch_lab,'global')
    ch = -1;
else
    ch = find(strcmp(chlab,ch_lab)); %....channel index
end

freq = spec_hc(1).plaw.freq; %.........frequency vector
N = length(freq); %....................# of frequencies

n_hc = length(spec_hc); %..............# of healthy controls
n_sz = length(spec_sz); %..............# of szhizophrenia patients

%% obtain HC data

% preallocations
hc_mixd = zeros(N,n_hc);
hc_frac = zeros(N,n_hc);
hc_osci = zeros(N,n_hc);

for n = 1:n_hc % iterate over HC participants
    tmp = spec_hc(n).plaw;
    if ch == -1
        hc_mixd(:,n) = mean(tmp.mixd,2);
        hc_frac(:,n) = mean(tmp.frac,2);
        hc_osci(:,n) = mean(log(tmp.mixd),2) - mean(log(tmp.frac),2);
    else
        hc_mixd(:,n) = tmp.mixd(:,ch1);
        hc_frac(:,n) = tmp.frac(:,ch1);
        hc_osci(:,n) = (log(tmp.mixd(:,ch1)) - log(tmp.frac(:,ch1)));
    end
end

%% obtain SZ data

% preallocations
sz_mixd = zeros(N,n_sz);
sz_frac = zeros(N,n_sz);
sz_osci = zeros(N,n_sz);

for n = 1:n_sz % iterate over SZ participants
    tmp = spec_sz(n).plaw;
    if ch == -1
        sz_mixd(:,n) = mean(tmp.mixd,2);
        sz_frac(:,n) = mean(tmp.frac,2);
        sz_osci(:,n) = mean(log(tmp.mixd),2) - mean(log(tmp.frac),2);
    else
        sz_mixd(:,n) = tmp.mixd(:,ch);
        sz_frac(:,n) = tmp.frac(:,ch);
        sz_osci(:,n) = (log(tmp.mixd(:,ch)) - log(tmp.frac(:,ch)));
    end
end

%% store output
output_hc = struct(...
    'freq', freq,...
    'mixd', hc_mixd,...
    'frac', hc_frac,...
    'osci', hc_osci);

output_sz = struct(...
    'freq', freq,...
    'mixd', sz_mixd,...
    'frac', sz_frac,...
    'osci', sz_osci);



%% plotting

xlims = [1 45];

if flag_fig
    fig1 = figure('color','w','units','normalized','outerposition',[0.35 0.05 0.3 0.9]);

    %% ploting spectra from HC

    v_mixd = hc_mixd;
    v_frac = hc_frac;
    v_osci = hc_osci;
    xlab = 'Log(frequency)';
    ylab = 'Log(amplitude)';

    m_mixd = mean(v_mixd,2); % mean raw (mixd) power spectrum
    m_mixd_p = m_mixd + std(v_mixd,[],2)./sqrt(n_hc); % mean + SEM
    m_mixd_n = m_mixd - std(v_mixd,[],2)./sqrt(n_hc); % mean - SEM

    m_frac = mean(v_frac,2); % mean fractal power spectrum
    % SEM range is not plotted for fractal power spectrum for clarity
   
    % main axis for mixd and frac spectra
    ax1_hc = subplot(2,1,1,'Position',[0.14 0.58, 0.76, 0.36]);

    % error ranges
    fill([freq; flipud(freq)],[m_mixd_n; flipud(m_mixd_p)],[0 0.45 1],'FaceAlpha',0.3,'EdgeAlpha',0.1)
    xlim(xlims)
    set(ax1_hc,'FontSize',18,'LineWidth',1.5,'XScale','log','YScale','log')
    hold on

    % plotting spectra
    pm = loglog(freq, m_mixd, 'color', [0, 0.45, 1], 'LineWidth', 3);
    pf = loglog(freq, m_frac, 'color', [0, 0.3, 0.5], 'LineWidth', 3);

    % frequency band boundaries
    ymin_tmp = min(get(gca,'ylim'));
    ymax_tmp = max(get(gca,'ylim'));
    plot([4 4],[ymin_tmp ymax_tmp],'k--','LineWidth',1.5)
    plot([8 8],[ymin_tmp ymax_tmp],'k--','LineWidth',1.5)
    plot([13 13],[ymin_tmp ymax_tmp],'k--','LineWidth',1.5)
    plot([25 25],[ymin_tmp ymax_tmp],'k--','LineWidth',1.5)
    ylim([ymin_tmp ymax_tmp])

    % legend
    legend([pm,pf],{'mixed','fractal'},'location','northeast','FontSize',20)

    % plot labels
    title(['Healthy - ' ch_lab],'FontSize',28)
    xlabel(xlab,'FontSize',24)
    ylabel(ylab,'FontSize',24)

    % inset plot for oscillatory spectrum
    ax2_hc = axes('Position',[0.19, 0.61, 0.2, 0.12]);

    m_osci = mean(v_osci,2); % mean oscillatory power
    m_osci_p = mean(v_osci,2) + std(v_osci,[],2)./sqrt(n_hc); % mean + SEM
    m_osci_n = mean(v_osci,2) - std(v_osci,[],2)./sqrt(n_hc); % mean - SEM

    % error ranges
    fill([freq; flipud(freq)],[m_osci_n; flipud(m_osci_p)],[0.6 0.6 0.75],'FaceAlpha',0.3,'EdgeAlpha',0.1)
    hold on

    % plot mean
    plot(freq, m_osci, 'color', [0 0 0.25], 'LineWidth', 1);
    xlim(xlims)
    set(ax2_hc,'XTick',[1 10 20 30 40],'XTickLabel',{'1','10','20','30','40'},'FontSize',12)

    % annotation for inset plot
    a_lo = annotation('textbox');
    a_lo.String = 'oscillatory';
    a_lo.FontSize = 20;
    a_lo.Color = 'k';
    a_lo.Position = [0.2, 0.66, 0.3, 0.1];
    a_lo.LineStyle = 'none';

    % get y limits
    ylim_mixd_hc = get(ax1_hc,'ylim');
    ylim_osci_hc = get(ax2_hc,'ylim');

    %% plotting spectra from SZ

    v_mixd = sz_mixd;
    v_frac = sz_frac;
    v_osci = sz_osci;
    xlab = 'Log(frequency)';
    ylab = 'Log(amplitude)';

    m_mixd = mean(v_mixd,2); % mean raw (mixd) power spectrum
    m_mixd_p = m_mixd + std(v_mixd,[],2)./sqrt(n_hc); % mean + SEM
    m_mixd_n = m_mixd - std(v_mixd,[],2)./sqrt(n_hc); % mean - SEM

    m_frac = mean(v_frac,2); % mean fractal power spectrum
    % SEM range is not plotted for fractal power spectrum for clarity

    % main axis for mixd and frac spectra
    ax1_sz = subplot(2,1,2,'Position',[0.14 0.08, 0.76, 0.36]);

    % error ranges
    fill([freq; flipud(freq)],[m_mixd_n; flipud(m_mixd_p)],[1 0.45 0],'FaceAlpha',0.3,'EdgeAlpha',0.1)
    xlim(xlims)
    set(ax1_sz,'FontSize',18,'LineWidth',1.5,'XScale','log','YScale','log')
    hold on

    % plotting spectra
    pm = loglog(freq, m_mixd, 'color', [1 0.45 0], 'LineWidth', 3);
    pf = loglog(freq, m_frac, 'color', [0.5 0.3 0], 'LineWidth', 3);

    % frequency band boundaries
    ymin_tmp = min(get(gca,'ylim'));
    ymax_tmp = max(get(gca,'ylim'));
    plot([4 4],[ymin_tmp ymax_tmp],'k--','LineWidth',1.5)
    plot([8 8],[ymax_tmp ymin_tmp],'k--','LineWidth',1.5)
    plot([13 13],[ymin_tmp ymax_tmp],'k--','LineWidth',1.5)
    plot([25 25],[ymax_tmp ymin_tmp],'k--','LineWidth',1.5)
    ylim([ymin_tmp ymax_tmp])

    % legend
    legend([pm,pf],{'mixed','fractal'},'location','northeast','FontSize',20)

    % plot labels
    title(['Schizophrenia - ' ch_lab],'FontSize',28)
    xlabel(xlab,'FontSize',24)
    ylabel(ylab,'FontSize',24)

    % inset plot for oscillatory spectrum
    ax2_sz = axes('Position',[0.19, 0.11, 0.2, 0.12]);

    m_osci = mean(v_osci,2); % mean oscillatory power
    m_osci_p = mean(v_osci,2) + std(v_osci,[],2)./sqrt(n_hc); % mean + SEM
    m_osci_n = mean(v_osci,2) - std(v_osci,[],2)./sqrt(n_hc); % mean - SEM

    % error ranges
    fill([freq; flipud(freq)],[m_osci_n; flipud(m_osci_p)],[0.75 0.6 0.6],'FaceAlpha',0.3,'EdgeAlpha',0.1)
    hold on

    % plot mean
    plot(freq, m_osci, 'color', [0.25 0 0], 'LineWidth', 1);
    xlim(xlims)
    set(ax2_sz,'XTick',[1 10 20 30 40],'XTickLabel',{'1','10','20','30','40'},'FontSize',12)

    % annotation for inset plot
    a_lo = annotation('textbox');
    a_lo.String = 'oscillatory';
    a_lo.FontSize = 20;
    a_lo.Color = 'k';
    a_lo.Position = [0.2, 0.16, 0.3, 0.1];
    a_lo.LineStyle = 'none';

    % get y limits
    ylim_mixd_sz = get(ax1_sz,'ylim');
    ylim_osci_sz = get(ax2_sz,'ylim');


    %% set y-axis limits
    set(ax1_hc,'ylim', [min([ylim_mixd_hc(:); ylim_mixd_sz(:)]), max([ylim_mixd_hc(:); ylim_mixd_sz(:)])])
    set(ax1_sz,'ylim', [min([ylim_mixd_hc(:); ylim_mixd_sz(:)]), max([ylim_mixd_hc(:); ylim_mixd_sz(:)])])
    set(ax2_hc,'ylim', [min([ylim_osci_hc(:); ylim_osci_sz(:)]), max([ylim_osci_hc(:); ylim_osci_sz(:)])])
    set(ax2_sz,'ylim', [min([ylim_osci_hc(:); ylim_osci_sz(:)]), max([ylim_osci_hc(:); ylim_osci_sz(:)])])

else
    fig1 = [];
end

end

