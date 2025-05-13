
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script plots main results obtained from EYES-CLOSED data after
% running script_02_statistics_BLP_ec.m
%
% Note that EEGLAB (or subfolder 'eeglab_functions') needs to be added to
% the path for generating topoplots.
% 
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 05/13/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% adding dependencies and loading outcomes
addpath(genpath('eeglab_functions'))
addpath(genpath('functions'))
addpath('miscellaneous')

load('ws_results_BLP_ec.mat')
load('ws_correlations_BLP_ec.mat')
load('ws_55ch_montage.mat')

chlab = M.lab;
Mr = func_map_to_RSN(M);
vn = Mr.lab{1};
sm = Mr.lab{2};
da = Mr.lab{3};
val = Mr.lab{4};
fp = Mr.lab{5};
dmn = Mr.lab{6};


%% plot IRASA illustration
[f00, ~, ~] = func_plot_irasa_illustration(struct_hc_ec, struct_sz_ec, 'global',1);
exportgraphics(f00,'figures/fig_IRASA_illustration_ec.png','resolution',300)
close(f00)

%% plot global spectra
[f01, f02, spec_hc_ec, spec_sz_ec, struct_stat_spec_ec] = func_plot_GA_spectra(struct_hc_ec,struct_sz_ec,0,1,1);
exportgraphics(f01,'figures/fig_GA_spectra_ec.png','resolution',300)
close(f01)
% exportgraphics(f02,'figures/fig_GA_matrices_ec.png','resolution',300)
close(f02)

%% plot BLP mean difference in the delta band
[fs1] = func_plot_blp_diff(struct_stat_blp_ec, table_sign_blp_ec, struct_stat_rsn_blp_ec, chanlocs, 'mean', 'delta', 'EC', 1);
exportgraphics(fs1,'figures/fig_mean_delta_BLP_ec.png','resolution',300)
close(fs1)

%% plot BLP mean difference in the theta band
[fs1] = func_plot_blp_diff(struct_stat_blp_ec, table_sign_blp_ec, struct_stat_rsn_blp_ec, chanlocs, 'mean', 'theta', 'EC', 1);
exportgraphics(fs1,'figures/fig_mean_theta_BLP_ec.png','resolution',300)
close(fs1)

%% plot BLP mean difference in the alpha band
[fs1] = func_plot_blp_diff(struct_stat_blp_ec, table_sign_blp_ec, struct_stat_rsn_blp_ec, chanlocs, 'mean', 'alpha', 'EC', 1);
exportgraphics(fs1,'figures/fig_mean_alpha_BLP_ec.png','resolution',300)
close(fs1)

%% plot BLP variance difference in the alpha band
[fs1] = func_plot_blp_diff(struct_stat_blp_ec, table_sign_blp_ec, struct_stat_rsn_blp_ec, chanlocs, 'var', 'alpha', 'EC', 1);
exportgraphics(fs1,'figures/fig_var_alpha_BLP_ec.png','resolution',300)
close(fs1)

%% plot BLP variance difference in the beta band
[fs1] = func_plot_blp_diff(struct_stat_blp_ec, table_sign_blp_ec, struct_stat_rsn_blp_ec, chanlocs, 'var', 'beta', 'EC',1);
exportgraphics(fs1,'figures/fig_var_beta_BLP_ec.png','resolution',300)
close(fs1)

%% combined figure
fc = figure('color','w','units','normalized','outerposition',[0.1,0.2,0.8,0.7]);

subplot(121)
sc3 = scatter(corr_panss_rsn_blp_neg_ec.res_eeg{1},corr_panss_rsn_blp_neg_ec.res_panss{1},140,'LineWidth',1.5,'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceColor',[0.7 0.9 1]);
ls = lsline;
ls.LineWidth = 3;
ls.Color = [0.8 0 0];
set(gca,'LineWidth',1,'FontSize',24,'box','on')
% annotations
a_lo = annotation('textbox');
a_lo.String = {['{\itr}=' num2str(corr_panss_rsn_blp_neg_ec.r(1),'%.4f')],['{\itp}=' num2str(corr_panss_rsn_blp_neg_ec.p(1),'%.4f')]};
a_lo.FontSize = 28;
a_lo.Position = [0.37, 0.7, 0.2, 0.2];
a_lo.LineStyle = 'none';
xlabel('Dorsal attention network \sigma(\beta-BLP^{EC}_{\it{mixd}})', 'FontSize', 28, 'visible', 'on')
ylabel('PANSS - Negative','FontSize',28,'visible','on')
title('PANSS-NEG ~ \sigma(\beta-BLP^{EC}_{\it{mixd}}) in DA','FontSize',28)

subplot(122)
sc4 = scatter(corr_panss_rsn_blp_sum_ec.res_eeg{1},corr_panss_rsn_blp_sum_ec.res_panss{1},140,'LineWidth',1.5,'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceColor',[0.7 0.9 1]);
ls = lsline;
ls.LineWidth = 3;
ls.Color = [0.8 0 0];
set(gca,'LineWidth',1,'FontSize',24,'box','on')
% annotations
a_lo = annotation('textbox');
a_lo.String = {['{\itr}=' num2str(corr_panss_rsn_blp_sum_ec.r(1),'%.4f')],['{\itp}=' num2str(corr_panss_rsn_blp_sum_ec.p(1),'%.4f')]};
a_lo.FontSize = 28;
a_lo.Position = [0.81, 0.7, 0.2, 0.2];
a_lo.LineStyle = 'none';
xlabel('Dorsal attention network \sigma(\beta-BLP^{EC}_{\it{mixd}})', 'FontSize', 28, 'visible', 'on')
ylabel('PANSS - Total','FontSize',28,'visible','on')
title('PANSS-SUM \sigma(\beta-BLP^{EC}_{\it{mixd}}) in DA','FontSize',28)

exportgraphics(fc,'fig_corr_panss_res_all.png','resolution',300)
close(fc)

