
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script plots main results obtained from EYES-OPEN vs. EYES CLOSED
% difference (reactivity) indices after running 
% script_02_statistics_BLP_vs_eo_ec.m and
% script_02_statistics_BLP_state_response.m
%
% Note that EEGLAB (or subfolder 'eeglab_functions') needs to be added to
% the path for generating topoplots.
% 
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 04/18/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% adding dependencies and loading outcomes
addpath(genpath('eeglab_functions'))
addpath(genpath('functions'))
addpath('miscellaneous')

load('ws_results_BLP_vs.mat')
load('ws_results_BLP_response.mat')
load('ws_correlations_BLP_ec.mat')
load('ws_correlations_BLP_vs.mat')
load('ws_55ch_montage.mat')

chlab = M.lab;
Mr = func_map_to_RSN(M);
vn = Mr.lab{1};
sm = Mr.lab{2};
da = Mr.lab{3};
val = Mr.lab{4};
fp = Mr.lab{5};
dmn = Mr.lab{6};


%% plot global spectra
[f01, spec_hc_vs, spec_sz_vs, struct_stat_spec_vs] = func_plot_GA_vs_spectra(struct_hc_vs,struct_sz_vs,0,1);
exportgraphics(f01,'figures/fig_GA_spectra_vs.png','resolution',400)
close(f01)

%% effect of EC-EO transition in the alpha band
[f03] = func_plot_EC_EO_BLP_response(struct_hc_stat_res, struct_sz_stat_res, chanlocs, 'mean', 'alpha', 1);
exportgraphics(f03,'figures/fig_mean_alpha_BLP_response.png','resolution',400)
close(f03)

%% effect of EC-EO transition in the gamma band
[fs02] = func_plot_EC_EO_BLP_response(struct_hc_stat_res, struct_sz_stat_res, chanlocs, 'mean', 'gamma', 1);
exportgraphics(fs02,'figures/fig_mean_gamma_BLP_response.png','resolution',400)
close(fs02)

%% plot BLP variance difference in the beta band
[fs1] = func_plot_blp_diff(struct_stat_blp_vs, table_sign_blp_vs, struct_stat_rsn_blp_vs, chanlocs, 'var', 'beta', 'dS', 1);
exportgraphics(fs1,'figures/fig_var_beta_BLP_vs.png','resolution',400)
close(fs1)


%% combined figure
f6 = figure('color','w','units','normalized','outerposition',[0,0.2,1,0.6]);
subplot(131)
sc3 = scatter(corr_panss_rsn_blp_neg_ec.res_eeg{1},corr_panss_rsn_blp_neg_ec.res_panss{1},140,'LineWidth',1.5,'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceColor',[0.7 0.9 1]);
ls = lsline;
ls.LineWidth = 3;
ls.Color = [0.8 0 0];
set(gca,'LineWidth',1,'FontSize',24,'box','on')
% annotations
a_lo = annotation('textbox');
a_lo.String = {['{\itr}=' num2str(corr_panss_rsn_blp_neg_ec.r(1),'%.4f')],['{\itp}=' num2str(corr_panss_rsn_blp_neg_ec.p(1),'%.4f')]};
a_lo.FontSize = 28;
a_lo.Position = [0.27, 0.7, 0.2, 0.2];
a_lo.LineStyle = 'none';
xlabel('Dorsal attention network \sigma(\betaBLP^{EC}_{\it{mixd}})', 'FontSize', 28, 'visible', 'on')
ylabel('PANSS - Negative','FontSize',28,'visible','on')
title('PANSS-NEG ~ \sigma(\betaBLP^{EC}_{\it{mixd}}) in DA','FontSize',28)

subplot(132)
sc4 = scatter(corr_panss_rsn_blp_gen_vs.res_eeg{1},corr_panss_rsn_blp_gen_vs.res_panss{1},140,'LineWidth',1.5,'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceColor',[0.7 0.9 1]);
ls = lsline;
ls.LineWidth = 3;
ls.Color = [0.8 0 0];
set(gca,'LineWidth',1,'FontSize',24,'box','on')
% annotations
a_lo = annotation('textbox');
a_lo.String = {['{\itr*}=' num2str(corr_panss_rsn_blp_gen_vs.r(1),'%.4f')],['{\itp*}=' num2str(corr_panss_rsn_blp_gen_vs.p(1),'%.4f')]};
a_lo.FontSize = 28;
a_lo.Position = [0.54, 0.7, 0.2, 0.2];
a_lo.LineStyle = 'none';
xlabel('Somatomotor network \sigma(\betaBLP^{dS}_{\it{osci}})', 'FontSize', 28, 'visible', 'on')
ylabel('PANSS - General','FontSize',28,'visible','on')
title('PANSS-GEN \sigma(\betaBLP^{dS}_{\it{osci}})','FontSize',28)

subplot(133)
sc5 = scatter(corr_panss_rsn_blp_pos_vs.res_eeg{1},corr_panss_rsn_blp_pos_vs.res_panss{1},140,'LineWidth',1.5,'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceColor',[0.7 0.9 1]);
ls = lsline;
ls.LineWidth = 3;
ls.Color = [0.8 0 0];
set(gca,'LineWidth',1,'FontSize',24,'box','on')
% annotations
a_lo = annotation('textbox');
a_lo.String = {['{\itr*}=' num2str(corr_panss_rsn_blp_pos_vs.r(1),'%.4f')],['{\itp*}=' num2str(corr_panss_rsn_blp_pos_vs.p(1),'%.4f')]};
a_lo.FontSize = 28;
a_lo.Position = [0.82, 0.7, 0.2, 0.2];
a_lo.LineStyle = 'none';
xlabel('Dorsal attention network \sigma(\betaBLP^{dS}_{\it{frac}})', 'FontSize', 28, 'visible', 'on')
ylabel('PANSS - Positive','FontSize',28,'visible','on')
title('PANSS-POS ~ \sigma(\betaBLP^{dS}_{\it{frac}})','FontSize',28)

exportgraphics(f6,'figures/fig_corr_panss_rsn_all.png','resolution',400)
close(f6)

