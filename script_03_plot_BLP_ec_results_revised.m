%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note that EEGLAB (or subfolder 'eeglab_functions') needs to be added to
% the path for generating topoplots.
%
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 12/28/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% adding dependencies and loading outcomes
addpath(genpath('eeglab_functions'))
addpath(genpath('functions'))
addpath('miscellaneous')

load('ws_results_BLP_ec_revised.mat')
load('ws_correlations_BLP_ec_revised.mat')
load('ws_55ch_montage.mat')

chlab = M.lab;
Mr = func_map_to_RSN(M);
vn = Mr.lab{1};
sm = Mr.lab{2};
da = Mr.lab{3};
val = Mr.lab{4};
fp = Mr.lab{5};
dmn = Mr.lab{6};

table_sign_blp_ec = table_blp_corr_ec(table_blp_corr_ec.h_FDR==1,:);


%% plot IRASA illustration
[f00, ~, ~] = func_plot_irasa_illustration(struct_hc_ec, struct_sz_ec, 'global',1);
exportgraphics(f00,'figures_revised/fig_IRASA_illustration_ec.png','resolution',400)
close(f00)

%% plot global spectra
[f01, f02, spec_hc_ec, spec_sz_ec, struct_stat_spec_ec] = func_plot_GA_spectra(struct_hc_ec,struct_sz_ec,0,1,0);
exportgraphics(f01,'figures_revised/fig_GA_spectra_ec.png','resolution',400)
close(f01)
exportgraphics(f02,'figures_revised/fig_GA_matrices_ec.png','resolution',400)
close(f02)


%% plot differences in total fractal power
[fs1] = func_plot_frac_diff_revised(table_all_corr_ec, table_all_rsn_corr_ec, chanlocs, 'EC', 1);
exportgraphics(fs1,'figures_revised/fig_TFP_ec.png','resolution',400)
close(fs1)

%% plot BLP variance difference in the alpha band
[fs1] = func_plot_blp_diff_revised(table_all_corr_ec, table_all_rsn_corr_ec, chanlocs, 'var', 'alpha', 'EC', 1);
exportgraphics(fs1,'figures_revised/fig_var_alpha_BLP_ec.png','resolution',400)
close(fs1)

%% plot BLP variance difference in the beta band
[fs1] = func_plot_blp_diff_revised(table_all_corr_ec, table_all_rsn_corr_ec, chanlocs, 'var', 'beta', 'EC', 1);
exportgraphics(fs1,'figures_revised/fig_var_beta_BLP_ec.png','resolution',400)
close(fs1)

%% plot spectral slope correlation, RSN-wise - DA and PANSS-NEG
corr_panns_rsn_blp = sortrows(table_panss_rsn_blp_all_ec,'p','ascend');

f5 = figure('color','w','units','normalized','outerposition',[0.25,0.1,0.5,0.7]);
sc1 = scatter(corr_panns_rsn_blp.res_eeg{1},corr_panns_rsn_blp.res_panss{1},140,'LineWidth',1.5,'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceColor',[0.7 0.9 1]);
ls = lsline;
ls.LineWidth = 3;
ls.Color = [0.8 0 0];
set(gca,'LineWidth',1,'FontSize',24,'box','on')
% annotations
a_lo = annotation('textbox');
a_lo.String = {['{\itr}=' num2str(corr_panns_rsn_blp.r(1),'%.4f')],['{\itp}=' num2str(corr_panns_rsn_blp.p(1),'%.4f')]};
a_lo.FontSize = 28;
a_lo.Position = [0.75, 0.7, 0.2, 0.2];
a_lo.LineStyle = 'none';
xlabel('Dorsal attention network \sigma(\betaBLP^{EC}_{\it{mixd}})', 'FontSize', 28, 'visible', 'on')
ylabel('PANSS - Negative','FontSize',28,'visible','on')
title('PANSS-NEG ~ \sigma(\betaBLP^{EC}_{\it{mixd}}) in DA','FontSize',28)
% exportgraphics(f5,'figures_revised/fig_corr_panss_rsn_neg_ec.png','resolution',400)
% close(f5)

%% plot spectral slope correlation, RSN-wise - DA and PANSS-SUM
f5 = figure('color','w','units','normalized','outerposition',[0.25,0.1,0.5,0.7]);
sc2 = scatter(corr_panns_rsn_blp.res_eeg{2},corr_panns_rsn_blp.res_panss{2},140,'LineWidth',1.5,'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceColor',[0.7 0.9 1]);
ls = lsline;
ls.LineWidth = 3;
ls.Color = [0.8 0 0];
set(gca,'LineWidth',1,'FontSize',24,'box','on')
% annotations
a_lo = annotation('textbox');
a_lo.String = {['{\itr}=' num2str(corr_panns_rsn_blp.r(2),'%.4f')],['{\itp}=' num2str(corr_panns_rsn_blp.p(2),'%.4f')]};
a_lo.FontSize = 28;
a_lo.Position = [0.75, 0.7, 0.2, 0.2];
a_lo.LineStyle = 'none';
xlabel('Dorsal attention network \sigma(\betaBLP^{EC}_{\it{mixd}})', 'FontSize', 28, 'visible', 'on')
ylabel('PANSS - Total','FontSize',28,'visible','on')
title('PANSS-SUM ~ \sigma(\betaBLP^{EC}_{\it{mixd}}) in DA','FontSize',28)
% exportgraphics(f5,'figures_revised/fig_corr_panss_rsn_sum_ec.png','resolution',400)
% close(f5)


