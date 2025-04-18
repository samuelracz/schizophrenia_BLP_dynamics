
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script collect all significant correlations found between BLP
% indices and PANSS scores (Supplementary Table 1).
%
% NOTE: this script utilizes output matlab workspaces that must be
% generated first using script_02_statistics_BLP_ec/eo/vs.
%
% This script evaluates both EYES CLOSED and EYES OPEN data
% 
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 04/18/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('miscellaneous')
load('ws_correlations_BLP_ec.mat')
load('ws_correlations_BLP_eo.mat')
load('ws_correlations_BLP_vs.mat')

table_corr_all = [];

% PANSS-gen, EC
table_tmp = func_parse_panss_table(corr_panss_rsn_blp_gen_ec,'EC','GEN');
table_corr_all = cat(1,table_corr_all,table_tmp);

% PANSS-neg, EC
table_tmp = func_parse_panss_table(corr_panss_rsn_blp_neg_ec,'EC','NEG');
table_corr_all = cat(1,table_corr_all,table_tmp);

% PANSS-pos, EC
table_tmp = func_parse_panss_table(corr_panss_rsn_blp_pos_ec,'EC','POS');
table_corr_all = cat(1,table_corr_all,table_tmp);

% PANSS-sum, EC
table_tmp = func_parse_panss_table(corr_panss_rsn_blp_sum_ec,'EC','SUM');
table_corr_all = cat(1,table_corr_all,table_tmp);

% PANSS-gen, EO
table_tmp = func_parse_panss_table(corr_panss_rsn_blp_gen_eo,'EO','GEN');
table_corr_all = cat(1,table_corr_all,table_tmp);

% PANSS-neg, EO
table_tmp = func_parse_panss_table(corr_panss_rsn_blp_neg_eo,'EO','NEG');
table_corr_all = cat(1,table_corr_all,table_tmp);

% PANSS-pos, EO
table_tmp = func_parse_panss_table(corr_panss_rsn_blp_pos_eo,'EO','POS');
table_corr_all = cat(1,table_corr_all,table_tmp);

% PANSS-sum, EO
table_tmp = func_parse_panss_table(corr_panss_rsn_blp_sum_eo,'EO','SUM');
table_corr_all = cat(1,table_corr_all,table_tmp);

% PANSS-gen, EC vs. EO
table_tmp = func_parse_panss_table(corr_panss_rsn_blp_gen_vs,'dS','GEN');
table_corr_all = cat(1,table_corr_all,table_tmp);

% PANSS-neg, EC vs. EO
table_tmp = func_parse_panss_table(corr_panss_rsn_blp_neg_vs,'dS','NEG');
table_corr_all = cat(1,table_corr_all,table_tmp);

% PANSS-pos, EC vs. EO
table_tmp = func_parse_panss_table(corr_panss_rsn_blp_pos_vs,'dS','POS');
table_corr_all = cat(1,table_corr_all,table_tmp);

% PANSS-sum, EC vs. EO
table_tmp = func_parse_panss_table(corr_panss_rsn_blp_sum_vs,'dS','SUM');
table_corr_all = cat(1,table_corr_all,table_tmp);

table_corr_all = sortrows(table_corr_all,'p','ascend');