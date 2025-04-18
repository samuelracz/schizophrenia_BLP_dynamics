
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script evaluates if CPZ equivalent dose is correlated to BLP indices
% obtained in the SZ group. In turn, this script should be run after all
% other statistical evaluations are completed and the workspaces
% ws_results_BLP_ec.mat, ws_results_BLP_eo.mat and ws_results_BLP_vs.mat
% exist.
%
% note that that is stored such that the first 30 subjects (sch_002 - 
% sch_050) are SZ patients, the remaining 31 subjects (sch_101 - sch_333)
% are healthy controls.
%
% NOTE: this script utilizes output matlab workspaces that must be
% generated first using script_01_analyze_BLP_ec_surrogate.m and
% script_01_analyze_BLP_eo_surrogate.m
%
% This script evaluates both EYES CLOSED and EYES OPEN data
% 
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 04/18/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('functions'))
addpath('miscellaneous')

load('ws_results_BLP_ec.mat')
load('ws_results_BLP_eo.mat')
load('ws_results_BLP_vs.mat')

load('ws_55ch_montage.mat')
load('ws_55ch_chanlocs.mat')

subj_list = {'002','003','004','005','006','007','008','009','013','014','016','018','021','022','023','026','027','029','030','031','033','034','035','036','039','040','044','045','046','050',...
    '101','102','103','104','105','106','107','108','109','111','113','115','116','117','118','120','121','122','123','124','125','126','127','128','129','131','136','137','138','140','333'};

hc_list = {'101','102','103','104','105','106','107','108','109','111','113','115','116','117','118','120','121','122','123','124','125','126','127','128','129','131','136','137','138','140','333'};
sz_list = {'002','003','004','005','006','007','008','009','013','014','016','018','021','022','023','026','027','029','030','031','033','034','035','036','039','040','044','045','046','050'};

%% analysis constants

% map 55-channel EEG to resting-state networks (RSN)
Mr = func_map_to_RSN(M);

chlab = M.lab; %....................................channel labels
rsnlab = Mr.labels; %...............................RSN labels
rsnlab = cat(2,rsnlab,'Global'); %..................add global average
Mr.labels = rsnlab;
Mr.lab = cat(2,Mr.lab,1:length(chlab));

ns = length(subj_list); %...total number of subjects
nch = length(chlab); %......total number of channels
nw = 23; %..................total number of sliding windows
n_sz = length(sz_list); %...number of schizophrenia patients
n_hc = length(hc_list); %...number of healthy controls
FDR_alpha = 0.05; %.........alpha for multiple comparisons adjustment

meanvars = {'mean','var'}; %..........................analysis types
spec_types = {'mixd','frac','osci'}; %................power spectra types
bands = {'delta','theta','alpha','beta','gamma'}; %...frequency bands
freqs = [1, 4; 4, 8; 8, 13; 13, 25; 25, 45]; %........frequency band limits
nrsn = length(rsnlab); %..............................# of RSNs

%% load and structure demographics and PANSS data

% HC demographics data
table_tmp_hc = readtable('demographics_table.xlsx','Sheet','HC');
table_dem_hc = table();
table_dem_hc.age = table_tmp_hc.age;
table_dem_hc.sex = table_tmp_hc.sex;
table_dem_hc.sexlab = table_tmp_hc.label;
table_dem_hc.edu_years = table_tmp_hc.edu_years;

% SZ demographics data
table_tmp_sz = readtable('demographics_table.xlsx','Sheet','SZ');
table_dem_sz = table();
table_dem_sz.age = table_tmp_sz.age;
table_dem_sz.sex = table_tmp_sz.sex;
table_dem_sz.sexlab = table_tmp_sz.label;
table_dem_sz.edu_years = table_tmp_sz.edu_years;
table_dem_sz.duration = table_tmp_sz.years;
table_dem_sz.CPZ = table_tmp_sz.CPZ;

% PANSS data
table_panss = readtable('demographics_table.xlsx','Sheet','SZ_correlation_analysis');

% define table of confounding variables - not CPZ!
table_conf = table();
table_conf.age = table_panss.age;
table_conf.sex = table_panss.sex;
table_conf.edu = table_panss.edu_years;
table_conf.cpz = table_panss.CPZ;
table_conf.dur = table_panss.duration;
z_conf = table2array(table_conf);

%% compute correlations with CPZ equivalent dose
chlab = struct_sz_ec.chlab;
rsnlab = struct_sz_rsn_ec.rsnlab;

corr_CPZ_blp_rsn_ec = func_corr_RSN_CPZ(table_sign_rsn_blp_ec, table_conf);
corr_CPZ_blp_rsn_vs = func_corr_RSN_CPZ(table_sign_rsn_blp_vs, table_conf);

corr_CPZ_blp_rsn_ec = corr_CPZ_blp_rsn_ec(contains(corr_CPZ_blp_rsn_ec.rsn,'Global'),:);
corr_CPZ_blp_rsn_vs = corr_CPZ_blp_rsn_vs(contains(corr_CPZ_blp_rsn_vs.rsn,'Global'),:);

%% compute correlations between CPZ and PANSS

corr_CPZ_PANSS_GEN = func_corr_PANSS_CPZ(table_panss,'GEN',table_conf,1);
corr_CPZ_PANSS_NEG = func_corr_PANSS_CPZ(table_panss,'NEG',table_conf,1);
corr_CPZ_PANSS_POS = func_corr_PANSS_CPZ(table_panss,'POS',table_conf,1);
corr_CPZ_PANSS_SUM = func_corr_PANSS_CPZ(table_panss,'SUM',table_conf,1);

corr_CPZ_PANSS_all = cat(1,corr_CPZ_PANSS_GEN,cat(1,corr_CPZ_PANSS_NEG,cat(1,corr_CPZ_PANSS_POS,corr_CPZ_PANSS_SUM)));

%% save results

save('miscellaneous/ws_results_CPZ_corr.mat','corr_CPZ_blp_rsn_ec','corr_CPZ_blp_rsn_vs','corr_CPZ_PANSS_all');

%% plot results

f = func_plot_CPZ_PANSS_corr(corr_CPZ_PANSS_all);
exportgraphics(f,'figures/fig_corr_cpz_panss_all.png','resolution',400)
close(f)
