
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script plots main results obtained from EYES-CLOSED data after
% running script_02_statistics_BLP_eo.m
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

load('ws_results_BLP_eo.mat')
load('ws_correlations_BLP_eo.mat')
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
[f01, f02, spec_hc_eo, spec_sz_eo, struct_stat_spec_eo] = func_plot_GA_spectra(struct_hc_eo,struct_sz_eo,0,1,1);
exportgraphics(f01,'figures/fig_GA_spectra_eo.png','resolution',400)
close(f01)
% exportgraphics(f02,'figures/fig_GA_matrices_eo.png','resolution',400)
close(f02)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: no significant difference at all between HC and SZ with eyes open!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



