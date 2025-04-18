
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script performs statistical evaluation of surrogate data testing and
% therefore requires the output of script_01_analyze_BLP_ec_surrogate.m and
% script_01_analyze_BLP_eo_surrogate.m, stored as workspaces for each
% participant in folders results_surrogate_ec and results_surrogate_eo.
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

fnames_ec = dir('results_surrogate_ec/*.mat');
fnames_eo = dir('results_surrogate_eo/*.mat');
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

ns = length(fnames_ec); %...total number of subjects
nch = length(chlab); %......total number of channels
nsur = 100; %...............number of surrogates
nw = 23; %..................total number of sliding windows
n_sz = length(sz_list); %...number of schizophrenia patients
n_hc = length(hc_list); %...number of healthy controls
FDR_alpha = 0.05; %.........alpha for multiple comparisons adjustment

survars = {'true','sur'}; %...........................true vs. surrogate
meanvars = {'mean','var'}; %..........................analysis types
irasa_vars = {'Beta_lo','Beta_hi','Beta_bb'}; %.......frequency ranges
spec_types = {'mixd','frac','osci'}; %................power spectra types
bands = {'delta','theta','alpha','beta','gamma'}; %...frequency bands
freqs = [1, 4; 4, 8; 8, 13; 13, 25; 25, 45]; %........frequency band limits
nrsn = length(rsnlab); %..............................# of RSNs

%% preallocations

% structures for storing structured analysis outputs (ch and RSN level)
struct_hc_ec = struct();
struct_sz_ec = struct();

% BLP analysis
for mv = 1:length(meanvars)
    % channel-level
    struct_hc_ec.(meanvars{mv}).true = zeros(n_hc, nch);
    struct_sz_ec.(meanvars{mv}).true = zeros(n_sz, nch);
    struct_hc_ec.(meanvars{mv}).sur = zeros(n_hc, nch, nsur);
    struct_sz_ec.(meanvars{mv}).sur = zeros(n_sz, nch, nsur);
end


% structures for storing structured analysis outputs (ch and RSN level)
struct_hc_eo = struct();
struct_sz_eo = struct();

% structures for storing statistical testing outputs
struct_stat_beta_eo = struct();
struct_stat_blp_eo = struct();

% BLP analysis
for mv = 1:length(meanvars)
    % channel-level
    struct_hc_eo.(meanvars{mv}).true = zeros(n_hc, nch);
    struct_sz_eo.(meanvars{mv}).true = zeros(n_sz, nch);
    struct_hc_eo.(meanvars{mv}).sur = zeros(n_hc, nch, nsur);
    struct_sz_eo.(meanvars{mv}).sur = zeros(n_sz, nch, nsur);
end

%% load and structure data

% HC
for subj = 1:n_hc % iterate over HC subjects
    load(['results_surrogate_ec/sch_' hc_list{subj} '_ec_surrogate_res.mat'],'results_seg_ec')
    load(['results_surrogate_eo/sch_' hc_list{subj} '_eo_surrogate_res.mat'],'results_seg_eo')


    % getting alpha BLP data
    for mv = 1:length(meanvars)
        struct_hc_ec.(meanvars{mv}).true(subj,:) = results_seg_ec.(meanvars{mv});
        struct_hc_ec.(meanvars{mv}).sur(subj,:,:) = results_seg_ec.([meanvars{mv} '_sur'])';
        struct_hc_eo.(meanvars{mv}).true(subj,:) = results_seg_eo.(meanvars{mv});
        struct_hc_eo.(meanvars{mv}).sur(subj,:,:) = results_seg_eo.([meanvars{mv} '_sur'])';
    end
end

% SZ
for subj = 1:n_sz % iterate over HC subjects
    load(['results_surrogate_ec/sch_' sz_list{subj} '_ec_surrogate_res.mat'],'results_seg_ec')
    load(['results_surrogate_eo/sch_' sz_list{subj} '_eo_surrogate_res.mat'],'results_seg_eo')


    % getting alpha BLP data
    for mv = 1:length(meanvars)
        struct_sz_ec.(meanvars{mv}).true(subj,:) = results_seg_ec.(meanvars{mv});
        struct_sz_ec.(meanvars{mv}).sur(subj,:,:) = results_seg_ec.([meanvars{mv} '_sur'])';
        struct_sz_eo.(meanvars{mv}).true(subj,:) = results_seg_eo.(meanvars{mv});
        struct_sz_eo.(meanvars{mv}).sur(subj,:,:) = results_seg_eo.([meanvars{mv} '_sur'])';
    end
end

% storing channel and RSN labels for convenience
struct_hc_ec.chlab = chlab;
struct_sz_ec.chlab = chlab;
struct_hc_eo.chlab = chlab;
struct_sz_eo.chlab = chlab;


%% run statistics 

p_mean_hc_sur_ec = func_pairwise_FDR(func_pairwise_surrogate(struct_hc_ec, chlab, 'mean'),0.05);
p_mean_sz_sur_ec = func_pairwise_FDR(func_pairwise_surrogate(struct_sz_ec, chlab, 'mean'),0.05);

p_var_hc_sur_ec = func_pairwise_FDR(func_pairwise_surrogate(struct_hc_ec, chlab, 'var'),0.05);
p_var_sz_sur_ec = func_pairwise_FDR(func_pairwise_surrogate(struct_sz_ec, chlab, 'var'),0.05);

p_mean_hc_sur_eo = func_pairwise_FDR(func_pairwise_surrogate(struct_hc_eo, chlab, 'mean'),0.05);
p_mean_sz_sur_eo = func_pairwise_FDR(func_pairwise_surrogate(struct_sz_eo, chlab, 'mean'),0.05);

p_var_hc_sur_eo = func_pairwise_FDR(func_pairwise_surrogate(struct_hc_eo, chlab, 'var'),0.05);
p_var_sz_sur_eo = func_pairwise_FDR(func_pairwise_surrogate(struct_sz_eo, chlab, 'var'),0.05);

%% saving data
save('miscellaneous/ws_results_surrogate_ec.mat',...
    'p_mean_hc_sur_ec', 'p_mean_sz_sur_ec',...
    'p_var_hc_sur_ec', 'p_var_sz_sur_ec')

save('miscellaneous/ws_results_surrogate_eo.mat',...
    'p_mean_hc_sur_eo', 'p_mean_sz_sur_eo',...
    'p_var_hc_sur_eo', 'p_var_sz_sur_eo')

%% plotting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: since the main between-group difference was found in alpha BLP
% variability in EYES CLOSED condition, surrogate data testing is plotted
% only for this case. Other cases can be plotted by uncommenting
% corresponding lines.

% [f1] = func_plot_surrogate_topoplot(p_mean_hc_sur_ec,struct_hc_ec,p_mean_sz_sur_ec,struct_sz_ec,chanlocs,'mean','EC');
% exportgraphics(f1,'figures/fig_mean_alpha_BLP_surrogate.png','resolution',400)
% close(f1)

[f2] = func_plot_surrogate_topoplot(p_var_hc_sur_ec,struct_hc_ec,p_var_sz_sur_ec,struct_sz_ec,chanlocs,'var','EC');
exportgraphics(f2,'figures/fig_var_alpha_BLP_surrogate.png','resolution',400)
close(f2)

% [f3] = func_plot_surrogate_topoplot(p_mean_hc_sur_eo,struct_hc_eo,p_mean_sz_sur_eo,struct_sz_eo,chanlocs,'mean','EO');
% exportgraphics(f3,'figures/fig_mean_alpha_BLP_surrogate_EO.png','resolution',400)
% close(f3)

% [f4] = func_plot_surrogate_topoplot(p_var_hc_sur_eo,struct_hc_eo,p_var_sz_sur_eo,struct_sz_eo,chanlocs,'var','EO');
% exportgraphics(f4,'figures/fig_var_alpha_BLP_surrogate_EO.png','resolution',400)
% close(f4)

