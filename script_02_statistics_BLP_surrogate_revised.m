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
struct_hc_alpha = struct();
struct_sz_alpha = struct();

% BLP analysis
for mv = 1:length(meanvars)
    % channel-level
    struct_hc_alpha.(meanvars{mv}).true = zeros(n_hc, nch);
    struct_sz_alpha.(meanvars{mv}).true = zeros(n_sz, nch);
    struct_hc_alpha.(meanvars{mv}).sur = zeros(n_hc, nch, nsur);
    struct_sz_alpha.(meanvars{mv}).sur = zeros(n_sz, nch, nsur);
end


% structures for storing structured analysis outputs (ch and RSN level)
struct_hc_beta = struct();
struct_sz_beta = struct();


% BLP analysis
for mv = 1:length(meanvars)
    % channel-level
    struct_hc_beta.(meanvars{mv}).true = zeros(n_hc, nch);
    struct_sz_beta.(meanvars{mv}).true = zeros(n_sz, nch);
    struct_hc_beta.(meanvars{mv}).sur = zeros(n_hc, nch, nsur);
    struct_sz_beta.(meanvars{mv}).sur = zeros(n_sz, nch, nsur);
end

%% load and structure data

% HC
for subj = 1:n_hc % iterate over HC subjects
    load(['results_surrogate_ec/sch_' hc_list{subj} '_ec_surrogate_res.mat'],'results_seg_alpha_ec', 'results_seg_beta_ec')

    % getting alpha and beta BLP data
    for mv = 1:length(meanvars)
        struct_hc_alpha.(meanvars{mv}).true(subj,:) = results_seg_alpha_ec.(meanvars{mv});
        struct_hc_alpha.(meanvars{mv}).sur(subj,:,:) = results_seg_alpha_ec.([meanvars{mv} '_sur'])';
        struct_hc_beta.(meanvars{mv}).true(subj,:) = results_seg_beta_ec.(meanvars{mv});
        struct_hc_beta.(meanvars{mv}).sur(subj,:,:) = results_seg_beta_ec.([meanvars{mv} '_sur'])';
    end
end

% SZ
for subj = 1:n_sz % iterate over HC subjects
    load(['results_surrogate_ec/sch_' sz_list{subj} '_ec_surrogate_res.mat'],'results_seg_alpha_ec','results_seg_beta_ec')

    % getting alpha and beta BLP data
    for mv = 1:length(meanvars)
        struct_sz_alpha.(meanvars{mv}).true(subj,:) = results_seg_alpha_ec.(meanvars{mv});
        struct_sz_alpha.(meanvars{mv}).sur(subj,:,:) = results_seg_alpha_ec.([meanvars{mv} '_sur'])';
        struct_sz_beta.(meanvars{mv}).true(subj,:) = results_seg_beta_ec.(meanvars{mv});
        struct_sz_beta.(meanvars{mv}).sur(subj,:,:) = results_seg_beta_ec.([meanvars{mv} '_sur'])';
    end
end

% storing channel and RSN labels for convenience
struct_hc_alpha.chlab = chlab;
struct_sz_alpha.chlab = chlab;
struct_hc_beta.chlab = chlab;
struct_sz_beta.chlab = chlab;

%% save data for confirmatory analysis

struct_hc_alpha_55ch = struct_hc_alpha;
struct_sz_alpha_55ch = struct_sz_alpha;
struct_hc_beta_55ch = struct_hc_beta;
struct_sz_beta_55ch = struct_sz_beta;

save('miscellaneous/ws_data_55ch.mat',...
    'struct_hc_alpha_55ch','struct_sz_alpha_55ch','struct_hc_beta_55ch','struct_sz_beta_55ch')

%% run statistics 

p_mean_hc_sur_alpha = func_pairwise_FDR(func_pairwise_surrogate(struct_hc_alpha, chlab, 'mean'),0.05);
p_mean_sz_sur_alpha = func_pairwise_FDR(func_pairwise_surrogate(struct_sz_alpha, chlab, 'mean'),0.05);

p_var_hc_sur_alpha = func_pairwise_FDR(func_pairwise_surrogate(struct_hc_alpha, chlab, 'var'),0.05);
p_var_sz_sur_alpha = func_pairwise_FDR(func_pairwise_surrogate(struct_sz_alpha, chlab, 'var'),0.05);

p_mean_hc_sur_beta = func_pairwise_FDR(func_pairwise_surrogate(struct_hc_beta, chlab, 'mean'),0.05);
p_mean_sz_sur_beta = func_pairwise_FDR(func_pairwise_surrogate(struct_sz_beta, chlab, 'mean'),0.05);

p_var_hc_sur_beta = func_pairwise_FDR(func_pairwise_surrogate(struct_hc_beta, chlab, 'var'),0.05);
p_var_sz_sur_beta = func_pairwise_FDR(func_pairwise_surrogate(struct_sz_beta, chlab, 'var'),0.05);

%% saving data
save('miscellaneous/ws_results_surrogate_ec_revised.mat',...
    'p_mean_hc_sur_alpha', 'p_mean_sz_sur_alpha',...
    'p_var_hc_sur_alpha', 'p_var_sz_sur_alpha',...
    'p_mean_hc_sur_beta', 'p_mean_sz_sur_beta',...
    'p_var_hc_sur_beta', 'p_var_sz_sur_beta')

%% plotting
% [f1] = func_plot_surrogate_topoplot_revised(p_mean_hc_sur_alpha,struct_hc_alpha,p_mean_sz_sur_alpha,struct_sz_alpha,chanlocs,'alpha','mean','EC');

[f2] = func_plot_surrogate_topoplot_revised(p_var_hc_sur_alpha,struct_hc_alpha,p_var_sz_sur_alpha,struct_sz_alpha,chanlocs,'alpha','var','EC');
% exportgraphics(f2,'figures_revised/fig_var_alpha_BLP_surrogate.png','resolution',400)
% close(f2)

% [f3] = func_plot_surrogate_topoplot_revised(p_mean_hc_sur_beta,struct_hc_beta,p_mean_sz_sur_beta,struct_sz_beta,chanlocs,'beta','mean','EC');

[f4] = func_plot_surrogate_topoplot_revised(p_var_hc_sur_beta,struct_hc_beta,p_var_sz_sur_beta,struct_sz_beta,chanlocs,'beta','var','EC');
% exportgraphics(f4,'figures_revised/supplementary_fig_var_beta_BLP_surrogate.png','resolution',400)
% close(f4)

%% correlation analyses

%% load and structure demographics and PANSS data

load('ws_55ch_montage.mat')
Mr = func_map_to_RSN(M);
rsnlab = Mr.labels; %...............................RSN labels
rsnlab = cat(2,rsnlab,'Global'); %..................add global average
Mr.labels = rsnlab;
Mr.lab = cat(2,Mr.lab,1:length(chlab));

% HC demographics data
dem_tmp_hc = readtable('demographics_table.xlsx','Sheet','HC');
dem_hc = table();
dem_hc.age = dem_tmp_hc.age;
dem_hc.sex = dem_tmp_hc.sex;
dem_hc.sexlab = dem_tmp_hc.label;
dem_hc.edu_years = dem_tmp_hc.edu_years;

% SZ demographics data
dem_tmp_sz = readtable('demographics_table.xlsx','Sheet','SZ');
dem_sz = table();
dem_sz.age = dem_tmp_sz.age;
dem_sz.sex = dem_tmp_sz.sex;
dem_sz.sexlab = dem_tmp_sz.label;
dem_sz.edu_years = dem_tmp_sz.edu_years;
dem_sz.duration = dem_tmp_sz.years;
dem_sz.CPZ = dem_tmp_sz.CPZ;

% PANSS data
data_panss = readtable('demographics_table.xlsx','Sheet','SZ_correlation_analysis');

% define table of confounding variables
data_conf = table();
data_conf.age = data_panss.age;
data_conf.sex = data_panss.sex;
data_conf.edu = data_panss.edu_years;
data_conf.cpz = data_panss.CPZ;
data_conf.dur = data_panss.duration;
z_conf = table2array(data_conf);

%% get EEG spectral data

var_alpha_hc_true = func_aggregate_RSNs_batch(struct_hc_alpha.var.true,Mr);
var_alpha_sz_true = func_aggregate_RSNs_batch(struct_sz_alpha.var.true,Mr);
var_beta_hc_true = func_aggregate_RSNs_batch(struct_hc_beta.var.true,Mr);
var_beta_sz_true = func_aggregate_RSNs_batch(struct_sz_beta.var.true,Mr);

var_alpha_hc_sur = func_aggregate_RSNs_batch(mean(struct_hc_alpha.var.sur,3),Mr);
var_alpha_sz_sur = func_aggregate_RSNs_batch(mean(struct_sz_alpha.var.sur,3),Mr);
var_beta_hc_sur = func_aggregate_RSNs_batch(mean(struct_hc_beta.var.sur,3),Mr);
var_beta_sz_sur = func_aggregate_RSNs_batch(mean(struct_sz_beta.var.sur,3),Mr);

table_alpha_true = func_pairwise_FDR(func_pairwise_RSN_stat(var_alpha_hc_true, var_alpha_sz_true, rsnlab, 'var_alpha_true'),0.05);
table_alpha_sur = func_pairwise_FDR(func_pairwise_RSN_stat(var_alpha_hc_sur, var_alpha_sz_sur, rsnlab, 'var_alpha_sur'),0.05);

table_beta_true = func_pairwise_FDR(func_pairwise_RSN_stat(var_beta_hc_true, var_beta_sz_true, rsnlab, 'var_alpha_true'),0.05);
table_beta_sur = func_pairwise_FDR(func_pairwise_RSN_stat(var_beta_hc_sur, var_beta_sz_sur, rsnlab, 'var_alpha_sur'),0.05);
