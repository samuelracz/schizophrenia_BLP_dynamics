
addpath(genpath('functions'))
addpath('miscellaneous')

fnames_ec = dir('results_phasecycle_ec/*.mat');
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

ns = length(fnames_ec); %...total number of subjects
nch = length(chlab); %......total number of channels
nw = 23; %..................total number of sliding windows
n_sz = length(sz_list); %...number of schizophrenia patients
n_hc = length(hc_list); %...number of healthy controls
FDR_alpha = 0.05; %.........alpha for multiple comparisons adjustment

freq_bins = [(1:44)',(2:45)'];
Nf = size(freq_bins,1); %.....number of frequency estimates
nrsn = length(rsnlab); %..............................# of RSNs

%% load and structure data

mean_hc = zeros(Nf,nch,n_hc);
var_hc = zeros(Nf,nch,n_hc);
mean_sz = zeros(Nf,nch,n_sz);
var_sz = zeros(Nf,nch,n_sz);

% HC
for subj = 1:n_hc % iterate over HC subjects
    load(['results_phasecycle_ec/sch_' hc_list{subj} '_ec_phasecycle_res.mat'])

    mean_hc(:,:,subj) = res_mean;
    var_hc(:,:,subj) = res_var;
end

% SZ
for subj = 1:n_sz % iterate over SZ subjects
    load(['results_phasecycle_ec/sch_' sz_list{subj} '_ec_phasecycle_res.mat'])

    mean_sz(:,:,subj) = res_mean;
    var_sz(:,:,subj) = res_var;
end

%% global statistics

freq = freq_bins(:,1);

mean_hc_glob = squeeze(mean(mean_hc,2))';
var_hc_glob = squeeze(mean(var_hc,2))';

mean_sz_glob = squeeze(mean(mean_sz,2))';
var_sz_glob = squeeze(mean(var_sz,2))';

table_mean = func_pairwise_FDR(func_pairwise_RSN_stat(mean_hc_glob,mean_sz_glob,num2cell(freq),'mean_PSD'),0.05);
table_var = func_pairwise_FDR(func_pairwise_RSN_stat(var_hc_glob,var_sz_glob,num2cell(freq),'mean_PSD'),0.05);

table_mean = renamevars(table_mean,["rsn","rsnID"],["freq","freqID"]);
table_var = renamevars(table_var,["rsn","rsnID"],["freq","freqID"]);

f0 = func_plot_phase_cycle(table_mean,table_var);
% exportgraphics(f0,'figures_revised/supplementary_fig_phase_cycle.png','resolution',400)
% close(f0)
