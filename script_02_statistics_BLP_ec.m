
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script first loads, structures and aggregates the outputs of the
% IRASA analysis, then performs the following statistical evaluations:
% - mixed, fractal and oscillatory BLP mean over time
% - mixed, fractal and oscillatory BLP variance over time
% - illustration of global EYES-CLSOED EEG power spectra in HC and SZ
% - demographics comparison between HC and SZ
% - correlation analysis between findings and PANSS scores in SZ
%
%
% note that that is stored such that the first 30 subjects (sch_002 - 
% sch_050) are SZ patients, the remaining 31 subjects (sch_101 - sch_333)
% are healthy controls.
%
% NOTE: this script utilizes output matlab workspaces that must be
% generated first using script_01_analyze_BLP_ec.m
%
% This script only evaluates EYES CLOSED data
% 
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 04/18/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('functions'))
addpath('miscellaneous')

fnames_ec = dir('results_BLP_ec/*.mat');
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

meanvars = {'mean','var'}; %..........................analysis types
irasa_vars = {'Beta_lo','Beta_hi','Beta_bb'}; %.......slope ranges (only for exploratory purposes)
spec_types = {'mixd','frac','osci'}; %................power spectra types
bands = {'delta','theta','alpha','beta','gamma'}; %...frequency bands
freqs = [1, 4; 4, 8; 8, 13; 13, 25; 25, 45]; %........frequency band limits
nrsn = length(rsnlab); %..............................# of RSNs

%% preallocations

% structures for storing structured analysis outputs (ch and RSN level)
struct_hc_ec = struct();
struct_sz_ec = struct();
struct_hc_rsn_ec = struct();
struct_sz_rsn_ec = struct();

% structures for storing statistical testing outputs
struct_stat_beta_ec = struct();
struct_stat_blp_ec = struct();
struct_stat_rsn_beta_ec = struct();
struct_stat_rsn_blp_ec = struct();

% IRASA variables
for iv = 1:length(irasa_vars)
    struct_hc_ec.(irasa_vars{iv}) = struct();
    struct_sz_ec.(irasa_vars{iv}) = struct();
    for mv = 1:length(meanvars)
        % channel-level
        struct_hc_ec.(irasa_vars{iv}).(meanvars{mv}) = zeros(n_hc, nch);
        struct_sz_ec.(irasa_vars{iv}).(meanvars{mv}) = zeros(n_sz, nch);
        struct_stat_beta_ec.(irasa_vars{iv}).(meanvars{mv}) = [];
        % RSN-level
        struct_hc_rsn_ec.(irasa_vars{iv}).(meanvars{mv}) = zeros(n_hc, nrsn);
        struct_sz_rsn_ec.(irasa_vars{iv}).(meanvars{mv}) = zeros(n_sz, nrsn);
        struct_stat_rsn_beta_ec.(irasa_vars{iv}).(meanvars{mv}) = [];
    end
end

% BLP variables
for st = 1:length(spec_types)
    struct_hc_ec.(spec_types{st}) = struct();
    struct_sz_ec.(spec_types{st}) = struct();
    for b = 1:length(bands)
        struct_hc_ec.(spec_types{st}).(bands{b}) = struct();
        struct_sz_ec.(spec_types{st}).(bands{b}) = struct();
        for mv = 1:length(meanvars)
            % channel-level
            struct_hc_ec.(spec_types{st}).(bands{b}).(meanvars{mv}) = zeros(n_hc, nch);
            struct_sz_ec.(spec_types{st}).(bands{b}).(meanvars{mv}) = zeros(n_sz, nch);
            struct_stat_blp_ec.(spec_types{st}).(bands{b}).(meanvars{mv}) = [];
            % RSN-level
            struct_hc_rsn_ec.(spec_types{st}).(bands{b}).(meanvars{mv}) = zeros(n_hc, nrsn);
            struct_sz_rsn_ec.(spec_types{st}).(bands{b}).(meanvars{mv}) = zeros(n_sz, nrsn);
            struct_stat_rsn_blp_ec.(spec_types{st}).(bands{b}).(meanvars{mv}) = [];
        end
    end
end

% for global spectra illustration
struct_hc_ec.spec = struct('subjID',[],'plaw',[]);
struct_hc_rsn_ec.spec = struct('subjID',[],'plaw',[]);
struct_sz_ec.spec = struct('subjID',[],'plaw',[]);
struct_sz_rsn_ec.spec = struct('subjID',[],'plaw',[]);

%% load and structure data

% HC
for subj = 1:n_hc % iterate over HC subjects
    load(['results_BLP_ec/sch_' hc_list{subj} '_ec_BLP_res.mat'])

    % IRASA analysis
    for iv = 1:length(irasa_vars) % iterate over fitting ranges
        v_tmp = zeros(nw, nch);
        for t = 1:nw % iterate over sliding windows
            v_tmp(t,:) = results_seg_ec(t).plaw.(irasa_vars{iv});
        end

        % aggregate over resting-state networks
        v_tmp_rsn = func_aggregate_RSNs_batch(v_tmp, Mr);

        struct_hc_ec.(irasa_vars{iv}).mean(subj,:) = mean(v_tmp,1);
        struct_hc_ec.(irasa_vars{iv}).var(subj,:) = std(v_tmp,[],1);
        struct_hc_rsn_ec.(irasa_vars{iv}).mean(subj,:) = mean(v_tmp_rsn,1);
        struct_hc_rsn_ec.(irasa_vars{iv}).var(subj,:) = std(v_tmp_rsn,[],1);
    end

    % BLP analysis
    for st = 1:length(spec_types)
        for b = 1:length(bands)
            v_tmp = zeros(nw, nch);
            for t = 1:nw
                freq = results_seg_ec(t).plaw.freq;
                if strcmp(spec_types{st},'mixd')
                    psd = results_seg_ec(t).plaw.mixd;
                    v_tmp(t,:) = log(sum(psd(freq>=freqs(b,1) & freq<freqs(b,2),:),1));
                elseif strcmp(spec_types{st},'frac')
                    psd = results_seg_ec(t).plaw.frac;
                    v_tmp(t,:) = log(sum(psd(freq>=freqs(b,1) & freq<freqs(b,2),:),1));
                elseif strcmp(spec_types{st},'osci')
                    psdm = results_seg_ec(t).plaw.mixd;
                    psdf = results_seg_ec(t).plaw.frac;
                    t_tmpm = log(sum(psdm(freq>=freqs(b,1) & freq<freqs(b,2),:),1));
                    t_tmpf = log(sum(psdf(freq>=freqs(b,1) & freq<freqs(b,2),:),1));
                    v_tmp(t,:) = t_tmpm - t_tmpf;
                end
            end

            % aggregate over resting-state networks
            v_tmp_rsn = func_aggregate_RSNs_batch(v_tmp, Mr);

            struct_hc_ec.(spec_types{st}).(bands{b}).mean(subj,:) = mean(v_tmp,1);
            struct_hc_ec.(spec_types{st}).(bands{b}).var(subj,:) = std(v_tmp,[],1);
            struct_hc_rsn_ec.(spec_types{st}).(bands{b}).mean(subj,:) = mean(v_tmp_rsn,1);
            struct_hc_rsn_ec.(spec_types{st}).(bands{b}).var(subj,:) = std(v_tmp_rsn,[],1);
        end
    end

    % store data, compute average power spectrum for RSNs (illustration)
    struct_hc_ec.spec(subj).subjID = fnames_ec(subj).name;
    struct_hc_rsn_ec.spec(subj).subjID = fnames_ec(subj).name;

    % compute average spectra (over sliding windows)
    [plaw, plaw_rsn] = func_average_RSN_results_irasa(results_seg_ec, chlab, Mr);
    struct_hc_ec.spec(subj).plaw = plaw;
    struct_hc_rsn_ec.spec(subj).plaw = plaw_rsn;
end

% SZ
for subj = 1:n_sz % iterate over SZ subjects
    load(['results_BLP_ec/sch_' sz_list{subj} '_ec_BLP_res.mat'])
    for iv = 1:length(irasa_vars) % iterate over fitting ranges
        v_tmp = zeros(nw, nch);
        for t = 1:nw % iterate over sliding windows
            v_tmp(t,:) = results_seg_ec(t).plaw.(irasa_vars{iv});
        end

        % aggregate over resting-state networks
        v_tmp_rsn = func_aggregate_RSNs_batch(v_tmp, Mr);

        struct_sz_ec.(irasa_vars{iv}).mean(subj,:) = mean(v_tmp,1);
        struct_sz_ec.(irasa_vars{iv}).var(subj,:) = std(v_tmp,[],1);
        struct_sz_rsn_ec.(irasa_vars{iv}).mean(subj,:) = mean(v_tmp_rsn,1);
        struct_sz_rsn_ec.(irasa_vars{iv}).var(subj,:) = std(v_tmp_rsn,[],1);
    end

    % BLP analysis
    for st = 1:length(spec_types)
        for b = 1:length(bands)
            v_tmp = zeros(nw, nch);
            for t = 1:nw
                freq = results_seg_ec(t).plaw.freq;
                if strcmp(spec_types{st},'mixd')
                    psd = results_seg_ec(t).plaw.mixd;
                    v_tmp(t,:) = log(sum(psd(freq>=freqs(b,1) & freq<freqs(b,2),:),1));
                elseif strcmp(spec_types{st},'frac')
                    psd = results_seg_ec(t).plaw.frac;
                    v_tmp(t,:) = log(sum(psd(freq>=freqs(b,1) & freq<freqs(b,2),:),1));
                elseif strcmp(spec_types{st},'osci')
                    psdm = results_seg_ec(t).plaw.mixd;
                    psdf = results_seg_ec(t).plaw.frac;
                    t_tmpm = log(sum(psdm(freq>=freqs(b,1) & freq<freqs(b,2),:),1));
                    t_tmpf = log(sum(psdf(freq>=freqs(b,1) & freq<freqs(b,2),:),1));
                    v_tmp(t,:) = t_tmpm - t_tmpf;
                end
            end

            % aggregate over resting-state networks
            v_tmp_rsn = func_aggregate_RSNs_batch(v_tmp, Mr);

            struct_sz_ec.(spec_types{st}).(bands{b}).mean(subj,:) = mean(v_tmp,1);
            struct_sz_ec.(spec_types{st}).(bands{b}).var(subj,:) = std(v_tmp,[],1);
            struct_sz_rsn_ec.(spec_types{st}).(bands{b}).mean(subj,:) = mean(v_tmp_rsn,1);
            struct_sz_rsn_ec.(spec_types{st}).(bands{b}).var(subj,:) = std(v_tmp_rsn,[],1);
        end
    end

    % store data, compute average power spectrum for RSNs (illustration)
    struct_sz_ec.spec(subj).subjID = fnames_ec(subj).name;
    struct_sz_rsn_ec.spec(subj).subjID = fnames_ec(subj).name;

    % compute average spectra (over sliding windows)
    [plaw, plaw_rsn] = func_average_RSN_results_irasa(results_seg_ec, chlab, Mr);
    struct_sz_ec.spec(subj).plaw = plaw;
    struct_sz_rsn_ec.spec(subj).plaw = plaw_rsn;
end

% storing channel and RSN labels for convenience
struct_hc_ec.chlab = chlab;
struct_sz_ec.chlab = chlab;
struct_hc_rsn_ec.rsnlab = rsnlab;
struct_sz_rsn_ec.rsnlab = rsnlab;


%% run statistics on BLP - channel-wise (exploratory)

% preallocate
table_sign_blp_ec = [];
table_sign_blp_unadjusted_ec = [];

for st = 1:length(spec_types) % iterate over spectra types
    for b = 1:length(bands) % iterate over frequency bands
        for mv = 1:length(meanvars) % iterate over mean and variance

            % selecting variables
            v_hc = struct_hc_ec.(spec_types{st}).(bands{b}).(meanvars{mv});
            v_sz = struct_sz_ec.(spec_types{st}).(bands{b}).(meanvars{mv});
            vname = [meanvars{mv} '_' bands{b} '_' spec_types{st}];

            % statistical testing
            p_corr = func_pairwise_FDR(func_pairwise_channel_stat(v_hc, v_sz, chlab, vname), FDR_alpha);

            % store outcomes
            struct_stat_blp_ec.(spec_types{st}).(bands{b}).(meanvars{mv}) = p_corr;
        end
    end
end

% collect significant differences
for st = 1:length(spec_types)
    for b = 1:length(bands)
        for mv = 1:length(meanvars)
            table_corr = struct_stat_blp_ec.(spec_types{st}).(bands{b}).(meanvars{mv});
            tmp_sign = table_corr(table_corr.h==1,:);
            tmp_sign_FDR = table_corr(table_corr.h_FDR==1,:);

            % unadjusted (exploratory)
            if isempty(tmp_sign)
                continue
            else
                if isempty(table_sign_blp_unadjusted_ec)
                    table_sign_blp_unadjusted_ec = tmp_sign;
                else
                    table_sign_blp_unadjusted_ec = cat(1,table_sign_blp_unadjusted_ec,tmp_sign);
                end
            end

            % adjusted
            if isempty(tmp_sign_FDR)
                continue
            else
                if isempty(table_sign_blp_ec)
                    table_sign_blp_ec = tmp_sign_FDR;
                else
                    table_sign_blp_ec = cat(1,table_sign_blp_ec,tmp_sign_FDR);
                end
            end
        end
    end
end


%% run BLP statistics - RSN-wise

% preallocate
table_sign_rsn_blp_ec = [];
table_sign_rsn_blp_unadjusted_ec = [];

for st = 1:length(spec_types) % iterate over spectra types
    for b = 1:length(bands) % iterate over frequency bands
        for mv = 1:length(meanvars) % iterate over mean and variance

            % select variables
            v_hc = struct_hc_rsn_ec.(spec_types{st}).(bands{b}).(meanvars{mv});
            v_sz = struct_sz_rsn_ec.(spec_types{st}).(bands{b}).(meanvars{mv});
            vname = [meanvars{mv} '_' bands{b} '_' spec_types{st}];


            % statistical testing
            p_corr = func_pairwise_FDR(func_pairwise_RSN_stat(v_hc, v_sz, rsnlab, vname), FDR_alpha);

            % store outcomes
            struct_stat_rsn_blp_ec.(spec_types{st}).(bands{b}).(meanvars{mv}) = p_corr;
        end
    end
end

% collect significant differences
for st = 1:length(spec_types)
    for b = 1:length(bands)
        for mv = 1:length(meanvars)
            table_corr = struct_stat_rsn_blp_ec.(spec_types{st}).(bands{b}).(meanvars{mv});
            tmp_sign = table_corr(table_corr.h==1,:);
            tmp_sign_FDR = table_corr(table_corr.h_FDR==1,:);

            % unadjusted (exploratory)
            if isempty(tmp_sign)
                continue
            else
                if isempty(table_sign_rsn_blp_unadjusted_ec)
                    table_sign_rsn_blp_unadjusted_ec = tmp_sign;
                else
                    table_sign_rsn_blp_unadjusted_ec = cat(1,table_sign_rsn_blp_unadjusted_ec,tmp_sign);
                end
            end

            % adjusted
            if isempty(tmp_sign_FDR)
                continue
            else
                if isempty(table_sign_rsn_blp_ec)
                    table_sign_rsn_blp_ec = tmp_sign_FDR;
                else
                    table_sign_rsn_blp_ec = cat(1,table_sign_rsn_blp_ec,tmp_sign_FDR);
                end
            end
        end
    end
end


%% saving structured data and between-group statistics
save('miscellaneous/ws_results_BLP_ec.mat', 'chanlocs', 'rsnlab', 'nch', 'nrsn',...
    'struct_hc_ec', 'struct_sz_ec', 'struct_stat_beta_ec', 'struct_stat_blp_ec',...
    'struct_hc_rsn_ec', 'struct_sz_rsn_ec','struct_stat_rsn_blp_ec',...
    'table_sign_blp_ec', 'table_sign_blp_unadjusted_ec',...
    'table_sign_rsn_blp_ec', 'table_sign_rsn_blp_unadjusted_ec')

%% statistics for global spectra
[~, ~, spec_hc, spec_sz, struct_stat_spec] = func_plot_GA_spectra(struct_hc_ec,struct_sz_ec,0,1,0);

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

% define table of confounding variables
table_conf = table();
table_conf.age = table_panss.age;
table_conf.sex = table_panss.sex;
table_conf.edu = table_panss.edu_years;
table_conf.cpz = table_panss.CPZ;
table_conf.dur = table_panss.duration;
z_conf = table2array(table_conf);

%% demographics analyses

%% sex distributions
n1 = sum(table_dem_hc.sexlab);
N1 = length(table_dem_hc.sexlab);
n2 = sum(table_dem_sz.sexlab);
N2 = length(table_dem_sz.sexlab);

p0 = (n1+n2)/(N1+N2);
n10 = N1*p0;
n20 = N2*p0;

observed = [n1, N1-n1, n2, N2-n2];
expected = [n10, N1-n10, n20, N2-n20];

[~,p,stats] = chi2gof([1,2,3,4],'freq',observed,'expected',expected,'nparams',2);
disp(['sex: p=' num2str(p,'%.4f') ', chi-square=' num2str(stats.chi2stat,'%.4f')])

%% age
v_hc = table_dem_hc.age;
v_sz = table_dem_sz.age;

if lillietest(v_hc) || lillietest(v_sz)
    [p,~,stats] = ranksum(v_hc,v_sz);
    ttype = 'ranksum';
disp(['age: p=' num2str(p,'%.4f') ', ' ttype ', z=' num2str(stats.zval, '%.4f')])
else
    [~,p,~,stats] = ttest2(v_hc,v_sz);
    ttype = 'ttest2';
disp(['age: HC=' num2str(mean(v_hc),'%.2f') '+-' num2str(std(v_hc),'%.2f') ', SZ=' num2str(mean(v_sz),'%.2f') '+-' num2str(std(v_sz),'%.2f') ', p=' num2str(p,'%.4f') ', ' ttype ', t=' num2str(stats.tstat, '%.4f')])
end

%% years in education
v_hc = table_dem_hc.edu_years;
v_sz = table_dem_sz.edu_years;

if lillietest(v_hc) || lillietest(v_sz)
    [p,~,stats] = ranksum(v_hc,v_sz);
    ttype = 'ranksum';
    disp(['years in education: HC=' num2str(median(v_hc),'%.2f') ' [' num2str(min(v_hc)) '; ' num2str(max(v_hc)) '], SZ=' num2str(median(v_sz),'%.2f') ' [' num2str(min(v_sz)) '; ' num2str(max(v_sz)) '], p=' num2str(p,'%.4f') ', ' ttype ', z=' num2str(stats.zval, '%.4f')])
else
    [~,p,~,stats] = ttest2(v_hc,v_sz);
    ttype = 'ttest2';
disp(['years in education: p=' num2str(p,'%.4f') ', ' ttype ', t=' num2str(stats.tstat, '%.4f')])
end

disp(['illness duration: ' num2str(mean(table_dem_sz.duration),'%.2f') '+-' num2str(std(table_dem_sz.duration),'%.2f') ])
disp(['CPZ: ' num2str(mean(table_dem_sz.CPZ),'%.2f') '+-' num2str(std(table_dem_sz.CPZ), '%.2f') ])
disp(['PANSS SUM: ' num2str(mean(table_panss.SUM),'%.2f') '+-' num2str(std(table_panss.SUM), '%.2f') ])
disp(['PANSS GEN: ' num2str(mean(table_panss.GEN),'%.2f') '+-' num2str(std(table_panss.GEN), '%.2f') ])
disp(['PANSS NEG: ' num2str(mean(table_panss.NEG),'%.2f') '+-' num2str(std(table_panss.NEG), '%.2f') ])
disp(['PANSS POS: ' num2str(mean(table_panss.POS),'%.2f') '+-' num2str(std(table_panss.POS), '%.2f') ])

%% compute PANSS correlations - BLP, channel-wise

corr_panss_blp_gen_ec = func_corr_channel_PANSS(table_sign_blp_ec, table_panss, table_conf, 'GEN');
corr_panss_blp_neg_ec = func_corr_channel_PANSS(table_sign_blp_ec, table_panss, table_conf, 'NEG');
corr_panss_blp_pos_ec = func_corr_channel_PANSS(table_sign_blp_ec, table_panss, table_conf, 'POS');
corr_panss_blp_sum_ec = func_corr_channel_PANSS(table_sign_blp_ec, table_panss, table_conf, 'SUM');

%% compute PANSS correlations - BLP, RSN-wise

corr_panss_rsn_blp_gen_ec = func_corr_RSN_PANSS(table_sign_rsn_blp_ec, table_panss, table_conf, 'GEN');
corr_panss_rsn_blp_neg_ec = func_corr_RSN_PANSS(table_sign_rsn_blp_ec, table_panss, table_conf, 'NEG');
corr_panss_rsn_blp_pos_ec = func_corr_RSN_PANSS(table_sign_rsn_blp_ec, table_panss, table_conf, 'POS');
corr_panss_rsn_blp_sum_ec = func_corr_RSN_PANSS(table_sign_rsn_blp_ec, table_panss, table_conf, 'SUM');

%% saving results
save('miscellaneous/ws_correlations_BLP_ec.mat', 'chanlocs', 'rsnlab', 'nch', 'nrsn',...
    'table_panss', 'table_conf', 'z_conf',...
    'corr_panss_blp_gen_ec', 'corr_panss_blp_neg_ec', 'corr_panss_blp_pos_ec', 'corr_panss_blp_sum_ec',...
    'corr_panss_rsn_blp_gen_ec', 'corr_panss_rsn_blp_neg_ec', 'corr_panss_rsn_blp_pos_ec', 'corr_panss_rsn_blp_sum_ec')


