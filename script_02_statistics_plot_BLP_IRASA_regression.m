
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script performs linear mixed effect modeling to see if spectral
% slope estimates are dependent on alpha BLP estimates, with adding subject
% as a random effect. In turn, this script utilizes output generated by
% script_01_analyze_BLP_ec and script_01_analyze_BLP_eo. Analyses are
% performed for both EYES OPEN and EYES CLOSED, for HC and SZ groups
% independently. Only correlations with 20-45 Hz spectral slope are tested,
% where significant between-group differences were found in our previous
% work (Racz et al., 2025).
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

fnames_ec = dir('results_BLP_ec/*.mat');
fnames_eo = dir('results_BLP_eo/*.mat');
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

ns = length(fnames_eo); %...total number of subjects
nch = length(chlab); %......total number of channels
nw = 23; %..................total number of sliding windows
n_sz = length(sz_list); %...number of schizophrenia patients
n_hc = length(hc_list); %...number of healthy controls
FDR_alpha = 0.05; %.........alpha for multiple comparisons adjustment

state_vars = {'EC','EO','dS'}; %..........................analysis types
irasa_vars = {'Beta_lo','Beta_hi','Beta_bb'}; %.......frequency ranges
spec_types = {'mixd','frac','osci'}; %................power spectra types
bands = {'delta','theta','alpha','beta','gamma'}; %...frequency bands
freqs = [1, 4; 4, 8; 8, 13; 13, 25; 25, 45]; %........frequency band limits
nrsn = length(rsnlab); %..............................# of RSNs

%% preallocations

% structures for storing structured analysis outputs (ch and RSN level)
struct_hc_reg = struct();
struct_sz_reg = struct();
struct_hc_rsn_reg = struct();
struct_sz_rsn_reg = struct();

% structures for storing statistical testing outputs
struct_stat_beta_reg = struct();
struct_stat_blp_reg = struct();
struct_stat_rsn_beta_reg = struct();
struct_stat_rsn_blp_reg = struct();

% IRASA variables
for iv = 1:length(irasa_vars)
    struct_hc_reg.(irasa_vars{iv}) = struct();
    struct_sz_reg.(irasa_vars{iv}) = struct();
    for sv = 1:length(state_vars)
        if strcmp(state_vars{sv},'dS')
            % channel-level
            struct_hc_reg.(irasa_vars{iv}).(state_vars{sv}) = zeros(n_hc, nch);
            struct_sz_reg.(irasa_vars{iv}).(state_vars{sv}) = zeros(n_sz, nch);
            struct_stat_beta_reg.(irasa_vars{iv}).(state_vars{sv}) = [];
            % RSN-level
            for nr = 1:nrsn
                struct_hc_rsn_reg.(irasa_vars{iv}).(state_vars{sv}).(Mr.labels{nr}) = zeros(n_hc, 1);
                struct_sz_rsn_reg.(irasa_vars{iv}).(state_vars{sv}).(Mr.labels{nr}) = zeros(n_sz, 1);
                struct_stat_rsn_beta_reg.(irasa_vars{iv}).(state_vars{sv}).(Mr.labels{nr}) = [];
            end
        else
            % channel-level
            struct_hc_reg.(irasa_vars{iv}).(state_vars{sv}) = zeros(n_hc, nw);
            struct_sz_reg.(irasa_vars{iv}).(state_vars{sv}) = zeros(n_sz, nw);
            struct_stat_beta_reg.(irasa_vars{iv}).(state_vars{sv}) = [];
            % RSN-level
            for nr = 1:nrsn
                struct_hc_rsn_reg.(irasa_vars{iv}).(state_vars{sv}).(Mr.labels{nr}) = zeros(n_hc, nw);
                struct_sz_rsn_reg.(irasa_vars{iv}).(state_vars{sv}).(Mr.labels{nr}) = zeros(n_sz, nw);
                struct_stat_rsn_beta_reg.(irasa_vars{iv}).(state_vars{sv}).(Mr.labels{nr}) = [];
            end
        end
    end
end

% BLP variables
for st = 1:length(spec_types)
    struct_hc_reg.(spec_types{st}) = struct();
    struct_sz_reg.(spec_types{st}) = struct();
    for b = 1:length(bands)
        struct_hc_reg.(spec_types{st}).(bands{b}) = struct();
        struct_sz_reg.(spec_types{st}).(bands{b}) = struct();
        for sv = 1:length(state_vars)
            if strcmp(state_vars{sv},'dS') % for EC vs. EO difference, no time resolution
                % channel-level
                struct_hc_reg.(spec_types{st}).(bands{b}).(state_vars{sv}) = zeros(n_hc, nch);
                struct_sz_reg.(spec_types{st}).(bands{b}).(state_vars{sv}) = zeros(n_sz, nch);
                struct_stat_blp_reg.(spec_types{st}).(bands{b}).(state_vars{sv}) = [];
                % RSN-level
                for nr = 1:nrsn
                    struct_hc_rsn_reg.(spec_types{st}).(bands{b}).(state_vars{sv}).(Mr.labels{nr}) = zeros(n_hc, 1);
                    struct_sz_rsn_reg.(spec_types{st}).(bands{b}).(state_vars{sv}).(Mr.labels{nr}) = zeros(n_sz, 1);
                    struct_stat_rsn_blp_reg.(spec_types{st}).(bands{b}).(state_vars{sv}).(Mr.labels{nr}) = [];
                end
            else % for EC and EO, take global average over all channels or RSNs
                % channel-level
                struct_hc_reg.(spec_types{st}).(bands{b}).(state_vars{sv}) = zeros(n_hc, nw);
                struct_sz_reg.(spec_types{st}).(bands{b}).(state_vars{sv}) = zeros(n_sz, nw);
                struct_stat_blp_reg.(spec_types{st}).(bands{b}).(state_vars{sv}) = [];
                % RSN-level
                for nr = 1:nrsn
                    struct_hc_rsn_reg.(spec_types{st}).(bands{b}).(state_vars{sv}).(Mr.labels{nr}) = zeros(n_hc, nw);
                    struct_sz_rsn_reg.(spec_types{st}).(bands{b}).(state_vars{sv}).(Mr.labels{nr}) = zeros(n_sz, nw);
                    struct_stat_rsn_blp_reg.(spec_types{st}).(bands{b}).(state_vars{sv}) = [];
                end
            end
        end
    end
end

%% load and structure data

% HC
for subj = 1:n_hc % iterate over HC subjects
    load(['results_BLP_ec/sch_' hc_list{subj} '_ec_BLP_res.mat'], 'results_seg_ec')
    load(['results_BLP_eo/sch_' hc_list{subj} '_eo_BLP_res.mat'], 'results_seg_eo')

    % IRASA analysis
    for iv = 1:length(irasa_vars) % iterate over fitting ranges
        v_tmp_ec = zeros(nw, nch);
        v_tmp_eo = zeros(nw, nch);
        for t = 1:nw % iterate over sliding windows
            v_tmp_ec(t,:) = results_seg_ec(t).plaw.(irasa_vars{iv});
            v_tmp_eo(t,:) = results_seg_eo(t).plaw.(irasa_vars{iv});
        end

        % aggregate over resting-state networks
        v_tmp_rsn_ec = func_aggregate_RSNs_batch(v_tmp_ec, Mr);
        v_tmp_rsn_eo = func_aggregate_RSNs_batch(v_tmp_eo, Mr);

        % store data
        struct_hc_reg.(irasa_vars{iv}).EC(subj,:) = mean(v_tmp_ec,2)';
        struct_hc_reg.(irasa_vars{iv}).EO(subj,:) = mean(v_tmp_eo,2)';
        struct_hc_reg.(irasa_vars{iv}).dS(subj,:) = mean(v_tmp_ec,1) - mean(v_tmp_eo,1);
        for nr = 1:nrsn
            struct_hc_rsn_reg.(irasa_vars{iv}).EC.(Mr.labels{nr})(subj,:) = v_tmp_rsn_ec(:,nr)';
            struct_hc_rsn_reg.(irasa_vars{iv}).EO.(Mr.labels{nr})(subj,:) = v_tmp_rsn_eo(:,nr)';
            struct_hc_rsn_reg.(irasa_vars{iv}).dS.(Mr.labels{nr})(subj) = mean(v_tmp_rsn_ec(:,nr),1) - mean(v_tmp_rsn_eo(:,nr),1);
        end
    end

    % BLP analysis
    for st = 1:length(spec_types)
        for b = 1:length(bands)
            v_tmp_ec = zeros(nw, nch);
            v_tmp_eo = zeros(nw, nch);
            for t = 1:nw
                freq_ec = results_seg_ec(t).plaw.freq;
                freq_eo = results_seg_eo(t).plaw.freq;
                if strcmp(spec_types{st},'mixd')
                    psd_ec = results_seg_ec(t).plaw.mixd;
                    psd_eo = results_seg_eo(t).plaw.mixd;
                    v_tmp_ec(t,:) = log(sum(psd_ec(freq_ec>=freqs(b,1) & freq_ec<freqs(b,2),:),1));
                    v_tmp_eo(t,:) = log(sum(psd_eo(freq_eo>=freqs(b,1) & freq_eo<freqs(b,2),:),1));
                elseif strcmp(spec_types{st},'frac')
                    psd_ec = results_seg_ec(t).plaw.frac;
                    psd_eo = results_seg_eo(t).plaw.frac;
                    v_tmp_ec(t,:) = log(sum(psd_ec(freq_ec>=freqs(b,1) & freq_ec<freqs(b,2),:),1));
                    v_tmp_eo(t,:) = log(sum(psd_eo(freq_eo>=freqs(b,1) & freq_eo<freqs(b,2),:),1));
                elseif strcmp(spec_types{st},'osci')
                    psdm_ec = results_seg_ec(t).plaw.mixd;
                    psdf_ec = results_seg_ec(t).plaw.frac;
                    psdm_eo = results_seg_eo(t).plaw.mixd;
                    psdf_eo = results_seg_eo(t).plaw.frac;
                    t_tmpm_ec = log(sum(psdm_ec(freq_ec>=freqs(b,1) & freq_ec<freqs(b,2),:),1));
                    t_tmpf_ec = log(sum(psdf_ec(freq_ec>=freqs(b,1) & freq_ec<freqs(b,2),:),1));
                    v_tmp_ec(t,:) = t_tmpm_ec - t_tmpf_ec;
                    t_tmpm_eo = log(sum(psdm_eo(freq_eo>=freqs(b,1) & freq_eo<freqs(b,2),:),1));
                    t_tmpf_eo = log(sum(psdf_eo(freq_eo>=freqs(b,1) & freq_eo<freqs(b,2),:),1));
                    v_tmp_eo(t,:) = t_tmpm_eo - t_tmpf_eo;
                end
            end

            % aggregate over resting-state networks
            v_tmp_rsn_ec = func_aggregate_RSNs_batch(v_tmp_ec, Mr);
            v_tmp_rsn_eo = func_aggregate_RSNs_batch(v_tmp_eo, Mr);

            % storing data
            struct_hc_reg.(spec_types{st}).(bands{b}).EC(subj,:) = mean(v_tmp_ec,2)';
            struct_hc_reg.(spec_types{st}).(bands{b}).EO(subj,:) = mean(v_tmp_eo,2)';
            struct_hc_reg.(spec_types{st}).(bands{b}).dS(subj,:) = mean(v_tmp_ec,1) - mean(v_tmp_eo,1);
            for nr = 1:nrsn
                struct_hc_rsn_reg.(spec_types{st}).(bands{b}).EC.(Mr.labels{nr})(subj,:) = v_tmp_rsn_ec(:,nr)';
                struct_hc_rsn_reg.(spec_types{st}).(bands{b}).EO.(Mr.labels{nr})(subj,:) = v_tmp_rsn_eo(:,nr)';
                struct_hc_rsn_reg.(spec_types{st}).(bands{b}).dS.(Mr.labels{nr})(subj) = mean(v_tmp_rsn_ec(:,nr),1) - mean(v_tmp_rsn_eo(:,nr),1);
            end
        end
    end

    % store data, compute average power spectrum for RSNs (illustration)
    struct_hc_reg.spec(subj).subjID = fnames_eo(subj).name;
    struct_hc_rsn_reg.spec(subj).subjID = fnames_eo(subj).name;

    % compute average spectra (over sliding windows)
    [plaw_ec, plaw_rsn_ec] = func_average_RSN_results_irasa(results_seg_ec, chlab, Mr);
    struct_hc_reg.spec_ec(subj).plaw = plaw_ec;
    struct_hc_rsn_reg.spec_ec(subj).plaw = plaw_rsn_ec;
    [plaw_eo, plaw_rsn_eo] = func_average_RSN_results_irasa(results_seg_eo, chlab, Mr);
    struct_hc_reg.spec_eo(subj).plaw = plaw_eo;
    struct_hc_rsn_reg.spec_eo(subj).plaw = plaw_rsn_eo;
end

% SZ
for subj = 1:n_sz % iterate over HC subjects
    load(['results_BLP_ec/sch_' sz_list{subj} '_ec_BLP_res.mat'], 'results_seg_ec')
    load(['results_BLP_eo/sch_' sz_list{subj} '_eo_BLP_res.mat'], 'results_seg_eo')

    % IRASA analysis
    for iv = 1:length(irasa_vars) % iterate over fitting ranges
        v_tmp_ec = zeros(nw, nch);
        v_tmp_eo = zeros(nw, nch);
        for t = 1:nw % iterate over sliding windows
            v_tmp_ec(t,:) = results_seg_ec(t).plaw.(irasa_vars{iv});
            v_tmp_eo(t,:) = results_seg_eo(t).plaw.(irasa_vars{iv});
        end

        % aggregate over resting-state networks
        v_tmp_rsn_ec = func_aggregate_RSNs_batch(v_tmp_ec, Mr);
        v_tmp_rsn_eo = func_aggregate_RSNs_batch(v_tmp_eo, Mr);

        % store data
        struct_sz_reg.(irasa_vars{iv}).EC(subj,:) = mean(v_tmp_ec,2)';
        struct_sz_reg.(irasa_vars{iv}).EO(subj,:) = mean(v_tmp_eo,2)';
        struct_sz_reg.(irasa_vars{iv}).dS(subj,:) = mean(v_tmp_ec,1) - mean(v_tmp_eo,1);
        for nr = 1:nrsn
            struct_sz_rsn_reg.(irasa_vars{iv}).EC.(Mr.labels{nr})(subj,:) = v_tmp_rsn_ec(:,nr)';
            struct_sz_rsn_reg.(irasa_vars{iv}).EO.(Mr.labels{nr})(subj,:) = v_tmp_rsn_eo(:,nr)';
            struct_sz_rsn_reg.(irasa_vars{iv}).dS.(Mr.labels{nr})(subj) = mean(v_tmp_rsn_ec(:,nr),1) - mean(v_tmp_rsn_eo(:,nr),1);
        end
    end

    % BLP analysis
    for st = 1:length(spec_types)
        for b = 1:length(bands)
            v_tmp_ec = zeros(nw, nch);
            v_tmp_eo = zeros(nw, nch);
            for t = 1:nw
                freq_ec = results_seg_ec(t).plaw.freq;
                freq_eo = results_seg_eo(t).plaw.freq;
                if strcmp(spec_types{st},'mixd')
                    psd_ec = results_seg_ec(t).plaw.mixd;
                    psd_eo = results_seg_eo(t).plaw.mixd;
                    v_tmp_ec(t,:) = log(sum(psd_ec(freq_ec>=freqs(b,1) & freq_ec<freqs(b,2),:),1));
                    v_tmp_eo(t,:) = log(sum(psd_eo(freq_eo>=freqs(b,1) & freq_eo<freqs(b,2),:),1));
                elseif strcmp(spec_types{st},'frac')
                    psd_ec = results_seg_ec(t).plaw.frac;
                    psd_eo = results_seg_eo(t).plaw.frac;
                    v_tmp_ec(t,:) = log(sum(psd_ec(freq_ec>=freqs(b,1) & freq_ec<freqs(b,2),:),1));
                    v_tmp_eo(t,:) = log(sum(psd_eo(freq_eo>=freqs(b,1) & freq_eo<freqs(b,2),:),1));
                elseif strcmp(spec_types{st},'osci')
                    psdm_ec = results_seg_ec(t).plaw.mixd;
                    psdf_ec = results_seg_ec(t).plaw.frac;
                    psdm_eo = results_seg_eo(t).plaw.mixd;
                    psdf_eo = results_seg_eo(t).plaw.frac;
                    t_tmpm_ec = log(sum(psdm_ec(freq_ec>=freqs(b,1) & freq_ec<freqs(b,2),:),1));
                    t_tmpf_ec = log(sum(psdf_ec(freq_ec>=freqs(b,1) & freq_ec<freqs(b,2),:),1));
                    v_tmp_ec(t,:) = t_tmpm_ec - t_tmpf_ec;
                    t_tmpm_eo = log(sum(psdm_eo(freq_eo>=freqs(b,1) & freq_eo<freqs(b,2),:),1));
                    t_tmpf_eo = log(sum(psdf_eo(freq_eo>=freqs(b,1) & freq_eo<freqs(b,2),:),1));
                    v_tmp_eo(t,:) = t_tmpm_eo - t_tmpf_eo;
                end
            end

            % aggregate over resting-state networks
            v_tmp_rsn_ec = func_aggregate_RSNs_batch(v_tmp_ec, Mr);
            v_tmp_rsn_eo = func_aggregate_RSNs_batch(v_tmp_eo, Mr);

            % storing data
            struct_sz_reg.(spec_types{st}).(bands{b}).EC(subj,:) = mean(v_tmp_ec,2)';
            struct_sz_reg.(spec_types{st}).(bands{b}).EO(subj,:) = mean(v_tmp_eo,2)';
            struct_sz_reg.(spec_types{st}).(bands{b}).dS(subj,:) = mean(v_tmp_ec,1) - mean(v_tmp_eo,1);
            for nr = 1:nrsn
                struct_sz_rsn_reg.(spec_types{st}).(bands{b}).EC.(Mr.labels{nr})(subj,:) = v_tmp_rsn_ec(:,nr)';
                struct_sz_rsn_reg.(spec_types{st}).(bands{b}).EO.(Mr.labels{nr})(subj,:) = v_tmp_rsn_eo(:,nr)';
                struct_sz_rsn_reg.(spec_types{st}).(bands{b}).dS.(Mr.labels{nr})(subj) = mean(v_tmp_rsn_ec(:,nr),1) - mean(v_tmp_rsn_eo(:,nr),1);
            end
        end
    end

    % store data, compute average power spectrum for RSNs (illustration)
    struct_sz_reg.spec(subj).subjID = fnames_eo(subj).name;
    struct_sz_rsn_reg.spec(subj).subjID = fnames_eo(subj).name;

    % compute average spectra (over sliding windows)
    [plaw_ec, plaw_rsn_ec] = func_average_RSN_results_irasa(results_seg_ec, chlab, Mr);
    struct_sz_reg.spec_ec(subj).plaw = plaw_ec;
    struct_sz_rsn_reg.spec_ec(subj).plaw = plaw_rsn_ec;
    [plaw_eo, plaw_rsn_eo] = func_average_RSN_results_irasa(results_seg_eo, chlab, Mr);
    struct_sz_reg.spec_eo(subj).plaw = plaw_eo;
    struct_sz_rsn_reg.spec_eo(subj).plaw = plaw_rsn_eo;
end

% storing channel and RSN labels for convenience
struct_hc_reg.chlab = chlab;
struct_sz_reg.chlab = chlab;
struct_hc_rsn_reg.rsnlab = rsnlab;
struct_sz_rsn_reg.rsnlab = rsnlab;


%% Linear mixed-effect modeling of spectral slope on alpha power

[lme_osci_ec_hc, lme_osci_table_ec_hc] = func_LME_slope_on_BLP(struct_hc_reg,'Beta_hi','osci','alpha','EC');
[lme_osci_ec_sz, lme_osci_table_ec_sz] = func_LME_slope_on_BLP(struct_sz_reg,'Beta_hi','osci','alpha','EC');

[lme_osci_eo_hc, lme_osci_table_eo_hc] = func_LME_slope_on_BLP(struct_hc_reg,'Beta_hi','osci','alpha','EO');
[lme_osci_eo_sz, lme_osci_table_eo_sz] = func_LME_slope_on_BLP(struct_sz_reg,'Beta_hi','osci','alpha','EO');

