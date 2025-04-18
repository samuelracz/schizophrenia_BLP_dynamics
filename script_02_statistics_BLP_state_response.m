
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script provides complimentary analysis in testing the state response
% (transitioning from EO to EC) in both HC and SZ groups independently.
% The following statistical evaluations are carried out:
% - difference in mixed, fractal and oscillatory BLP mean over time
% - difference in mixed, fractal and oscillatory BLP variance over time
% - illustration of global EYES-OPEN vs. EYES CLOSED EEG power spectra
% - demographics comparison between HC and SZ
% - correlation analysis between findings and PANSS scores in SZ
%
%
% note that that is stored such that the first 30 subjects (sch_002 - 
% sch_050) are SZ patients, the remaining 31 subjects (sch_101 - sch_333)
% are healthy controls.
%
% NOTE: this script utilizes aggregate output generated previously via
% script_02_statistics_BLP_ec.m and script_02_statistics_BLP_eo.m in
% loading ws_results_BLP_ec.mat and ws_results_BLP_eo.mat
%
% This script evaluates EYES OPEN vs. EYES CLOSED difference measures
% 
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 04/18/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('functions'))
addpath('miscellaneous')

load('ws_results_BLP_ec.mat','struct_hc_ec','struct_sz_ec')
load('ws_results_BLP_eo.mat','struct_hc_eo','struct_sz_eo')
load('ws_55ch_montage.mat')
load('ws_55ch_chanlocs.mat')

subj_list = {'002','003','004','005','006','007','008','009','013','014','016','018','021','022','023','026','027','029','030','031','033','034','035','036','039','040','044','045','046','050',...
    '101','102','103','104','105','106','107','108','109','111','113','115','116','117','118','120','121','122','123','124','125','126','127','128','129','131','136','137','138','140','333'};

hc_list = {'101','102','103','104','105','106','107','108','109','111','113','115','116','117','118','120','121','122','123','124','125','126','127','128','129','131','136','137','138','140','333'};
sz_list = {'002','003','004','005','006','007','008','009','013','014','016','018','021','022','023','026','027','029','030','031','033','034','035','036','039','040','044','045','046','050'};

chlab = M.lab; %....................................channel labels

ns = length(subj_list); %...total number of subjects
nch = length(chlab); %......total number of channels
nw = 23; %..................total number of sliding windows
n_sz = length(sz_list); %...number of schizophrenia patients
n_hc = length(hc_list); %...number of healthy controls
FDR_alpha = 0.05; %.........alpha for multiple comparisons adjustment

meanvars = {'mean','var'}; %..........................analysis types
irasa_vars = {'Beta_lo','Beta_hi'}; %.......frequency ranges
spec_types = {'mixd','frac','osci'}; %................power spectra types
bands = {'delta','theta','alpha','beta','gamma'}; %...frequency bands
freqs = [1, 4; 4, 8; 8, 13; 13, 25; 25, 45]; %........frequency band limits

%% preallocations

% structures for storing structured analysis outputs (ch level)
struct_hc_stat_res = struct();
struct_sz_stat_res = struct();

% IRASA analysis
for iv = 1:length(irasa_vars)
    struct_hc_stat_res.(irasa_vars{iv}) = struct();
    struct_sz_stat_res.(irasa_vars{iv}) = struct();
    for mv = 1:length(meanvars)
        % channel-level
        struct_hc_stat_res.(irasa_vars{iv}).(meanvars{mv}) = [];
        struct_sz_stat_res.(irasa_vars{iv}).(meanvars{mv}) = [];
    end
end

% BLP analysis
for st = 1:length(spec_types)
    struct_hc_stat_res.(spec_types{st}) = struct();
    struct_sz_stat_res.(spec_types{st}) = struct();
    for b = 1:length(bands)
        struct_hc_stat_res.(spec_types{st}).(bands{b}) = struct();
        struct_sz_stat_res.(spec_types{st}).(bands{b}) = struct();
        for mv = 1:length(meanvars)
            % channel-level
            struct_hc_stat_res.(spec_types{st}).(bands{b}).(meanvars{mv}) = [];
            struct_sz_stat_res.(spec_types{st}).(bands{b}).(meanvars{mv}) = [];
        end
    end
end

%% spectral slope analyses
for iv = 1:length(irasa_vars)
    ivn = irasa_vars{iv};
    for mv = 1:length(meanvars)
        mvn = meanvars{mv};
        % Healthy controls
        struct_hc_stat_res.(ivn).(mvn) =...
            func_pairwise_FDR(func_dS_slope(struct_hc_ec,struct_hc_eo,chlab,ivn,mvn),0.05);
        struct_hc_stat_res.(ivn).([mvn '_global']) =...
            func_pairwise_FDR(func_dS_slope(struct_hc_ec,struct_hc_eo,[],ivn,mvn),0.05);
        % Schizophrenia patients
        struct_sz_stat_res.(ivn).(mvn) =...
            func_pairwise_FDR(func_dS_slope(struct_sz_ec,struct_sz_eo,chlab,ivn,mvn),0.05);
        struct_sz_stat_res.(ivn).([mvn '_global']) =...
            func_pairwise_FDR(func_dS_slope(struct_sz_ec,struct_sz_eo,[],ivn,mvn),0.05);
    end
end

%% BLP analyses
for st = 1:length(spec_types)
    stn = spec_types{st};
    for b = 1:length(bands)
        bn = bands{b};
        for mv = 1:length(meanvars)
            mvn = meanvars{mv};
            % Healthy controls
            struct_hc_stat_res.(stn).(bn).(mvn) = ...
                func_pairwise_FDR(func_dS_BLP(struct_hc_ec,struct_hc_eo,chlab,stn,bn,mvn),0.05);
            struct_hc_stat_res.(stn).(bn).([mvn '_global']) = ...
                func_pairwise_FDR(func_dS_BLP(struct_hc_ec,struct_hc_eo,[],stn,bn,mvn),0.05);
            % Schizophrenia patients
            struct_sz_stat_res.(stn).(bn).(mvn) = ...
                func_pairwise_FDR(func_dS_BLP(struct_sz_ec,struct_sz_eo,chlab,stn,bn,mvn),0.05);
            struct_sz_stat_res.(stn).(bn).([mvn '_global']) = ...
                func_pairwise_FDR(func_dS_BLP(struct_sz_ec,struct_sz_eo,[],stn,bn,mvn),0.05);
        end
    end
end 

%% saving data
save('miscellaneous/ws_results_BLP_response.mat', 'struct_hc_stat_res', 'struct_sz_stat_res')