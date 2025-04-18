
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis script to perform surrogate data testing for alpha BLP dynamics.
% Only for eyes-closed data
%
% NOTE: this script loads the pre-processed data from MATLAB workspaces,
% however corresponding .edf files are also provided.
%
% Reference for the surrogate data testing:
% Theiler et al. "Testing for nonlinearity in time series: the method of 
% surrogate data." Physica D: Nonlinear Phenomena 58.1-4 (1992): 77-94.
% 
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 04/17/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('functions'))

OD = cd;
fnames_ec = dir('data_ec_mara_avg_256Hz_mat/*.mat');

%% analysis constants
nch = 55; %...................number of channels
fs = 256; %...................frequency of sampling
ns = length(fnames_ec); %.....total number of participants

W = 8; %......................sliding window analysis window size (s)
wstep = 1; %..................sliding window analysis step size (s)
nseg = (30-W)/wstep+1; %......total number of windows
nsur = 100;

%% analysis
for subj = 1:ns
    tic
    % load data
    load(['data_ec_mara_avg_256Hz_mat/' fnames_ec(subj).name])
    % [hdr, data_ec] = edfread([fnames_ec(subj).folder '/' fnames_ec(subj).name]);
    eeg_ec = data_ec(1:nch,1:30*fs)'; % remove trigger channel

    % z-score data for BLP analysis
    eeg_ec = zscore(eeg_ec);
    
    chlab = cell(1,nch);
    for ch = 1:nch
        chlab{ch} = hdr.label{ch};
    end

    % pre-allocate 
    results_seg_ec = struct('ID',[],'mean',[],'var',[],'mean_sur',[],'var_sur',[]);

    %% analyze real data
    % reshape data
    eeg_seg_ec = zeros(W*fs,nch,nseg);
    for seg = 1:nseg
        eeg_seg_ec(:,:,seg) = eeg_ec((seg-1)*wstep*fs+1:(seg-1)*wstep*fs+W*fs,:);
    end
    
    % get alpha BLP from windows (8s segments, 75% overlap)
    seg_alpha = zeros(nseg,nch);
    for seg = 1:nseg
        [A,f] = pwelch(squeeze(eeg_seg_ec(:,:,seg)), [], [], 1:0.1:45, fs);
        seg_alpha(seg,:) = log(sum(A(f>=8 & f<=13,:),1));
    end

    % store output
    results_seg_ec.ID = num2str(subj,'%.3d');
    results_seg_ec.mean = mean(seg_alpha,1);
    results_seg_ec.var = std(seg_alpha,[],1);

    %% surrogate data testing
    sur_mean = zeros(nsur,nch);
    sur_var = zeros(nsur,nch);
    for sur = 1:nsur
        % generate surrogate data
        eeg_ec_sur = func_surrogate_EEG(eeg_ec);

        % reshape data
        eeg_seg_ec_sur = zeros(W*fs,nch,nseg);
        for seg = 1:nseg
            eeg_seg_ec_sur(:,:,seg) = eeg_ec_sur((seg-1)*wstep*fs+1:(seg-1)*wstep*fs+W*fs,:);
        end

        % get alpha BLP from windows (8s segments, 87.5% overlap)
        seg_alpha = zeros(nseg,nch);
        for seg = 1:nseg
            [A,f] = pwelch(squeeze(eeg_seg_ec_sur(:,:,seg)), [], [], 1:0.1:45, fs);
            seg_alpha(seg,:) = log(sum(A(f>=8 & f<=13,:),1));
        end

        % store output
        sur_mean(sur,:) = mean(seg_alpha,1);
        sur_var(sur,:) = std(seg_alpha,[],1);
    end

    % store output
    results_seg_ec.mean_sur = sur_mean;
    results_seg_ec.var_sur = sur_var;

    % saving results
    if ~exist('results_surrogate_ec', 'dir')
       mkdir('results_surrogate_ec')
    end
    subj_saved = subj;
    hdr_saved = hdr;
    fname_ec = [fnames_ec(subj).name(1:end-8) '_surrogate_res.mat'];
    save(['results_surrogate_ec/' fname_ec], 'results_seg_ec', 'hdr_saved', 'subj_saved');
    disp([fname_ec, ' || time elapsed: ', num2str(toc)])
end