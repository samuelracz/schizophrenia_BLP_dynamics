
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis script to compute spectral indices in a dynamic manner using
% Irregular Resampling Auto-Spectral Analysis (IRASA) by Wen & Liu (2016).
% Analysis is performed on the same datasets evaluated in Racz et al.
% (2025); only for eyes-open data
%
% NOTE: this script loads the pre-processed data from MATLAB workspaces,
% however corresponding .edf files are also provided.
%
% Reference for the IRASA method:
% Wen H. and Liu Z. (2015) Separating Fractal and Oscillatory Components in
% the Power Spectrum of Neurophysiological Signals
% 
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 04/17/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('functions'))

OD = cd;
fnames_eo = dir('data_eo_mara_avg_256Hz_mat/*.mat');

%% analysis constants
nch = 55; %...................number of channels
fs = 256; %...................frequency of sampling
ns = length(fnames_eo); %.....total number of participants

W = 8; %......................sliding window analysis window size (s)
wstep = 1; %..................sliding window analysis step size (s)
nseg = (30-W)/wstep+1; %......total number of windows

frange_bb = [1 45]; %.........broadband frequency range
frange_lo = [1, 4]; %.........low-frequency range
frange_hi = [20 45];%.........high-frequency range
hset = 1.1:0.1:2.6; %.........set of resampling factors

%% analysis
for subj = 1
    tic
    % load data
    load(['data_eo_mara_avg_256Hz_mat/' fnames_eo(subj).name])
    % [hdr, data_eo] = edfread([fnames_eo(subj).folder '/' fnames_eo(subj).name]);
    eeg_eo = data_eo(1:nch,1:30*fs)'; % remove trigger channel
    
    % z-score data for BLP analysis
    eeg_eo = zscore(eeg_eo);

    chlab = cell(1,nch);
    for ch = 1:nch
        chlab{ch} = hdr.label{ch};
    end

    % IRASA from full signal (30s) - only for sanity check
    spec_eo = amri_sig_fractal(eeg_eo,fs,'hset',hset);
    plaw_eo = func_IRASA_plawfit_multi(spec_eo, frange_bb, frange_lo, frange_hi);
    
    % reshape data
    eeg_seg_eo = zeros(W*fs,nch,nseg);
    for seg = 1:nseg
        eeg_seg_eo(:,:,seg) = eeg_eo((seg-1)*wstep*fs+1:(seg-1)*wstep*fs+W*fs,:);
    end
    
    % IRASA from segmented data (8s segments, 75% overlap)
    results_seg_eo = struct('ID',[],'spec',[],'plaw',[],'flag_af',[]);
    for seg = 1:nseg
        spec_seg = amri_sig_fractal(squeeze(eeg_seg_eo(:,:,seg)),fs,'hset',hset);
        plaw_seg = func_IRASA_plawfit_multi(spec_seg,frange_bb,frange_lo,frange_hi);

        % store output
        results_seg_eo(seg).ID = seg;
        results_seg_eo(seg).spec = spec_seg;
        results_seg_eo(seg).plaw = plaw_seg;
        results_seg_eo(seg).flag_af = 0;
    end

    % saving results
    if ~exist('results_BLP_eo', 'dir')
       mkdir('results_BLP_eo')
    end
    subj_saved = subj;
    hdr_saved = hdr;
    fname_eo = [fnames_eo(subj).name(1:end-8) '_BLP_res.mat'];
    save(['results_BLP_eo/' fname_eo], 'spec_eo', 'plaw_eo', 'results_seg_eo', 'hdr_saved', 'subj_saved');
    disp([fname_eo, ' || time elapsed: ', num2str(toc)])
end