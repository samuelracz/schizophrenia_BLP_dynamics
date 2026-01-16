addpath(genpath('functions'))
addpath('miscellaneous')

OD = cd;
fnames_ec = dir('data_ec_mara_avg_256Hz/*.edf');

%% analysis constants
nch = 55; %...................number of channels
fs = 256; %...................frequency of sampling
ns = length(fnames_ec); %.....total number of participants
freq_bins = [(1:44)',(2:45)'];
Nc = 20; %....................number of phase cycles per frequency
Nf = size(freq_bins,1); %.....number of frequency estimates

%% analysis
for subj = 1:ns
    tic
    % load data
    % load(['data_ec_mara_avg_256Hz_mat/' fnames_ec(subj).name])
    [hdr, data_ec] = edfread([fnames_ec(subj).folder '/' fnames_ec(subj).name]);
    eeg_ec = data_ec(1:nch,1:30*fs)'; % remove trigger channel

    % z-score data for BLP analysis
    eeg_ec = zscore(eeg_ec);
    
    % preallocate
    res_mean = zeros(Nf,nch);
    res_var = zeros(Nf,nch);

    for f = 1:Nf
        f_lo = freq_bins(f,1);
        f_hi = freq_bins(f,2);
        W = floor((1/f_lo)*Nc*fs); % window size
        wstep = floor((1/8)*W); % window overlap
        nseg = floor((30*fs-W)/wstep)+1; %......total number of windows


        % reshape data
        eeg_seg_ec = zeros(W,nch,nseg);
        for seg = 1:nseg
            eeg_seg_ec(:,:,seg) = eeg_ec((seg-1)*wstep+1:(seg-1)*wstep+W,:);
        end
        
        PSD = zeros(nseg,nch);

        for seg = 1:nseg
            [A,fseg] = pwelch(eeg_seg_ec(:,:,seg),[],[],f_lo:0.1:f_hi,fs);
            PSD(seg,:) = log(mean(A,1));
        end

        res_mean(f,:) = mean(PSD,1);
        res_var(f,:) = std(PSD,[],1);

    end

    % saving results
    if ~exist('results_phasecycle_ec', 'dir')
       mkdir('results_phasecycle_ec')
    end
    subj_saved = subj;
    hdr_saved = hdr;
    fname_ec = [fnames_ec(subj).name(1:end-8) '_phasecycle_res.mat'];
    save(['results_phasecycle_ec/' fname_ec], 'res_mean', 'res_var', 'freq_bins', 'hdr_saved', 'subj_saved');
    disp([fname_ec, ' || time elapsed: ', num2str(toc)])
end