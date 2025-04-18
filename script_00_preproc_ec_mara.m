
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a template/example script %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make sure to add the recent version of EEGLAB (Delorme & Makeig, 2004) to
% path. Also make sure that the loadcurry plugin is installed.
%
% This script is only illustrating the data pre-processing pipeline, but is
% not intended to be used as is. The code is derived from the history
% script generated by EEGLAB.
%
% The beginning and end positions (in seconds) of the 30-second
% epochs analyzed in Racz et al. (2025) is stored in the Matlab workspace
% 'fnames_times_ec.mat'. In the article, we used a zero-phase Butterworth
% filter of order 4 instead of the standard FIR filter implemented in
% EEGLAB. 
%
% NOTE: in some versions, the processMARA() might crash with options set to
% [0,0,0,0,1] (no visualization, automatic rejection)
% 
% Author: F. Samuel Racz, 2024, The University of Texas at Austin
% email: fsr324@austin.utexas.edu

addpath('eeglab2024.0')
load('ws_epoch_locations_EC_EO.mat')
subj_list = {'002','003','004','005','006','007','008','009','013','014','016','018','021','022','023','026','027','029','030','031','033','034','035','036','039','040','044','045','046','050',...
    '101','102','103','104','105','106','107','108','109','111','113','115','116','117','118','120','121','122','123','124','125','126','127','128','129','131','136','137','138','140','333'};
ns = length(subj_list);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define path to the subset of 61 participants included in final sample %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path_to_raw_eeg = 'sch_resting_64ch_ec';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define output path to save pre-processed EEG and MARA statistics %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path_to_pre_processed_eeg = 'data_ec_mara_avg_256Hz';
path_to_mara_statistics = 'dataset_ec_MARA';

for subj = 1:ns
    %% load data
    subjID = subj_list{subj};
    fname_in = ['/sch_' subjID '_ec.dat'];
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    % eeglab;
    EEG = loadcurry([path_to_raw_eeg, '/', fname_in]);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off');
    EEG.etc.eeglabvers = '2021.0'; % this tracks which version of EEGLAB is being used, you may ignore it
    EEG = eeg_checkset(EEG);
    EEG=pop_chanedit(EEG, 'lookup',[cd,'\\eeglab2024.0\\plugins\\dipfit\\standard_BESA\\standard-10-5-cap385.elp']);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
    EEG.setname = '0. raw data';
    
    %% select channels
    EEG = pop_select( EEG,'channel',{'FP1' 'FPZ' 'FP2' 'AF3' 'AF4' 'F7' 'F5' 'F3' 'FZ' 'F4' 'F6' 'F8' 'FT7' 'FC5' 'FC3' 'FC1' 'FC2' 'FC4' 'FC6' 'FT8' 'T7' 'C5' 'C3' 'C1' 'CZ' 'C2' 'C4' 'C6' 'T8' 'CP5' 'CP3' 'CP1' 'CPZ' 'CP2' 'CP4' 'CP6' 'P7' 'P5' 'P3' 'P1' 'PZ' 'P2' 'P4' 'P6' 'P8' 'PO7' 'PO5' 'PO3' 'POZ' 'PO4' 'PO6' 'PO8' 'O1' 'OZ' 'O2'});
    EEG = eeg_checkset( EEG );
    EEG.setname = '1. raw data (55 ch)';
     
    %% segment data (select clean epoch)
    EEG = pop_select( EEG, 'time',[fnames_ec(subj).start fnames_ec(subj).stop-1/1000] );
    EEG = eeg_checkset( EEG );
    EEG.setname = '2. 30s data';
       
    %% high-pass filtering
    EEG = pop_eegfiltnew(EEG, 'locutoff',0.1);
    EEG = eeg_checkset( EEG );
    EEG.setname = '3. after HP';

    %% low-pass filtering (plus notch filters)
    EEG = pop_eegfiltnew(EEG, 'locutoff',49.5,'hicutoff',50.5,'revfilt',1);
    EEG = eeg_checkset( EEG );
    EEG = pop_eegfiltnew(EEG, 'locutoff',99.5,'hicutoff',100.5,'revfilt',1);
    EEG = eeg_checkset( EEG );
    EEG = pop_eegfiltnew(EEG,'hicutoff',128);
    EEG = eeg_checkset( EEG );
    EEG.setname = '4. filtered data';
    
    %% resample data
    EEG = pop_resample( EEG, 256);
    EEG = eeg_checkset( EEG );
    EEG.setname = '5. downsampled data';

    %% ICA
    EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on');
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
    EEG.setname = '6. ICA';
    
    %% artefact removal with MARA
    options = [0,0,0,0,1];
    [ALLEEG,EEG,CURRENTSET] = processMARA(ALLEEG,EEG,CURRENTSET,options);
    EEG = eeg_checkset( EEG );
    EEG.setname = '7. MARA';
    
    %% re-referencing to common average
    EEG = pop_reref( EEG, []);
    EEG = eeg_checkset( EEG );
    EEG.setname = '8. re-reference';
    
    %% saving data
    fname_out = [fname_in(1:end-4), '_30s.edf'];
    pop_writeeeg(EEG, [path_to_pre_processed_eeg, '\', fname_out], 'TYPE','EDF');

    %% saving MARA statistics
    fname_out = [fname_in(1:end-4), '_MARA.set'];
    EEG = pop_saveset( EEG, 'filename' ,fname_out, 'filepath', path_to_mara_statistics);

    EEG = eeg_checkset( EEG );
    save([path_to_mara_statistics '/' fname_in(1:end-4) '_MARA.mat'])
end