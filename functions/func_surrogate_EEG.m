function [eeg_sur] = func_surrogate_EEG(eeg)

% Function to generate surrogate time series from EEG data (timestamp*nch)
% using the phase randomization technique introduced by Theiler and
% colleagues. Nonlinear structure is eliminated via randomizing the phases,
% while linear properties are preserved.
%
% Reference for the surrogate data testing:
% Theiler et al. "Testing for nonlinearity in time series: the method of 
% surrogate data." Physica D: Nonlinear Phenomena 58.1-4 (1992): 77-94.
%
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 04/17/2025

% ensure Nsample*Nch structure
if size(eeg,1) < size(eeg,2)
    eeg = eeg';
end

[N, nch] = size(eeg);

% ensure number of samples is odd
if mod(N,2) == 0
    N = N-1;
    eeg = eeg(1:N,:);
end

% get helper variables
N2 = (N-1)/2;
seg1 = 2:N2+1;
seg2 = N2+2:N;

% Fourier transform
A = fft(eeg);
A_sur = A;

% randomizing phases
phase_rand = rand([N2, nch]);
phase_rand1 = exp(2*pi*1i*phase_rand);
phase_rand2 = conj(flipud(phase_rand1));

A_sur(seg1,:) = A(seg1,:).*phase_rand1;
A_sur(seg2,:) = A(seg2,:).*phase_rand2;

% reconstruct data
eeg_sur = [eeg(1,:); real(ifft(A_sur))];

end