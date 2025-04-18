function [ lme, data_table ] = func_LME_slope_on_BLP(struct_data, vbeta, stype, fband, vstate, varargin)

% This function carries out a linear mixed effect modeling such as:
% predictor variable: BLP (spectrum type + frequency band)
% response variable: spectral slope (low- or high-range)
% grouping variable (random effect): study participant ID
% 
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 04/18/2025


% varargin for RSN-level analysis

% get variables
if isempty(varargin) % channel-level analysis (global average)
    var_beta = struct_data.(vbeta).(vstate);
    var_blp = struct_data.(stype).(fband).(vstate);
elseif length(varargin) == 1 % RSN-level analysis: provide RSN label
    var_beta = struct_data.(vbeta).(vstate).(varargin{1});
    var_blp = struct_data.(stype).(fband).(vstate).(varargin{1});
else
    error('too many input arguments')
end

ns = size(var_beta,1);
nw = size(var_beta,2);

% structure data
data_table = table(zeros(ns*nw,1), zeros(ns*nw,1), cell(ns*nw,1), 'VariableNames', {'beta','BLP','subjID'});
data_table.beta = reshape(var_beta,[ns*nw,1]);
data_table.BLP = reshape(var_blp,[ns*nw,1]);

% fill subject ID column
for subj = 1:ns
    subjID = ['subj' num2str(subj,'%.3d')];
    for w = 1:nw
        ind = (subj-1)*nw + w;
        data_table.subjID{ind} = subjID;
    end
end

data_table.subjID = categorical(data_table.subjID);

lme = fitlme(data_table, 'beta ~ BLP + (1|subjID) + (-1 + BLP|subjID)');

end