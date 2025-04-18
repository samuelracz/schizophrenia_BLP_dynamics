function [ output_table ] = func_corr_RSN_CPZ( input_table, conf_table)

% Function to compute correlations between spectral indices and CPZ
% equivalent dose, after regressing out the effect of possible confounding
% variables. Analysis is performed on the level of resting-state networks.
% 
% Note that this analysis is exploratory, and therefore we did
% not adjust for multiple comparisons.
%
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 02/20/2025

%% preallocations and constants

% if no prior significant difference, return empty
if isempty(input_table)
    output_table = [];
    return
end

ntest = size(input_table,1); %..........# of comparisons
y_cpz = conf_table.cpz; %........CPZ equivalent dose
conf_table_red = removevars(conf_table,'cpz'); %...remove CPZ from confounders
z_conf = table2array(conf_table_red); %.....table of confounding variables

output_struct = struct(...
    'meas', [],...
    'rsn', [],...
    'rsnID', [],...
    'r', [],...
    'p', [],...
    'h', [],...
    'p_FDR', [],...
    'h_FDR', [],...
    'trend', [],...
    'ctype', [],...
    'res_eeg', [],...
    'res_cpz', []);

%% compute correlations
for n = 1:ntest

    % EEG index
    v_sz = input_table.v_sz{n};

    % regressing out confounders from EEG index
    lm_eeg = fitlm(z_conf, v_sz);
    res_eeg = lm_eeg.Residuals.Raw;

    % regressing out confounders from PANSS
    lm_cpz = fitlm(z_conf, y_cpz);
    res_cpz = lm_cpz.Residuals.Raw;

    % test for normality
    if lillietest(res_eeg) || lillietest(res_cpz)
        % non-normal: Spearman correlation
        [r, p] = corr(res_eeg, res_cpz, 'type', 'Spearman');
        ctype = 'Spearman';
    else
        % normal: Pearson correlation
        [r, p] = corr(res_eeg, res_cpz, 'type', 'Pearson');
        ctype = 'Pearson';
        output_struct(n).ctype = 'Pearson';
    end

    % derive significance
    if p<0.05
        output_struct(n).h = 1;
    else
        output_struct(n).h = 0;
    end

    % storing outcomes
    output_struct(n).meas = input_table.meas{n};
    output_struct(n).rsn = input_table.rsn{n};
    output_struct(n).rsnID = input_table.rsnID(n);
    output_struct(n).r = r;
    output_struct(n).p = p;
    output_struct(n).p_FDR = p;
    output_struct(n).h_FDR = 0;
    output_struct(n).trend = sign(r);
    output_struct(n).ctype = ctype;
    output_struct(n).res_eeg = res_eeg;
    output_struct(n).res_cpz = res_cpz;

end

% stucture results
if isempty(output_struct)
    output_table = [];
elseif length(output_struct) == 1
    output_table = table();
    vnames = fieldnames(output_struct);
    for i = 1:length(vnames)
        if ~contains(vnames{i},{'res_eeg','res_panss'})
            output_table.(vnames{i}) = output_struct.(vnames{i});
        else
            output_table.(vnames{i}) = {output_struct.(vnames{i})};
        end
    end
else
    output_table = struct2table(output_struct);
    output_table = sortrows(output_table, 'p', 'ascend');
end

end