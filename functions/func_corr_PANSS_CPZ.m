function [ output_table ] = func_corr_PANSS_CPZ( input_table, panss, conf_table, flag_conf)

% Function to compute correlations between PANSS scores and CPZ
% equivalent dose, after regressing out the effect of possible confounding
% variables. Analysis is performed on the level of resting-state networks.
% 
% Note that this analysis is exploratory, and therefore we did
% not adjust for multiple comparisons.
%
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 04/02/2025

%% preallocations and constants

% if no prior significant difference, return empty
if isempty(input_table)
    output_table = [];
    return
end

y_cpz = conf_table.cpz; %........CPZ equivalent dose
conf_table_red = removevars(conf_table,'cpz'); %...remove CPZ from confounders
z_conf = table2array(conf_table_red); %.....table of confounding variables

output_struct = struct(...
    'meas', [],...
    'r', [],...
    'p', [],...
    'h', [],...
    'p_FDR', [],...
    'h_FDR', [],...
    'trend', [],...
    'ctype', [],...
    'res_panss', [],...
    'res_cpz', []);

%% compute correlations

    % EEG index
    v_sz = [input_table.(panss)];

    if flag_conf
        % regressing out confounders from PANSS
        lm_eeg = fitlm(z_conf, v_sz);
        res_panss = lm_eeg.Residuals.Raw;

        % regressing out confounders from CPZ
        lm_cpz = fitlm(z_conf, y_cpz);
        res_cpz = lm_cpz.Residuals.Raw;
    else
        res_panss = v_sz;
        res_cpz = y_cpz;
    end

    % test for normality
    if lillietest(res_panss) || lillietest(res_cpz)
        % non-normal: Spearman correlation
        [r, p] = corr(res_panss, res_cpz, 'type', 'Spearman');
        ctype = 'Spearman';
    else
        % normal: Pearson correlation
        [r, p] = corr(res_panss, res_cpz, 'type', 'Pearson');
        ctype = 'Pearson';
    end

    % derive significance
    if p<0.05
        output_struct.h = 1;
    else
        output_struct.h = 0;
    end

    % storing outcomes
    output_struct.meas = panss;
    output_struct.r = r;
    output_struct.p = p;
    output_struct.p_FDR = p;
    output_struct.h_FDR = 0;
    output_struct.trend = sign(r);
    output_struct.ctype = {ctype};
    output_struct.res_panss = res_panss;
    output_struct.res_cpz = res_cpz;

% stucture results
if isempty(output_struct)
    output_table = [];
elseif length(output_struct) == 1
    output_table = table();
    vnames = fieldnames(output_struct);
    for i = 1:length(vnames)
        if ~contains(vnames{i},{'res_panss','res_cpz'})
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