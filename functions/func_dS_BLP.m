function [ output_table ] = func_dS_BLP(data_ec, data_eo, chlab, stype, frange, meas)

% Function compare spectral indices obtained from EYES OPEN and EYES CLOSED
% states (paired comparison design)
%
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 02/25/2025

%% preallocation and analysis constants

if isempty(chlab) % global analysis
    nch = 1;
    chlab = {'global'};
else % all channels individually
    nch = length(chlab); % # of channels
end

output_table = table(cell(nch,1), cell(nch,1), zeros(nch,1),...
    zeros(nch,1), zeros(nch,1),...
    zeros(nch,1), zeros(nch,1),...
    cell(nch,1), zeros(nch,1), zeros(nch,1),...
    zeros(nch,1), zeros(nch,1),...
    cell(nch,1), cell(nch,1),...
    'VariableNames', {'meas','ch', 'chID', 'E_ec','E_eo','p','h','ttype','tstat','ES','p_FDR','h_FDR','v_ec','v_eo'});

%% statistical tests

for ch = 1:nch % iterate over channels

    % select variables
    if nch == 1
        v_ec = mean(data_ec.(stype).(frange).(meas),2);
        v_eo = mean(data_eo.(stype).(frange).(meas),2);
    else
        v_ec = data_ec.(stype).(frange).(meas)(:,ch);
        v_eo = data_eo.(stype).(frange).(meas)(:,ch);
    end

    % test for normality
    if lillietest(v_ec) || lillietest(v_eo)
        % nonparametric comparison - Wilcoxon signed rank test
        [p,h,s] = signrank(v_ec, v_eo);
        E_ec = median(v_ec);
        E_eo = median(v_eo);
        ttype = 'signrank';
        statvalue = s.zval;
        ES = abs(statvalue)/sqrt(length(v_ec));
    else
        % parametric comparison - paired t-test
        [h,p,~,s] = ttest(v_ec, v_eo);
        E_ec = mean(v_ec);
        E_eo = mean(v_eo);
        ttype = 'ttest';
        statvalue = s.tstat;
        ES_tmp = meanEffectSize(v_ec-v_eo,'Effect','cohen');
        ES = abs(ES_tmp.Effect('CohensD'));
    end

    % store outcomes
    output_table.meas{ch} = meas;
    output_table.ch{ch} = chlab{ch};
    output_table.chID(ch) = ch;
    output_table.E_ec(ch) = E_ec;
    output_table.E_eo(ch) = E_eo;
    output_table.p(ch) = p;
    output_table.h(ch) = h;
    output_table.ttype{ch} = ttype;
    output_table.tstat(ch) = statvalue;
    output_table.ES(ch) = ES;
    output_table.p_FDR(ch) = p;
    output_table.h_FDR(ch) = 0;
    output_table.v_ec{ch} = v_ec;
    output_table.v_eo{ch} = v_eo;
end

% sorting outcomes in descending order or significance (1-p)
output_table = sortrows(output_table, 'p', 'ascend');

end