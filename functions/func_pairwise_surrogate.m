function [output_table] = func_pairwise_surrogate(struct_data, chlab, meas)

% Function contrasts spectral indices to those obtained from surrogate time 
% series.
%
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 04/18/2025

%% preallocation and analysis constants

nch = length(chlab); % # of channels

output_table = table(cell(nch,1), cell(nch,1), zeros(nch,1),...
    zeros(nch,1), zeros(nch,1),...
    zeros(nch,1), zeros(nch,1),...
    cell(nch,1), zeros(nch,1),...
    zeros(nch,1), zeros(nch,1),...
    zeros(nch,1),...
    cell(nch,1), cell(nch,1), cell(nch,1),...
    'VariableNames', {'meas','ch', 'chID', 'E_tr','E_su','p','h','ttype','tstat','p_FDR','h_FDR','nZ','v_tr','v_su','v_zs'});

%% statistical tests

for ch = 1:nch % iterate over channels

    % select variables
    v_true = struct_data.(meas).true(:,ch);
    v_sur_all = squeeze(struct_data.(meas).sur(:,ch,:));
    v_sur = mean(v_sur_all,2);

    % test for normality
    if lillietest(v_true) || lillietest(v_sur)
        % nonparametric comparison - Wilcoxon signed rank test
        [p,h,s] = signrank(v_true, v_sur);
        E_tr = median(v_true);
        E_su = median(v_sur);
        ttype = 'signrank';
        statvalue = s.signedrank;
    else
        % parametric comparison - paired t-test
        [h,p,~,s] = ttest(v_true, v_sur);
        E_tr = mean(v_true);
        E_su = mean(v_sur);
        ttype = 'ttest';
        statvalue = s.tstat;
    end

    % z-testing
    ns = length(v_true);
    p_vec = zeros(ns,1);
    h_vec = zeros(ns,1);
    for subj = 1:ns
        [hz,pz] = ztest(v_true(subj),mean(v_sur_all(subj,:)),std(v_sur_all(subj,:)));
        p_vec(subj) = pz;
        h_vec(subj) = hz;
    end


    % store outcomes
    output_table.meas{ch} = meas;
    output_table.ch{ch} = chlab{ch};
    output_table.chID(ch) = ch;
    output_table.E_tr(ch) = E_tr;
    output_table.E_su(ch) = E_su;
    output_table.ttype{ch} = ttype;
    output_table.tstat(ch) = statvalue;
    output_table.p(ch) = p;
    output_table.h(ch) = h;
    output_table.p_FDR(ch) = p;
    output_table.h_FDR(ch) = 0;
    output_table.nZ(ch) = sum(h_vec);
    output_table.v_tr{ch} = v_true;
    output_table.v_su{ch} = v_sur;
    output_table.v_zs{ch} = p_vec;
end

% sorting outcomes in descending order or significance (1-p)
output_table = sortrows(output_table, 'p', 'ascend');

end