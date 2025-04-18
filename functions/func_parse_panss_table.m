function [ table_out ] = func_parse_panss_table(table_in, pstate, panss)

% function to select physiological state (EO vs. EC) and PANSS score (GEN,
% NEG, POS, SUM) pairs.
% 
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 04/18/2025


if isempty(table_in)
    table_out = [];
    return
end

table_sign = table_in(table_in.h==1,:);
if isempty(table_sign)
    table_out = [];
    return
end

N = size(table_sign,1);
table_out = table(cell(N,1),cell(N,1),cell(N,1),cell(N,1),zeros(N,1),zeros(N,1),zeros(N,1),cell(N,1),...
    'VariableNames',{'meas','state','RSN','PANSS','r','p','h','ctype'});

for n = 1:N
    table_out.meas{n} = table_sign.meas{n};
    table_out.state{n} = pstate;
    table_out.RSN{n} = table_sign.rsn{n};
    table_out.PANSS{n} = panss;
    table_out.r(n) = table_sign.r(n);
    table_out.p(n) = table_sign.p(n);
    table_out.h(n) = table_sign.h(n);
    table_out.ctype{n} = table_sign.ctype{n};
end

end