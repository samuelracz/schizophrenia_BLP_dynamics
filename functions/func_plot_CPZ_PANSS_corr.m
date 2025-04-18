function [ f ] = func_plot_CPZ_PANSS_corr(input_table)

% function to visualize correlations between PANSS scores and CPZ
% equivalent dose.
% 
% Author: F. Samuel Racz, The University of Texas at Austin
% email: fsr324@austin.utexas.edu
% last modified: 04/18/2025


f = figure('color','w','units','normalized','outerposition',[0.2,0,0.6,0.95]);

% PANSS-GEN
ind_panss = 1;
y_cpz = input_table.res_cpz{ind_panss};
y_panss = input_table.res_panss{ind_panss};
r = input_table.r(ind_panss);
p = input_table.p(ind_panss);
h = input_table.h(ind_panss);

subplot(221)
scatter(y_cpz,y_panss,140,'LineWidth',1.5,'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceColor',[0.7 0.9 1]);
ls = lsline;
ls.LineWidth = 3;
set(gca,'LineWidth',1,'FontSize',24,'box','on')
% annotations
a_panss = annotation('textbox');
a_panss.String = {['{\itr}=' num2str(r,'%.4f')],['{\itp}=' num2str(p,'%.4f')]};
if h == 1 
    a_panss.FontWeight = 'bold';
    ls.Color = [0.8 0 0];
else
    ls.Color = [0.2 0.2 0.2];
    ls.LineStyle = '--';
end
a_panss.FontSize = 28;
a_panss.Position = [0.34, 0.475, 0.2, 0.2];
a_panss.LineStyle = 'none';
xlabel('CPZ equivalent dose', 'FontSize', 28, 'visible', 'on')
ylabel('PANSS - General','FontSize',28,'visible','on')
title('CPZ ~ PANSS-GEN','FontSize',28)


% PANSS-NEG
ind_panss = 2;
y_cpz = input_table.res_cpz{ind_panss};
y_panss = input_table.res_panss{ind_panss};
r = input_table.r(ind_panss);
p = input_table.p(ind_panss);
h = input_table.h(ind_panss);

subplot(222)
scatter(y_cpz,y_panss,140,'LineWidth',1.5,'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceColor',[0.7 0.9 1]);
ls = lsline;
ls.LineWidth = 3;
set(gca,'LineWidth',1,'FontSize',24,'box','on')
% annotations
a_panss = annotation('textbox');
a_panss.String = {['{\itr}=' num2str(r,'%.4f')],['{\itp}=' num2str(p,'%.4f')]};
if h == 1 
    a_panss.FontWeight = 'bold';
    ls.Color = [0.8 0 0];
else
    ls.Color = [0.2 0.2 0.2];
    ls.LineStyle = '--';
end
a_panss.FontSize = 28;
a_panss.Position = [0.775, 0.475, 0.2, 0.2];
a_panss.LineStyle = 'none';
xlabel('CPZ equivalent dose', 'FontSize', 28, 'visible', 'on')
ylabel('PANSS - Negative','FontSize',28,'visible','on')
title('CPZ ~ PANSS-NEG','FontSize',28)


% PANSS-POS
ind_panss = 3;
y_cpz = input_table.res_cpz{ind_panss};
y_panss = input_table.res_panss{ind_panss};
r = input_table.r(ind_panss);
p = input_table.p(ind_panss);
h = input_table.h(ind_panss);

subplot(223)
scatter(y_cpz,y_panss,140,'LineWidth',1.5,'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceColor',[0.7 0.9 1]);
ls = lsline;
ls.LineWidth = 3;
set(gca,'LineWidth',1,'FontSize',24,'box','on')
% annotations
a_panss = annotation('textbox');
a_panss.String = {['{\itr}=' num2str(r,'%.4f')],['{\itp}=' num2str(p,'%.4f')]};
if h == 1 
    a_panss.FontWeight = 'bold';
    ls.Color = [0.8 0 0];
else
    ls.Color = [0.2 0.2 0.2];
    ls.LineStyle = '--';
end
a_panss.FontSize = 28;
a_panss.Position = [0.34, 0, 0.2, 0.2];
a_panss.LineStyle = 'none';
xlabel('CPZ equivalent dose', 'FontSize', 28, 'visible', 'on')
ylabel('PANSS - Positive','FontSize',28,'visible','on')
title('CPZ ~ PANSS-POS','FontSize',28)


% PANSS-SUM
ind_panss = 4;
y_cpz = input_table.res_cpz{ind_panss};
y_panss = input_table.res_panss{ind_panss};
r = input_table.r(ind_panss);
p = input_table.p(ind_panss);
h = input_table.h(ind_panss);

subplot(224)
scatter(y_cpz,y_panss,140,'LineWidth',1.5,'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceColor',[0.7 0.9 1]);
ls = lsline;
ls.LineWidth = 3;
set(gca,'LineWidth',1,'FontSize',24,'box','on')
% annotations
a_panss = annotation('textbox');
a_panss.String = {['{\itr}=' num2str(r,'%.4f')],['{\itp}=' num2str(p,'%.4f')]};
if h == 1 
    a_panss.FontWeight = 'bold';
    ls.Color = [0.8 0 0];
else
    ls.Color = [0.2 0.2 0.2];
    ls.LineStyle = '--';
end
a_panss.FontSize = 28;
a_panss.Position = [0.775, 0, 0.2, 0.2];
a_panss.LineStyle = 'none';
xlabel('CPZ equivalent dose', 'FontSize', 28, 'visible', 'on')
ylabel('PANSS - Total','FontSize',28,'visible','on')
title('CPZ ~ PANSS-SUM','FontSize',28)


end