function [ f ] = func_plot_phase_cycle(table_mean, table_var)

table_mean = sortrows(table_mean,'freqID','ascend');
table_var = sortrows(table_var,'freqID','ascend');
N = size(table_mean,1);
freq = cell2mat(table_mean.freq);

% get variables
m_hc_m = zeros(N,1);
m_hc_sem = zeros(N,1);
m_sz_m = zeros(N,1);
m_sz_sem = zeros(N,1);

v_hc_m = zeros(N,1);
v_hc_sem = zeros(N,1);
v_sz_m = zeros(N,1);
v_sz_sem = zeros(N,1);

for n = 1:N
    m_hc_m(n) = mean(table_mean.v_hc{n});
    m_hc_sem(n) = std(table_mean.v_hc{n})/sqrt(length(table_mean.v_hc{n}));
    m_sz_m(n) = mean(table_mean.v_sz{n});
    m_sz_sem(n) = std(table_mean.v_sz{n})/sqrt(length(table_mean.v_sz{n}));
    v_hc_m(n) = mean(table_var.v_hc{n});
    v_hc_sem(n) = std(table_var.v_hc{n})/sqrt(length(table_var.v_hc{n}));
    v_sz_m(n) = mean(table_var.v_sz{n});
    v_sz_sem(n) = std(table_var.v_sz{n})/sqrt(length(table_var.v_sz{n}));
end

f = figure('color','w','units','normalized','outerposition',[0.1,0.1,0.8,0.6]);
subplot(121)
hold on
fill([freq; flipud(freq)],[(m_hc_m-m_hc_sem); flipud(m_hc_m+m_hc_sem)],'b','FaceColor','b','EdgeAlpha',0,'FaceAlpha',0.25);
p1 = plot(freq,m_hc_m,'bs-','LineWidth',2);

fill([freq; flipud(freq)],[(m_sz_m-m_sz_sem); flipud(m_sz_m+m_sz_sem)],'r','FaceColor','r','EdgeAlpha',0,'FaceAlpha',0.25);
p2 = plot(freq,m_sz_m,'rs-','LineWidth',2);

ymin_m = min(get(gca,'ylim'));
ymax_m = max(get(gca,'ylim'));
for n = 1:n
    if table_mean.h_FDR(n)==1
        plot(freq(n),ymax_m,'k*','MarkerSize',10,'LineWidth',1.5);
    end
end

set(gca,'XScale','log','box','on','xlim',[freq(1) freq(end)],'ylim',[ymin_m ymax_m],'FontSize',16)
legend([p1,p2],{'HC','SZ'},'FontSize',24,'location','southwest')
xlabel('Log(frequency)','FontSize',20)
ylabel('Log(amplitude)','FontSize',20)
title('Global \mu(\it{P_{mixd}(f))}','FontSize',24)

subplot(122)
hold on
fill([freq; flipud(freq)],[(v_hc_m-v_hc_sem); flipud(v_hc_m+v_hc_sem)],'b','FaceColor','b','EdgeAlpha',0,'FaceAlpha',0.25);
p1 = plot(freq,v_hc_m,'bo-','LineWidth',2);

fill([freq; flipud(freq)],[(v_sz_m-v_sz_sem); flipud(v_sz_m+v_sz_sem)],'r','FaceColor','r','EdgeAlpha',0,'FaceAlpha',0.25);
p2 = plot(freq,v_sz_m,'ro-','LineWidth',2);

ymin_v = min(get(gca,'ylim'));
ymax_v = max(get(gca,'ylim'));
for n = 1:n
    if table_var.h_FDR(n)==1
        plot(freq(n),ymax_v,'k*','MarkerSize',10,'LineWidth',1.5);
    end
end

set(gca,'XScale','log','box','on','xlim',[freq(1) freq(end)],'ylim',[ymin_v 1.1*ymax_v],'FontSize',16)
legend([p1,p2],{'HC','SZ'},'FontSize',24,'location','southeast')
xlabel('Log(frequency)','FontSize',20)
ylabel('Log(amplitude)','FontSize',20)
title('Global \sigma(\it{P_{mixd}(f))}','FontSize',24)



end

%% helper functions
function [] = helper_plot_var(x,y,d,col)

hold on
plot([x-d x+d],[median(y) median(y)],'-','color',col,'LineWidth',2); % plot median
fill([x-d x-d x+d x+d],[quantile(y,0.25) quantile(y,0.75) quantile(y,0.75) quantile(y,0.25)],col,'EdgeColor',col,'LineWidth',1,'FaceAlpha',0); % plot IQR
plot([x x],[quantile(y,0.75) quantile(y,0.95)],'-','color',col,'LineWidth',1)
plot([x x],[quantile(y,0.05) quantile(y,0.25)],'-','color',col,'LineWidth',1)


end