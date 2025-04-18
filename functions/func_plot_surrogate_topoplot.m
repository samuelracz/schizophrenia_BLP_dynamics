function [ f ] = func_plot_surrogate_topoplot(table_hc, struct_hc, table_sz, struct_sz, chanlocs, meas, pstate)

if strcmp(meas,'mean')
    str_title = ['\mu(\alpha-BLP^{' pstate '}_{mixd})'];
elseif strcmp(meas,'var')
    str_title = ['\sigma(\alpha-BLP^{' pstate '}_{mixd})'];

else
    error('wrong input (mean/var).')
end

nch = length(struct_hc.chlab);

% ensure that data is sorted in original channel order
table_hc = sortrows(table_hc,'chID','ascend');
table_sz = sortrows(table_sz,'chID','ascend');

% find global minima and maxima for setting the color scale
min_glob = min([table_hc.nZ; table_sz.nZ]);
max_glob = max([table_hc.nZ; table_sz.nZ]);

%% plotting
f = figure('color','w','units','normalized','outerposition',[0.05 0.05 0.9 0.9]);

%% HC topoplots

% low-frequency regime (1-4 Hz)
subplot(2,2,1)
mu_lo = table_hc.nZ;
topoplot(mu_lo, chanlocs, 'emarker', {'.','k',20,2}, 'whitebk', 'on');
clim([min_glob max_glob])
colormap parula
cb1 = colorbar;
cb1.Label.String = '# of subj. > chance';
cb1.Label.Rotation = 270;
cb1.FontSize = 16;
cb1.Label.FontSize = 28;
cb1.Label.VerticalAlignment = 'bottom';
title(['Healthy, ' str_title],'FontSize',28)

%% HC line plots
subplot(223)
p1 = plot(mean(struct_hc.(meas).true,1),'ro-','LineWidth',2);
hold on
p2 = plot(mean(mean(struct_hc.(meas).sur,3),1),'b-','LineWidth',1.5);
plot(mean(mean(struct_hc.(meas).sur,3),1) + mean(std(struct_hc.(meas).sur,[],3),1),'b:','LineWidth',1)
plot(mean(mean(struct_hc.(meas).sur,3),1) - mean(std(struct_hc.(meas).sur,[],3),1),'b:','LineWidth',1)
set(gca,'XTick',1:nch,'XTickLabel',struct_hc.chlab,'xlim',[0 nch+1],'ylim',[0.1 0.4],'FontSize',12)
legend([p1 p2],{'True signal','Surrogates'},'FontSize',14,'location','northeast')
xlabel('Channels','FontSize',24)
ylabel(str_title,'FontSize',24)
title(['Healthy, ' str_title],'FontSize',28)


%% SZ topoplots

% low-frequency regime (1-4 Hz)
subplot(2,2,2)
mu_lo = table_sz.nZ;
topoplot(mu_lo, chanlocs, 'emarker', {'.','k',20,2}, 'whitebk', 'on');
clim([min_glob max_glob])
colormap parula
cb4 = colorbar;
cb4.Label.String = '# of subj. > chance';
cb4.Label.Rotation = 270;
cb4.FontSize = 16;
cb4.Label.FontSize = 28;
cb4.Label.VerticalAlignment = 'bottom';
title(['Schizophrenia, ' str_title],'FontSize',28)

%% HC line plots
subplot(224)
p3 = plot(mean(struct_sz.(meas).true,1),'ro-','LineWidth',2);
hold on
p4 = plot(mean(mean(struct_sz.(meas).sur,3),1),'b-','LineWidth',1.5);
plot(mean(mean(struct_sz.(meas).sur,3),1) + mean(std(struct_sz.(meas).sur,[],3),1),'b:','LineWidth',1)
plot(mean(mean(struct_sz.(meas).sur,3),1) - mean(std(struct_sz.(meas).sur,[],3),1),'b:','LineWidth',1)
set(gca,'XTick',1:nch,'XTickLabel',struct_sz.chlab,'xlim',[0 nch+1],'ylim',[0.1 0.4],'FontSize',12)
legend([p3 p4],{'True signal','Surrogates'},'FontSize',14,'location','northeast')
xlabel('Channels','FontSize',24)
ylabel(str_title,'FontSize',24)
title(['Schizophrenia, ' str_title],'FontSize',28)

end