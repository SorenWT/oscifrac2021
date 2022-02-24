% Figures
labs = {'a','b','c','d','e'};

%% Figure 1: schema pieces

% oscifrac schema

specs = load('~/Desktop/Data/Oscifrac_ISC_figures/212823_mmp_grpdata.mat_IRASA_specs.mat');
specs = specs.EEG;

frange = find((specs.freq > 2) & (specs.freq < 85));
freqs = specs.freq(frange);

pris = prism;

figure

set(gcf,'units','normalized','position',[0.15 0.3 0.7 0.7])

p = panel('no-manage-font');
p.pack('v',{1/2 1/2})
p(1).pack('h',{25 50 25})
p(2).pack('h',{50 50})

p(1,2).select()
loglog(freqs,10.^mean(log10(specs.mixd(frange,2)),2),'linewidth',5,'color',pris(6,:))
xlabel('Log power')
ylabel('Log frequency')
FixAxes(gca,18)
set(gca,'XTick',[],'YTick',[])

p(2,1).select()
plotdata = nanmean(log10(specs.osci(frange,2)),2);
plotdata(imag(plotdata)~=0) = NaN;
plotdata = medfilt1(plotdata,25,'omitnan');
loglog(freqs(4:end),10.^plotdata(4:end),'linewidth',5,'color',[0 0 1])
yl = get(gca,'YLim');
hold on
line([4 4],yl,'LineWidth',3,'linestyle','--','color','k')
line([8 8],yl,'LineWidth',3,'linestyle','--','color','k')
line([13 13],yl,'LineWidth',3,'linestyle','--','color','k')
line([30 30],yl,'LineWidth',3,'linestyle','--','color','k')
%patch([4 4 8 8],[yl(1) yl(2) yl(2) yl(1)],[0 0 0],'FaceAlpha',0.15,'EdgeColor','none')
%patch([13 13 30 30],[yl(1) yl(2) yl(2) yl(1)],[0 0 0],'FaceAlpha',0.15,'EdgeColor','none')
set(gca,'YLim',yl)
xlabel('Log power')
ylabel('Log frequency')
FixAxes(gca,18)
set(gca,'XTick',[],'YTick',[])

p(2,2).select()
loglog(freqs,10.^nanmean(log10(specs.frac(frange,2)),2),'linewidth',5,'color',[1 0 0]);
tmp = specs; tmp.frac = specs.frac(:,2);%tmp.frac = 10.^nanmean(log10(specs.frac),2);
beta = amri_sig_plawfit(tmp,[2 85]);
hold on
loglog(linspace(1,100,1000),10.^(-beta.Beta*log10(linspace(1,100,1000))+beta.Cons),...
    'linestyle','--','color','k','linewidth',3)

xlabel('Log power')
ylabel('Log frequency')
FixAxes(gca,18)
yl = get(gca,'YLim'); yl(2) = yl(2)+0.2*(abs(yl(1)-yl(2)));
set(gca,'XTick',[],'YTick',[],'XLim',[1 100],'YLim',yl)
p(1).marginbottom = 50;

savefig('Fig1_oscifrac_schema.fig'); export_fig('Fig1_oscifrac_schema.png','-m4')

% rest-task change figure

figure
p = panel('no-manage-font')
p.pack('h',{1/3 1/3 1/3})
set(gcf,'units','normalized','position',[0.15 0.4 0.7 0.6])

plotdata = zeros(362,1);
plotdata(chansind) = avgrest.data(5,chansind,2);
[p,~,ax,cbar] = brains4(plotdata,inflated,atlas,'panel',p,'panelindx',{2},'colormap',monored);
cla(ax(1)); cla(ax(3)); cla(ax(4));
Set_Clim(ax,[min(plotdata(plotdata>0)) max(plotdata(plotdata > 0))]);
delete(cbar)

plotdata = zeros(362,1);
plotdata(chansind) = taskoutputs{1}.data(5,chansind,2);
[p,~,ax,cbar]=brains4(plotdata,inflated,atlas,'panel',p,'panelindx',{1},'colormap',monored);
cla(ax(1)); cla(ax(3)); cla(ax(4));
Set_Clim(ax,[min(plotdata(plotdata>0)) max(plotdata(plotdata > 0))]);
delete(cbar)

plotdata = zeros(362,1);
plotdata(chansind) = taskoutputs{1}.data(5,chansind,2)-avgrest.data(5,chansind,2);
[p,~,ax,cbar]=brains4(plotdata,inflated,atlas,'panel',p,'panelindx',{3},'colormap',lkcmap2);
cla(ax(1)); cla(ax(3)); cla(ax(4));
Set_Clim(ax,[min(plotdata(plotdata~=0)) max(plotdata(plotdata~= 0))]);
delete(cbar)

p.marginright = 20;
savefig('Fig1b_rtchange_schema.fig'); export_fig('Fig1b_rtchange_schema.png','-m4')

%intra-subject variability schema

figure
p = panel('no-manage-font')
p.pack('h',{1/3 1/3 1/3})
set(gcf,'units','normalized','position',[0.15 0.4 0.7 0.6])

plotdata = zeros(362,1);
plotdata(chansind) = tresdata.data(1,chansind,2,1);
[p,~,ax,cbar] = brains4(plotdata,inflated,atlas,'panel',p,'panelindx',{1},'colormap',monored);
cla(ax(1)); cla(ax(3)); cla(ax(4));
Set_Clim(ax,[min(plotdata(plotdata>0)) max(plotdata(plotdata > 0))]);
delete(cbar)

plotdata = zeros(362,1);
plotdata(chansind) = tresdata.data(1,chansind,2,25);
[p,~,ax,cbar]=brains4(plotdata,inflated,atlas,'panel',p,'panelindx',{2},'colormap',monored);
cla(ax(1)); cla(ax(3)); cla(ax(4));
Set_Clim(ax,[min(plotdata(plotdata>0)) max(plotdata(plotdata > 0))]);
delete(cbar)

plotdata = zeros(362,1);
plotdata(chansind) = tresdata.data(1,chansind,2,50);
[p,~,ax,cbar]=brains4(plotdata,inflated,atlas,'panel',p,'panelindx',{3},'colormap',monored);
cla(ax(1)); cla(ax(3)); cla(ax(4));
Set_Clim(ax,[min(plotdata(plotdata~=0)) max(plotdata(plotdata~= 0))]);
delete(cbar)

p.marginright = 20;
savefig('Fig1b_intrasub_schema.fig'); export_fig('Fig1b_intrasub_schema.png','-m4')

%inter-subject variability schema


figure
p = panel('no-manage-font')
p.pack('v',{1/3 1/3 1/3})
set(gcf,'units','normalized','position',[0.3 0 0.6 1])

plotdata = zeros(362,1);
plotdata(chansind) = avgrest.data(1,chansind,2);
[p,~,ax,cbar] = brains4(plotdata,inflated,atlas,'panel',p,'panelindx',{1},'colormap',monored);
cla(ax(1)); cla(ax(3)); cla(ax(4));
Set_Clim(ax,[min(plotdata(plotdata>0)) max(plotdata(plotdata > 0))]);
delete(cbar)

plotdata = zeros(362,1);
plotdata(chansind) = avgrest.data(2,chansind,2);
[p,~,ax,cbar]=brains4(plotdata,inflated,atlas,'panel',p,'panelindx',{2},'colormap',monored);
cla(ax(1)); cla(ax(3)); cla(ax(4));
Set_Clim(ax,[min(plotdata(plotdata>0)) max(plotdata(plotdata > 0))]);
delete(cbar)

plotdata = zeros(362,1);
plotdata(chansind) = avgrest.data(3,chansind,2);
[p,~,ax,cbar]=brains4(plotdata,inflated,atlas,'panel',p,'panelindx',{3},'colormap',monored);
cla(ax(1)); cla(ax(3)); cla(ax(4));
Set_Clim(ax,[min(plotdata(plotdata~=0)) max(plotdata(plotdata~= 0))]);
delete(cbar)

p.marginright = 20;
savefig('Fig1b_intersub_schema.fig'); export_fig('Fig1b_intersub_schema.png','-m4')


%% figure 2: oscillatory vs fractal

hot2 = hot;
hot2 = hot2(1:192,:);

% v1: percent change

for i = 1:3
    figure
    
    p = panel('no-manage-font');
    
    set(gcf,'units','normalized','position',[0.1 0.2 0.8 0.8])
    
    p.pack('v',{50 50})
    p(1).pack('h',{10 40 40 10})
    p(2).pack('h',{30 40 30})
    p.fontsize = 14;
    
    % fractal
    tmp = [];
    for c = 1:2
        tmp(c,:) = horz(nanmedian(prcchange_sub{i}(:,:,c),1)).*horz(stats_resttask{i,c}.mask);
    end
    plotdata = zeros(nchans,1);
    plotdata(chansind) = sum(abs(tmp)/2,1);
    [p,~,ax1,cbar] = brains4(plotdata,inflated,atlas,'panel',p,'panelindx',{1,2},'mask',plotdata > 0);
    delete(cbar)
    p(1,2).margin = [10 15 5 5];
    p(1,2).de.margin = [30 15 3 3];
    
    %oscillatory
    tmp = [];
    for c = 1:5
        tmp(c,:) = horz(nanmedian(prcchange_sub{i}(:,:,c+2),1)).*horz(stats_resttask{i,c+2}.mask);
    end
    plotdata = zeros(nchans,1);
    plotdata(chansind) = sum(abs(tmp)/5,1);
    [p,~,ax2,cbar] = brains4(plotdata,inflated,atlas,'panel',p,'panelindx',{1,3},'mask',plotdata > 0);
    p(1,3).de.margin = [30 15 3 3];
    p(1,3).marginleft = 50;
    cbar.FontSize = 16;
    cbar.Label.String = 'Median absolute percent change';
    cbar.Label.FontSize = 20;
    
    Normalize_Clim([ax1 ax2],0)
    colormap(hot2)
    
    p(2,1).select()
    b1 = bar([1 2],prcchange(i,1:2),'FaceColor',palecol([1 0 0],0.3));
    hold on
    e = errorbar([1 2],prcchange(i,1:2),...
        prcchange(i,1:2)-prcchange_ci(i,1:2,1),prcchange_ci(i,1:2,2)-prcchange(i,1:2),...
        'LineStyle','none','LineWidth',2,'Color','k');
    ylabel('Fractal summed percent change')
    set(gca,'XTick',[1 2],'XTickLabel',measnames_short(1:2))
    FixAxes(gca,20)
    
    p(2,3).select()
    b2 = bar(1:5,prcchange(i,3:end),'FaceColor',palecol([0 0 1],0.3));
    hold on
    e = errorbar(1:5,prcchange(i,3:end),...
        prcchange(i,3:end)-prcchange_ci(i,3:end,1),prcchange_ci(i,3:end,2)-prcchange(i,3:end),...
        'LineStyle','none','LineWidth',2,'Color','k');
    set(gca,'XTick',1:5,'XTickLabel',measnames_short(3:end))
    ylabel('Oscillatory summed percent change')
    FixAxes(gca,20)
    p(2,3).marginleft = 42;
    p(2,2).marginleft = 5;
    
    Normalize_Ylim([p(2,1).axis, p(2,3).axis],0)
    
    p(2,2).select()
    [~,~,cbar]=imagesc_stars(vert(prcchange(i,:))-horz(prcchange(i,:)),pdif_prcchange(:,:,i),measnames_short,22);
    cbar.FontSize = 14;
    cbar.Label.String = 'Summed percent change difference (row - column)';
    cbar.Label.FontSize = 18;
    FixAxes(gca,16)
    title('Difference in rest-task change magnitude')
    
    p.marginleft = 34;
    p.marginbottom = 25;
    
    pause
    savefig(['Fig2' labs{i} '_prcchange_oscifrac_abs_mmpgroup.fig'])
    export_fig(['Fig2' labs{i} '_prcchange_oscifrac_abs_mmpgroup.png'],'-m5')
    close
end

% v2: pseudo effect size

for i = 1:3
    figure
    
    p = panel('no-manage-font');
    
    set(gcf,'units','normalized','position',[0.1 0.2 0.8 0.8])
    
    p.pack('v',{50 50})
    p(1).pack('h',{10 40 40 10})
    p(2).pack('h',{30 40 30})
    p.fontsize = 14;
    
    % fractal
    tmp = [];
    for c = 1:2
        tmp(c,:) = horz(nanmean(pseudoeffsize_sub{i}(:,:,c),1)).*horz(stats_resttask{i,c}.mask);
    end
    plotdata = zeros(nchans,1);
    plotdata(chansind) = sum(abs(tmp)/2,1);
    [p,~,ax1,cbar] = brains4(plotdata,inflated,atlas,'panel',p,'panelindx',{1,2},'mask',plotdata > 0);
    delete(cbar)
    p(1,2).margin = [10 15 5 5];
    p(1,2).de.margin = [30 15 3 3];
    
    %oscillatory
    tmp = [];
    for c = 1:5
        tmp(c,:) = horz(nanmean(pseudoeffsize_sub{i}(:,:,c+2),1)).*horz(stats_resttask{i,c+2}.mask);
    end
    plotdata = zeros(nchans,1);
    plotdata(chansind) = sum(abs(tmp)/5,1);
    [p,~,ax2,cbar] = brains4(plotdata,inflated,atlas,'panel',p,'panelindx',{1,3},'mask',plotdata > 0);
    p(1,3).de.margin = [30 15 3 3];
    p(1,3).marginleft = 50;
    cbar.FontSize = 16;
    cbar.Label.String = 'Average absolute effect size';
    cbar.Label.FontSize = 20;
    
    Normalize_Clim([ax1 ax2],0)
    colormap(hot)
    
    p(2,1).select()
    b1 = bar([1 2],pseudoeffsize(i,1:2),'FaceColor',palecol([1 0 0],0.3));
    hold on
    e = errorbar([1 2],pseudoeffsize(i,1:2),...
        pseudoeffsize(i,1:2)-pseudoeffsize_ci(i,1:2,1),pseudoeffsize_ci(i,1:2,2)-pseudoeffsize(i,1:2),...
        'LineStyle','none','LineWidth',2,'Color','k');
    ylabel('Fractal summed effect size')
    set(gca,'XTick',[1 2],'XTickLabel',measnames_short(1:2))
    FixAxes(gca,20)
    
    p(2,3).select()
    b2 = bar(1:5,pseudoeffsize(i,3:end),'FaceColor',palecol([0 0 1],0.3));
    hold on
    e = errorbar(1:5,pseudoeffsize(i,3:end),...
        pseudoeffsize(i,3:end)-pseudoeffsize_ci(i,3:end,1),pseudoeffsize_ci(i,3:end,2)-pseudoeffsize(i,3:end),...
        'LineStyle','none','LineWidth',2,'Color','k');
    set(gca,'XTick',1:5,'XTickLabel',measnames_short(3:end))
    ylabel('Oscillatory summed effect size')
    FixAxes(gca,20)
    p(2,3).marginleft = 42;
    p(2,2).marginleft = 5;
    
    Normalize_Ylim([p(2,1).axis, p(2,3).axis],0)
    
    p(2,2).select()
    [~,~,cbar]=imagesc_stars(vert(pseudoeffsize(i,:))-horz(pseudoeffsize(i,:)),pdif_pseudoeff(:,:,i),measnames_short,22);
    cbar.Label.String = 'Summed effect size difference (row - column)';
    cbar.FontSize = 14;
    cbar.Label.FontSize = 18;
    FixAxes(gca,16)
    
    p.marginleft = 34;
    p.marginbottom = 25;
    
    pause
    savefig(['Fig2' labs{i} '_effsize_oscifrac_abs.fig'])
    export_fig(['Fig2' labs{i} '_effsize_oscifrac_abs.png'],'-m5')
    close
end

% v3: Cohen's d


for i = 1:3
    figure
    
    p = panel('no-manage-font');
    
    set(gcf,'units','normalized','position',[0.1 0.2 0.8 0.8])
    
    p.pack('v',{50 50})
    p(1).pack('h',{10 40 40 10})
    p(2).pack('h',{30 40 30})
    p.fontsize = 14;
    
    % fractal
    tmp = [];
    for c = 1:2
        tmp(c,:) = horz(nanmean(cohend{i}(:,:,c),1)).*horz(stats_resttask{i,c}.mask);
    end
    plotdata = zeros(nchans,1);
    plotdata(chansind) = sum(abs(tmp)/2,1);
    [p,~,ax1,cbar] = brains4(plotdata,inflated,atlas,'panel',p,'panelindx',{1,2},'mask',plotdata > 0);
    delete(cbar)
    p(1,2).margin = [10 15 5 5];
    p(1,2).de.margin = [30 15 3 3];
    
    %oscillatory
    tmp = [];
    for c = 1:5
        tmp(c,:) = horz(nanmean(cohend{i}(:,:,c+2),1)).*horz(stats_resttask{i,c+2}.mask);
    end
    plotdata = zeros(nchans,1);
    plotdata(chansind) = sum(abs(tmp)/5,1);
    [p,~,ax2,cbar] = brains4(plotdata,inflated,atlas,'panel',p,'panelindx',{1,3},'mask',plotdata > 0);
    p(1,3).de.margin = [30 15 3 3];
    p(1,3).marginleft = 50;
    cbar.FontSize = 16;
    cbar.Label.String = 'Average absolute effect size';
    cbar.Label.FontSize = 20;
    
    Normalize_Clim([ax1 ax2],0)
    colormap(hot)
    
    p(2,1).select()
    b1 = bar([1 2],cohend_sum(i,1:2),'FaceColor',palecol([1 0 0],0.3));
    hold on
    e = errorbar([1 2],cohend_sum(i,1:2),...
        cohend_sum(i,1:2)-cohend_ci(i,1:2,1),cohend_ci(i,1:2,2)-cohend_sum(i,1:2),...
        'LineStyle','none','LineWidth',2,'Color','k');
    ylabel('Fractal summed effect size')
    set(gca,'XTick',[1 2],'XTickLabel',measnames_short(1:2))
    FixAxes(gca,20)
    
    p(2,3).select()
    b2 = bar(1:5,cohend_sum(i,3:end),'FaceColor',palecol([0 0 1],0.3));
    hold on
    e = errorbar(1:5,cohend_sum(i,3:end),...
        cohend_sum(i,3:end)-cohend_ci(i,3:end,1),cohend_ci(i,3:end,2)-cohend_sum(i,3:end),...
        'LineStyle','none','LineWidth',2,'Color','k');
    set(gca,'XTick',1:5,'XTickLabel',measnames_short(3:end))
    ylabel('Oscillatory summed effect size')
    FixAxes(gca,20)
    p(2,3).marginleft = 42;
    p(2,2).marginleft = 5;
    
    Normalize_Ylim([p(2,1).axis, p(2,3).axis],0)
    
    p(2,2).select()
    [~,~,cbar]=imagesc_stars(vert(cohend_sum(i,:))-horz(cohend_sum(i,:)),pdif_cohend(:,:,i),measnames_short,22);
    cbar.FontSize = 16;
    cbar.Label.String = 'Summed effect size difference (row - column)';
    cbar.Label.FontSize = 18;
    FixAxes(gca,16)

    
    p.marginleft = 34;
    p.marginbottom = 25;
    
    pause
    savefig(['Fig2' labs{i} '_cohen_oscifrac_abs.fig'])
    export_fig(['Fig2' labs{i} '_cohen_oscifrac_abs.png'],'-m5')
    close
end


%% figure 3: intercept vs ple

% v1: percent change

sp = spring; wi = winter;

difmap = customcolormap([0 0.5 1],[sp(1,:);[1 1 1];wi(end,:)]);
%difmap = customcolormap_preset('pink-white-green');

for i = 1:3
    figure
    
    p = panel();
    
    set(gcf,'units','normalized','position',[0.2 0 0.6 1])
    
    p.pack('v',{30 40 30})
    p(1).pack('h',{10 40 40 10})
    p(2).pack('h',{25 50 25})
    p(3).pack('h',{30 40 30})
    
    p.fontsize = 16;
    
    % intercept
    plotdata = zeros(nchans,1);
    plotdata(chansind) = squeeze(nanmedian(prcchange_sub{i}(:,:,1),1));
    mask = zeros(nchans,1);
    mask(chansind) = stats_resttask{i,1}.mask;
    [p,~,~,cbar] = brains4(plotdata,inflated,atlas,'panel',p,'panelindx',{1,2},'mask',mask);
    delete(cbar)
    p(1,2).margin = [10 15 5 5];
    p(1,2).de.margin = [30 15 3 3];
    
    % PLE
    
    plotdata = zeros(nchans,1);
    plotdata(chansind) = squeeze(nanmedian(prcchange_sub{i}(:,:,2),1));
    mask = zeros(nchans,1);
    mask(chansind) = stats_resttask{i,2}.mask;
    [p,~,~,cbar] = brains4(plotdata,inflated,atlas,'panel',p,'panelindx',{1,3},'mask',mask);
    p(1,3).de.margin = [30 15 3 3];
    p(1,3).marginleft = 50;
    cbar.Label.String = 'Rest-task percent change (%)';
    cbar.Label.FontSize = 18;
    
    Normalize_Clim(gcf,1)
    colormap(lkcmap2)
    
    p(1).marginbottom = 40;
    
    % difference between intercept and PLE
    mask = zeros(nchans,1);
    mask(chansind) = pprc_intvsple(i,:)<0.05;
    plotdata = zeros(nchans,1);
    plotdata(chansind) = squeeze(nanmedian(prcchange_sub{i}(:,:,1)-prcchange_sub{i}(:,:,2),1));
    [p,~,ax,cbar] = brains4(plotdata,inflated,atlas,'panel',p,'panelindx',{2,2},'mask',mask,'colormap',difmap);
    cbar.Label.String = 'Rest-task percent change difference (%)';
    cbar.Label.FontSize = 18;
    p(2,2).de.margin = [20 12 3 3];
    Normalize_Clim(ax,1)
    
    % unique activations for intercept
    mask = zeros(nchans,1);
    mask(chansind) = stats_resttask{i,1}.mask & (~stats_resttask{i,2}.mask);
    plotdata = mask;
    
    p(2,1).pack('v',{20 60 20})
    [p,~,~,cbar] = brains4(plotdata,inflated,atlas,'panel',p,'panelindx',{2,1,2},'mask',mask,'colormap',difmap);
    delete(cbar)
    p(2,1).margintop = 30;
    p(2,1).de.margin = [8 8 3 3];
    
    % unique activations for PLE
    mask = zeros(nchans,1);
    mask(chansind) = (~stats_resttask{i,1}.mask) & stats_resttask{i,2}.mask;
    plotdata = -mask;
    
    p(2,3).pack('v',{20 60 20})
    [p,~,~,cbar] = brains4(plotdata,inflated,atlas,'panel',p,'panelindx',{2,3,2},'mask',mask,'colormap',difmap);
    delete(cbar)
    p(2,3).margintop = 30;
    p(2,3).de.margin = [8 8 3 3];
    
    
    p(3,2).select()
    %b1 = bar(1,prcchange(i,1),'FaceColor',difmap(end-9,:));
    %hold on
    %b2 = bar(2,prcchange(i,2),'FaceColor',difmap(10,:));
    %e = errorbar([1 2],prcchange(i,1:2),...
    %    prcchange(i,1:2)-prcchange_ci(i,1:2,1),prcchange_ci(i,1:2,2)-prcchange(i,1:2),...
    %    'LineStyle','none','LineWidth',2,'Color','k');
    tmpp = signrank(median(abs(prcchange_sub{i}(:,:,1)),2),median(abs(prcchange_sub{i}(:,:,2)),2));
    pairsamp_lineplot([median(abs(prcchange_sub{i}(:,:,1)),2),median(abs(prcchange_sub{i}(:,:,2)),2)],[difmap(end,:);difmap(1,:)],palecol([0 0 1;1 0 0;0 0 0]))
    %ylabel('Summed absolute percent change')
    ylabel('Median parcel absolute percent change')
    set(gca,'XTick',[1 2],'XTickLabel',measnames_short(1:2))
    %FixAxes(gca,14)
    %sigstar({[1 2]},pdif_prcchange(2,1),0,18);
    sigstar({[1 2]},tmpp,0,18);
    
    p.marginleft = 20;
    p.marginbottom = 12;
    
    pause
    savefig(['Fig3' labs{i} '_prcchange_intple_pairsplot_mmpgroup.fig']);
    export_fig(['Fig3' labs{i} '_prcchange_intple_pairsplot_mmpgroup.png'],'-m4')
    close
end

% v2: pseudo effect size

%difmap = customcolormap_preset('pink-white-green');

for i = 1:3
    figure
    
    p = panel();
    
    set(gcf,'units','normalized','position',[0.2 0 0.6 1])
    
    p.pack('v',{30 40 30})
    p(1).pack('h',{10 40 40 10})
    p(2).pack('h',{25 50 25})
    p(3).pack('h',{30 40 30})
    
    p.fontsize = 16;
    
    % intercept
    plotdata = zeros(nchans,1);
    plotdata(chansind) = squeeze(mean(pseudoeffsize_sub{i}(:,:,1),1));
    mask = zeros(nchans,1);
    mask(chansind) = stats_resttask{i,1}.mask;
    [p,~,~,cbar] = brains4(plotdata,inflated,atlas,'panel',p,'panelindx',{1,2},'mask',mask);
    delete(cbar)
    p(1,2).margin = [10 15 5 5];
    p(1,2).de.margin = [30 15 3 3];
    
    % PLE
    
    plotdata = zeros(nchans,1);
    plotdata(chansind) = squeeze(mean(pseudoeffsize_sub{i}(:,:,2),1));
    mask = zeros(nchans,1);
    mask(chansind) = stats_resttask{i,2}.mask;
    [p,~,~,cbar] = brains4(plotdata,inflated,atlas,'panel',p,'panelindx',{1,3},'mask',mask);
    p(1,3).de.margin = [30 15 3 3];
    p(1,3).marginleft = 50;
    cbar.Label.String = 'Pseudo effect size';
    cbar.Label.FontSize = 18;
    
    Normalize_Clim(gcf,1)
    colormap(lkcmap2)
    
    p(1).marginbottom = 40;
    
    % difference between intercept and PLE
    mask = zeros(nchans,1);
    mask(chansind) = peff_intvsple(i,:)<0.05;
    plotdata = zeros(nchans,1);
    plotdata(chansind) = squeeze(mean(pseudoeffsize_sub{i}(:,:,1),1)-mean(pseudoeffsize_sub{i}(:,:,2),1));
    [p,~,ax,cbar] = brains4(plotdata,inflated,atlas,'panel',p,'panelindx',{2,2},'mask',mask,'colormap',difmap);
    cbar.Label.String = 'Effect size difference';
    cbar.Label.FontSize = 14;
    p(2,2).de.margin = [20 12 3 3];
    Normalize_Clim(ax,1)
    
    % unique activations for intercept
    mask = zeros(nchans,1);
    mask(chansind) = stats_resttask{i,1}.mask & (~stats_resttask{i,2}.mask);
    plotdata = mask;
    
    p(2,1).pack('v',{20 60 20})
    [p,~,~,cbar] = brains4(plotdata,inflated,atlas,'panel',p,'panelindx',{2,1,2},'mask',mask,'colormap',difmap);
    delete(cbar)
    p(2,1).margintop = 30;
    p(2,1).de.margin = [8 8 3 3];
    
    % unique activations for PLE
    mask = zeros(nchans,1);
    mask(chansind) = (~stats_resttask{i,1}.mask) & stats_resttask{i,2}.mask;
    plotdata = -mask;
    
    p(2,3).pack('v',{20 60 20})
    [p,~,~,cbar] = brains4(plotdata,inflated,atlas,'panel',p,'panelindx',{2,3,2},'mask',mask,'colormap',difmap);
    delete(cbar)
    p(2,3).margintop = 30;
    p(2,3).de.margin = [8 8 3 3];
    
    
    p(3,2).select()
%     b1 = bar(1,pseudoeffsize(i,1),'FaceColor',difmap(end-9,:));
%     hold on
%     b2 = bar(2,pseudoeffsize(i,2),'FaceColor',difmap(10,:));
%     e = errorbar([1 2],pseudoeffsize(i,1:2),...
%         pseudoeffsize(i,1:2)-pseudoeffsize_ci(i,1:2,1),pseudoeffsize_ci(i,1:2,2)-pseudoeffsize(i,1:2),...
%         'LineStyle','none','LineWidth',2,'Color','k');
    tmpp = signrank(median(abs(pseudoeffsize_sub{i}(:,:,1)),2),median(abs(pseudoeffsize_sub{i}(:,:,2)),2));
    pairsamp_lineplot([median(abs(pseudoeffsize_sub{i}(:,:,1)),2),median(abs(pseudoeffsize_sub{i}(:,:,2)),2)],[difmap(end,:);difmap(1,:)],palecol([0 0 1;1 0 0;0 0 0]))
    ylabel('Median parcel absolute effect size')
    set(gca,'XTick',[1 2],'XTickLabel',measnames_short(1:2))
    FixAxes(gca,14)
    sigstar({[1 2]},tmpp,0,18);
    %sigstar({[1 2]},pdif_pseudoeff(2,1),0,18);
    
    p.marginleft = 20;
    p.marginbottom = 12;
    
    pause
    savefig(['Fig3' labs{i} '_effsize_intple_pairsplot.fig']);
    export_fig(['Fig3' labs{i} '_effsize_intple_pairsplot.png'],'-m4')
    close
end

%% figure 4: Inter- and intra-subject CV

% v1. HCP dataset

monored = customcolormap_preset('red-white-blue');
monoblue = flipud(monored(1:128,:));
monored = monored(129:end,:);
monopink = difmap(129:end,:);
monogreen = flipud(difmap(1:128,:));

cvs = {'cvinter','cvintra','cvrat'};
cv_labs = {'Inter-subject CV','Intra-subject CV','Normalized inter-subject CV'};

for i = 1:3
    cvdata = eval(cvs{i});
    cv_ci = eval([cvs{i} '_ci;']);
    cv_pdiff = eval([cvs{i} '_pdiff;']);
    
    figure
    
    p = panel('no-manage-font');
    
    set(gcf,'units','normalized','position',[0.15 0 0.7 1])
    
    p.pack('v',{30 20 50})
    p.marginleft = 20; p.de.marginleft = 15; p.marginright = 10; p.de.marginright = 5;
    p.marginbottom = 22; p.de.marginbottom = 15;
    p(1).pack('h',{5 25 70})
    p(1,3).pack('v',{70 30})
    p(1,3,1).pack('h',{25 25 25 25})
    p(3).pack('h',{5 25 10 60})
    p(3,2).pack('v',{50 50})
    p(1,2).marginleft = 16; 
    p(3,2).marginleft = 10;
    p(1,3).marginleft = 50;  p(1,3).de.marginleft = 5;
    
    p(2).select()
    b1 = bar(1:2,cvdata(1:2),'FaceColor',palecol([1 0 0],0.3));
    hold on
    b2 = bar(4:8,cvdata(3:end),'FaceColor',palecol([0 0 1],0.3));
    set(gca,'XTick',[1 2 4:8],'XTickLabel',measnames_short)
    e = errorbar([1 2 4:8],cvdata,...
        vert(cvdata)-cv_ci(:,1),cv_ci(:,2)-vert(cvdata),...
        'LineStyle','none','LineWidth',2,'Color','k');
    ylabel(cv_labs{i})
    FixAxes(gca,14)
    
    if i == 1
        plotdata = nanmean(geocv(avgrest.data(:,:,1:2),1),3);
    elseif i == 2
        plotdata = nanmean(nanmedian(geocv(tresdata.data(:,:,1:2,:),4),1),3);
    elseif i == 3
        plotdata = nanmean(geocv(avgrest.data(:,:,1:2),1),3)./nanmean(nanmedian(geocv(tresdata.data(:,:,1:2,:),4),1),3);
    end
    mask = ones(nchans,1);
    mask(except(1:nchans,chansind)) = 0;
    [p,~,ax,cbar] = brains4(plotdata,inflated,atlas,'panel',p,'panelindx',{1 2},'mask',mask,'colormap',monored)
    cbar.Label.String = 'CV'; cbar.Label.FontSize = 12; cbar.FontSize = 12;
    p(1,2).de.margin = [20 8 5 5];
        Set_Clim(ax,[min(plotdata) max(plotdata)]);

    
    if i == 1
        plotdata = nanmean(geocv(avgrest.data(:,:,3:end),1),3);
    elseif i == 2
        plotdata = nanmean(nanmedian(geocv(tresdata.data(:,:,3:end,:),4),1),3);
    elseif i == 3
        plotdata = nanmean(geocv(avgrest.data(:,:,3:end),1),3)./nanmean(nanmedian(geocv(tresdata.data(:,:,3:end,:),4),1),3);
    end
    mask = ones(nchans,1);
    mask(except(1:nchans,chansind)) = 0;
    extrafig = figure;
    [~,~,ax,cbar] = brains4(plotdata,inflated,atlas,'mask',mask,'colormap',monoblue)
    cbar.Label.String = 'CV'; cbar.Label.FontSize = 12; cbar.FontSize = 12;
    for q = 1:length(ax)
        p(1,3,1,q).select(ax(q))
        ax(q) = gca;
    end
    p(1,3,1).de.margin = [40 5 5 5];
    %Normalize_Clim(ax,0)
    Set_Clim(ax,[min(plotdata) max(plotdata)]);
    delete(extrafig)
    
    
    if i == 1
        plotdata = nanmean(geocv(avgrest.data(:,:,1),1),3);
    elseif i == 2
        plotdata = nanmean(nanmedian(geocv(tresdata.data(:,:,1,:),4),1),3);
    elseif i == 3
        plotdata = nanmean(geocv(avgrest.data(:,:,1),1),3)./nanmean(nanmedian(geocv(tresdata.data(:,:,1,:),4),1),3);
    end
    mask = ones(nchans,1);
    mask(except(1:nchans,chansind)) = 0;
    [p,~,ax,cbar] = brains4(plotdata,inflated,atlas,'panel',p,'panelindx',{3 2 1},'mask',mask,'colormap',monopink)
    cbar.Label.String = 'CV'; cbar.Label.FontSize = 12; cbar.FontSize = 12;
        Set_Clim(ax,[min(plotdata) max(plotdata)]);

    
    if i == 1
        plotdata = nanmean(geocv(avgrest.data(:,:,2),1),3);
    elseif i == 2
        plotdata = nanmean(nanmedian(geocv(tresdata.data(:,:,2,:),4),1),3);
    elseif i == 3
        plotdata = nanmean(geocv(avgrest.data(:,:,2),1),3)./nanmean(nanmedian(geocv(tresdata.data(:,:,2,:),4),1),3);
    end
    mask = ones(nchans,1);
    mask(except(1:nchans,chansind)) = 0;
    [p,~,ax,cbar] = brains4(plotdata,inflated,atlas,'panel',p,'panelindx',{3 2 2},'mask',mask,'colormap',monogreen)
    cbar.Label.String = 'CV'; cbar.Label.FontSize = 12; cbar.FontSize = 12;
    Set_Clim(ax,[min(plotdata) max(plotdata)]);

    
    p(3,4).select()
    [~,~,cbar]=imagesc_stars(cvdata'-cvdata,cv_pdiff,measnames_short);
    FixAxes(gca,14)
    cbar.Label.String = 'CV difference';
    cbar.Label.FontSize = 16; 
    
    pause
    savefig(['Fig5' labs{i} '_' cvs{i} '_mmpgroup.fig']); 
    export_fig(['Fig5' labs{i} '_' cvs{i} '_mmpgroup.png'],'-m4');
    close
end


%% Figure 6: 3D and 2D variability plots
figure

set(gcf,'units','normalized','position',[0.3 0 0.4 1])

p = panel('no-manage-font');

p.pack('v',{1/2 1/2})

p.marginleft = 20;

p(1).select()
cols = cat(1,palecol([1 0 0],0.3),palecol([1 0 0],0.3),...
        repmat(palecol([0 0 1],0.3),5,1));
msum = mean(prcchange,1);
%msum = mean(pseudoeffsize,1);
sp = spheres_plot(cvinter,cvintra,msum,sqrt(vecnorm(NormOntoRange([cvrat',cvintra',msum'],[0.1 1.1]),2,2))*0.09,cols);
% move around and set camlight
material([0.65 0.55 0.4 20 1])
legend(measnames)
xlabel('Inter-subject variability')
ylabel('Intra-subject variability')
zlabel('Task-responsiveness')
FixAxes(gca,14)
set(gca,'XGrid','on','YGrid','on','ZGrid','on')
p(1).margin = [30 30 5 5];

p(2).select()

meanintra = (zscore(msum)+zscore(cvintra))/2;
scatter(cvinter(1:2),meanintra(1:2),144,palecol([1 0 0],0.3),'filled')
hold on
scatter(cvinter(3:7),meanintra(3:7),144,palecol([0 0 1],0.3),'filled')
labelpoints(cvinter,meanintra,measnames_short,'fontsize',12,'buffer',0.15)
FixAxes(gca,16)
xlabel('Inter-subject CV')
ylabel('Aggregated intra-subject variability')

p(1).marginbottom = 40;
p.marginbottom = 20; p.marginleft = 24;

savefig('Fig6_intra-inter.fig'); export_fig('Fig6_intra-inter.png','-m5')

% old version using pseudoeffsize

scatter(cvinter(1:2),msum(1:2),144,palecol([1 0 0],0.3),'filled')
hold on
scatter(cvinter(3:7),msum(3:7),144,palecol([0 0 1],0.3),'filled')
labelpoints(cvrat,mean(pseudoeffsize,1)',measnames_short,'fontsize',12)
FixAxes(gca,16)
xlabel('Normalized inter-subject variability')
ylabel('Normalized task-responsiveness')
% 
scatter(cvinter(1:2),cvintra(1:2),144,palecol([1 0 0],0.3),'filled')
hold on
scatter(cvinter(3:7),cvintra(3:7),144,palecol([0 0 1],0.3),'filled')
labelpoints(cvrat,mean(pseudoeffsize,1)',measnames_short,'fontsize',12)
FixAxes(gca,16)
xlabel('Normalized inter-subject variability')
ylabel('Normalized task-responsiveness')


scatter(cvrat(1:2),mean(pseudoeffsize(:,1:2),1)',144,palecol([1 0 0],0.3),'filled')
hold on
scatter(cvrat(3:7),mean(pseudoeffsize(:,3:7),1)',144,palecol([0 0 1],0.3),'filled')
labelpoints(cvrat,mean(pseudoeffsize,1)',measnames_short,'fontsize',12,'buffer',0.15)
FixAxes(gca,14)
xlabel('Normalized inter-subject variability')
ylabel('Normalized task-responsiveness')


savefig('Fig6_intra-inter_norm.fig'); export_fig('FigS6_intra-inter_norm.png','-m4')