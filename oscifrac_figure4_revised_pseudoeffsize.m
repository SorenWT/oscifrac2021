
addpath(genpath('/project/def-gnorthof/sorenwt/MATLAB'))
addpath('/project/def-gnorthof/sorenwt/fieldtrip-master')
ft_defaults
addpath([toolboxdir('signal'),'/signal'])

load('settings_camcan_1Hz.mat')

cd /home/soren/Documents/camcan/source/task/epoched/

files = dir('*tf_source.mat');

m = matfile(files(1).name);

fields = {'mixd','osci','frac'};

for c = 1:length(fields)
    meandata.(fields{c}) = m.(fields{c});
end

%% Reading in data
donefiles = [];
for i = 1:length(files)
    try
        %if ~ismember({files(i).name},donefiles)
            m = matfile(files(i).name);
            for c = 1:3
                tmp = m.(fields{c});
                tmp.fourierspctrm(find(tmp.fourierspctrm < 0)) = NaN;
                meandata.(fields{c}).fourierspctrm(i,:,:,:) = nanmean(tmp.fourierspctrm,1);
            end
            donefiles = [donefiles {files(i).name}];
        %end
    catch
    end
end
save('oscifrac_data_revision1.mat','-v7.3')

for c = 1:length(fields)
   meandata.(fields{c}).fourierspctrm = permute(meandata.(fields{c}).fourierspctrm,[1 3 2 4]); 
   meandata.(fields{c}).fourierspctrm_dimord = meandata.(fields{c}).dimord;
   %meandata.(fields{c}).grad = settings.datasetinfo.grad;
   meandata.(fields{c}).fourierspctrm = meandata.(fields{c}).fourierspctrm(1:length(files),:,:,:);
end


%% Analysis and statistics

for c = 1:3
    post = find(meandata.(fields{c}).time >= 0);
    meanpost.(fields{c}) = meandata.(fields{c});
    meanpost.(fields{c}).fourierspctrm = meandata.(fields{c}).fourierspctrm(:,:,:,post);
    meanpost.(fields{c}).time = meandata.(fields{c}).time(post);
    
    bl = find(meandata.(fields{c}).time < 0);
    %bl = bl((end-length(post)+1):end);
    meanbl.(fields{c}) = meandata.(fields{c});
    meanbl.(fields{c}).fourierspctrm = meandata.(fields{c}).fourierspctrm(:,:,:,bl);
    meanbl.(fields{c}).time = meandata.(fields{c}).time(bl);
end

meanblflat = meanpost; 

for i = 1:length(fields)
   meanblflat.(fields{i}).fourierspctrm = repmat(nanmean(meanbl.(fields{i}).fourierspctrm,4),...
       1,1,1,length(meanpost.(fields{i}).time));
end

post = cell(1,size(meanpost.frac.fourierspctrm,1));
bl = post;
postint = post;
blint = bl;
postindx = find(meandata.(fields{c}).time >= 0);
    blindx = find(meandata.mixd.time < 0);
    %blindx = blindx((end-length(postindx)+1):end);


frange = intersect(find(meanpost.frac.freq > settings.tfparams.fbands{2}(1)),find(meanpost.frac.freq < settings.tfparams.fbands{end}(2)));
parfor i = 1:size(meanpost.frac.fourierspctrm,1)
    for ii = 1:size(meanpost.frac.fourierspctrm,2)
        for iii = 1:size(meanpost.frac.fourierspctrm,4)
            tmp = log10(meanpost.frac.freq(frange));
            pow = log10(squeeze(meanpost.frac.fourierspctrm(i,ii,frange,iii)));
            lintmp = vert(linspace(tmp(1),tmp(end),length(tmp)));
            pow = interp1(tmp,pow,lintmp);
            p = polyfit(lintmp,pow,1);
            post{i}(1,ii,iii) = -p(1);
            postint{i}(1,ii,iii) = p(2);
        end
        for iii = 1:size(meanbl.frac.fourierspctrm,4)
            tmp = log10(meanbl.frac.freq(frange));
            pow = log10(squeeze(meanbl.frac.fourierspctrm(i,ii,frange,iii)));
            lintmp = vert(linspace(tmp(1),tmp(end),length(tmp)));
            pow = interp1(tmp,pow,lintmp);
            p = polyfit(lintmp,pow,1);
            bl{i}(1,ii,iii) = -p(1);
            blint{i}(1,ii,iii) = p(2);
        end
    end
end
save('oscifrac_data_revision1.mat','-v7.3')

pledata.post = cat(1,post{:});
pledata.bl = cat(1,bl{:});

intdata.post = cat(1,postint{:});
intdata.bl = cat(1,blint{:});


load('~/Documents/oscifrac2021/dependencies/mmp_atlas_grouped.mat')
neighbs = mmpgroup.neighbs;
%datasetinfo.atlas = mmpgroup;
%datasetinfo.atlasname = 'mmp';

%neighbs = ecc_getneighbs(datasetinfo);

for i = 1:length(fields)
meandata.(fields{i}).label = mmpgroup.parcellationlabel;
meanpost.(fields{i}).label = mmpgroup.parcellationlabel;
meanbl.(fields{i}).label = mmpgroup.parcellationlabel;
meanblflat.(fields{i}).label = mmpgroup.parcellationlabel;
end

% NaNify the voxels not part of the atlas
for i = 1:length(fields)
    meanpost.(fields{i}).fourierspctrm(:,[1 24],:,:) = NaN.*meanpost.(fields{i}).fourierspctrm(:,[1 24],:,:);
    meanbl.(fields{i}).fourierspctrm(:,[1 24],:,:) = NaN.*meanbl.(fields{i}).fourierspctrm(:,[1 24],:,:);
    meanblflat.(fields{i}).fourierspctrm(:,[1 24],:,:) = NaN.*meanblflat.(fields{i}).fourierspctrm(:,[1 24],:,:);
    meandata.(fields{i}).fourierspctrm(:,[1 24],:,:) = NaN.*meandata.(fields{i}).fourierspctrm(:,[1 24],:,:);
    intdata.post(:,[1 24],:) = intdata.post(:,[1 24],:).*NaN;
    intdata.bl(:,[1 24],:) = intdata.bl(:,[1 24],:).*NaN;
    pledata.post(:,[1 24],:) = pledata.post(:,[1 24],:).*NaN;
    pledata.bl(:,[1 24],:) = pledata.bl(:,[1 24],:).*NaN;
end

for c = 1:3
    cfg = []; cfg.channel = 'all';%cfg.channel = [2:23 25:46]; 
    cfg.avgoverchan = 'no'; cfg.latency = [0 0.75];
    cfg.frequency = [2 85]; cfg.method = 'montecarlo'; cfg.statistic = 'ft_statfun_fast_signrank';%cfg.statistic = 'ft_statfun_actvsblT';
    cfg.correctm = 'cluster'; cfg.clusteralpha = 0.05; cfg.clusterstatistic = 'maxsum';
    cfg.tail = 0; cfg.clustertail = 0; cfg.alpha = 0.025; cfg.numrandomization = 2000;
    cfg.parameter = 'fourierspctrm';
    
    ntrials = size(meanpost.(fields{c}).fourierspctrm,1);
    design  = zeros(2,2*ntrials);
    design(1,1:ntrials) = 1;
    design(1,ntrials+1:2*ntrials) = 2;
    design(2,1:ntrials) = [1:ntrials];
    design(2,ntrials+1:2*ntrials) = [1:ntrials];
    
    cfg.design = design;
    cfg.ivar = 1;
    cfg.uvar = 2;
    
    cfg.parpool = 48;
    
    cfg.neighbours = neighbs;
    
    %meanbl.(fields{c}).time = meanpost.(fields{c}).time;
    
    stats_signrank{c} = ft_freqstatistics(cfg,meanpost.(fields{c}),meanblflat.(fields{c}));
end

% Cluster stats for PLE and broadband

datasetinfo = struct('atlas',mmpgroup,'atlasname','mmp');
datasetinfo.neighbs = neighbs;
datasetinfo.label = mmpgroup.parcellationlabel;
%datasetinfo = settings.datasetinfo;

opts = []; opts.nrand = 2000; opts.parpool = 48; opts.minnbchan = 0;
%opts.external = cfg.design;

stats_ple = EasyClusterCorrect({permute(pledata.post,[2 1 3]) repmat(mean(pledata.bl,3)',1,1,size(pledata.post,3))},...
    datasetinfo,'ft_statfun_fast_signrank',opts);

%stats_bb = EasyClusterCorrect({permute(squeeze(mean(meanpost.frac.fourierspctrm,3)),[2 1 3]) repmat(mean(squeeze(mean(meanbl.frac.fourierspctrm,3)),3)',1,1,38)},...
%    datasetinfo,'ft_statfun_fast_signrank',opts);

stats_int = EasyClusterCorrect({permute(intdata.post,[2 1 3]) repmat(mean(intdata.bl,3)',1,1,size(pledata.post,3))},...
    datasetinfo,'ft_statfun_fast_signrank',opts);

%settings = NA_alpha_pf(settings);
settings.nfreqs = 6;
settings.tfparams.fbands = repmat({[] [2 4] [4 8] [8 13] [13 30] [30 85]},size(meandata.mixd.fourierspctrm,1),1);

freq = meandata.mixd.freq;
fbands = cell(1,7);
for i = 2:settings.nfreqs
    tmp = meandata.osci.fourierspctrm;
    %tmp = abs(tmp);
    tmp(tmp<=0) = NaN;
    for q = 1:size(meandata.osci.fourierspctrm,1)
        tmpfrange = intersect(find(meandata.osci.freq >= settings.tfparams.fbands{q,i}(1)),...
            find(meandata.osci.freq <= settings.tfparams.fbands{q,i}(2)));
        fbands{i}.post(q,:,:) = squeeze(simps(freq(tmpfrange),tmp(q,:,tmpfrange,postindx),3));
        fbands{i}.bl(q,:,:) = squeeze(simps(freq(tmpfrange),tmp(q,:,tmpfrange,blindx),3));
    end
    stats_fbands{i} = EasyClusterCorrect({permute(fbands{i}.post,[2 1 3]) repmat(mean(fbands{i}.bl,3)',1,1,size(pledata.post,3))},...
        datasetinfo,'ft_statfun_fast_signrank',opts);
end

load('/group/northoff/share/camcan/Rest/Pseudo_epochs/IRASAtf/camcan_irasatf_rest_tresoutputs.mat');
outputs.meas = outputs.meas([2 1 3:7]); outputs.data = outputs.data(:,:,[2 1 3:7],:);

%% Plotting the figure

frange = find(ismember(meandata.mixd.freq,stats{1}.freq));
for i = 1:3
    plotdata = squeeze(nanmean(10.*log10(meanpost.(fields{i}).fourierspctrm)-10.*log10(meanblflat.(fields{i}).fourierspctrm),1));
    plotdata = plotdata(:,frange,:);
    p = tf_clustplot(stats_signrank{i},mmpgroup,[-1 1],plotdata);
    save(['Fig4a_' num2str(i) '_erspplot.mat'],'p','plotdata','stats_signrank')
    export_fig(['Fig4a_' num2str(i) '_erspplot.mat'])
end

savefig('Fig4a_erspplot.fig'); export_fig('Fig4a_erspplot.png','-m5')

%part 2

settings.datatype = 'source';
settings.layout = mmpgroup;

datatemplate = zeros(size(pledata.post)); datatemplate = permute(datatemplate,[2 3 1]);

figure

p = panel('no-manage-font');

p.pack('h',{1/2 1/2})
p(1).pack('h',{1/2 1/2})
p(2).pack(2,3)
p.marginbottom = 22; p.marginleft = 22;
p.marginright = 12;
p(2).marginleft = 25;
p.margintop = 10;
set(gcf,'units','normalized','position',[0 0 1 1])

%p(1,4).pack('v',{50 50})

plotdata = effsizetc.ple;
%plotdata = 100*(pledata.post-nanmean(pledata.bl,3))./nanmean(pledata.bl,3);
%plotdata = permute(plotdata,[2 3 1]);
%pleindx = nansum(nansum(plotdata.*repmat(permute(stats_ple.mask,[3 1 2]),49,1,1),2),3);

plot_tc_topo(meanpost.mixd.time*1000,plotdata+rand(size(datatemplate))*0.001,settings,p,{1 2},'color',{'r'})
%stdshade(meanpost.mixd.time,squeeze(mean(pledata.post,2)),'k',0.15,2,'sem')
%Plot_sigmask(gca,stats_ple.mask,'cmapline','LineWidth',4)
p(1,2,1).select()
hold on
xlabel('Time (s)')
%ylabel('Percent change from baseline (%)')
ylabel('Pseudo effect size')
title('Fractal PLE')
FixAxes(gca,16)
set(gca,'XLim',[min(meanpost.mixd.time) max(meanpost.mixd.time)]*1000)
clust_sigmasks(gca,stats_ple,1)
%H = sigstar({[min(meanpost.mixd.time) max(meanpost.mixd.time)]*1000},stats_ple.negclusters.prob,0,18);
%delete(H(1));


%plotdata = 100*(10.^intdata.post-nanmean(10.^intdata.bl,3))./nanmean(10.^intdata.bl,3);
%plotdata = permute(plotdata,[2 3 1]);
plotdata = effsizetc.int;
%intindx = nansum(nansum(plotdata.*repmat(permute(stats_int.mask,[3 1 2]),49,1,1),2),3);

plot_tc_topo(meanpost.mixd.time*1000,plotdata+rand(size(datatemplate))*0.001,settings,p,{1 1},'color',{'r'})
%Plot_sigmask(gca,stats_int.mask,'cmapline','LineWidth',4)
p(1,1,1).select()
FixAxes(gca,16)
xlabel('Time (s)')
%ylabel('Percent change from prestim (%)')
title('Fractal intercept')
%ylabel('Percent change from baseline (%)')
ylabel('Pseudo effect size')
set(gca,'XLim',[min(meanpost.mixd.time) max(meanpost.mixd.time)]*1000)
clust_sigmasks(gca,stats_int,1)
%H = sigstar({[min(meanpost.mixd.time) max(meanpost.mixd.time)]*1000},stats_int.negclusters.prob,0,18);
%delete(H(1));

%yl=Normalize_Ylim([p(1,1,1).axis,p(1,2,1).axis],0);
%Set_Ylim([p(1,1,1).axis,p(1,2,1).axis],[0 yl(2)]);


% 
tmpindx = [0 1 2 3 1 2 3];
for i = 2:settings.nfreqs
   plotdata = effsizetc.(fields{i+1});
   %plotdata = 100*(fbands{i}.post-nanmean(fbands{i}.bl(:,:,end-4:end),3))./nanmean(fbands{i}.bl(:,:,end-4:end),3);
   %plotdata = permute(plotdata,[2 3 1]);
   %fbandindx(:,i) = nansum(nansum(plotdata.*repmat(permute(stats_fbands{i}.mask,[3 1 2]),49,1,1),2),3);
   plot_tc_topo(meanpost.mixd.time*1000,plotdata+rand(size(datatemplate))*0.001,settings,p,{2,ceil((i-1)/3),tmpindx(i)},'avgfunc',@nanmedian);
   p(2,ceil((i-1)/3),tmpindx(i),1).select()
   FixAxes(gca,12)
   xlabel('Time (s)')
   if i ==2 || i == 5
    %ylabel('Percent change from baseline (%)')
    ylabel('Pseudo effect size')
   end
   title(settings.tfparams.fbandnames{i})
   set(gca,'XLim',[min(meanpost.mixd.time) max(meanpost.mixd.time)]*1000)
   ax = gca;
   clust_sigmasks(gca,stats_fbands{i},1)
   %H = sigstar({[min(meanpost.mixd.time) max(meanpost.mixd.time)]*1000},stats_int.negclusters.prob,0,18);
   %delete(H(1));
end

colormap(lkcmap2)


save('Fig4b_param_tc_pseudoeffsize.mat','p')
savefig('Fig4b_param_tc_pseudoeffsize.fig'); export_fig('Fig4b_param_tc_prcchange.png','-m5')


load('Oscifrac_data_camcan_final.mat','tresdata')
ccids = cellstr(extractBefore(extractfield(files,'name'),'_'));
ccids_all = cellstr(extractBefore(tresdata.sub,'_'));

tmp = find(contains(ccids_all,ccids));
tresdata.sub = tresdata.sub(tmp); tresdata.data = tresdata.data(tmp,:,:,:);
for i = 1:length(mmpgroup.parcellationlabel)
   tresdata.grpdata(:,i,:,:) = nanmean(tresdata.data(:,mmpgroup.roiassign==i,:,:),2); 
end

fbandnames = {'delta','theta','alpha','beta','gamma'};

%effsizefun = @(post,bl,mask) nansum(nansum((nanmean(post-nanmean(bl,3),1)./nanstd(post-nanmean(bl,3),[],1)).*mean(mask,1)));
%effsizefun = @(post,bl,mask) nansum(nansum((post-nanmean(bl,3))./nanstd(bl,[],3)));
%effsizefun = @(post,bl,mask) nansum(nansum(abs(nanmedian(100*(post-nanmean(bl,3))./nanmean(bl,3),1)).*mean(mask,1),3),2);
effsizefun = @(post,bl,mask) nansum(nansum(abs(nanmean(log10(post)-nanmean(log10(bl),3),1)./((repmat(nanstd(nanmean(log10(bl),3),[],1),1,1,size(post,3))+nanstd(log10(post),[],1))/2)).*mean(mask,1),3),2);
%effsizefun = @(post,bl,mask) nansum(nansum(nanmean(abs(post-nanmean(bl,3)./nanstd(bl,[],3)),1).*mean(mask,1),3),2);
effsizefun = @(post,bl,reststd,mask) nansum(nansum(abs(nanmean((log10(post)-nanmean(log10(bl),3))./reststd,1)).*mean(mask,1),3),2);

effsizedifffun = @(effsizefun,post1,bl1,rest1,mask1,post2,bl2,rest2,mask2) effsizefun(post1,bl1,rest1,mask1)-effsizefun(post2,bl2,rest2,mask2);

effsizetc = struct;
%effsizetc.int = nanmean(10.^intdata.post-nanmean(10.^intdata.bl,3),1)./nanstd(10.^intdata.bl,[],3);

%effsizetc.int = nanmean(10.^intdata.post-nanmean(10.^intdata.bl,3),1)./...
%    ((repmat(nanstd(nanmean(10.^intdata.bl,3),[],1),1,1,size(pledata.post,3))+nanstd(10.^intdata.post,[],1))/2);
%effsizetc.ple = nanmean(10.^pledata.post-nanmean(10.^pledata.bl,3),1)./...
%    ((repmat(nanstd(nanmean(10.^pledata.bl,3),[],1),1,1,size(pledata.post,3))+nanstd(10.^pledata.post,[],1))/2);

% Cohen's d version
effsizetc.int = nanmean(intdata.post-nanmean(intdata.bl,3),1)./...
    ((repmat(nanstd(nanmean(intdata.bl,3),[],1),1,1,size(intdata.post,3))+nanstd(intdata.post,[],1))/2);
effsizetc.ple = nanmean(pledata.post-nanmean(pledata.bl,3),1)./...
    ((repmat(nanstd(nanmean(pledata.bl,3),[],1),1,1,size(pledata.post,3))+nanstd(pledata.post,[],1))/2);
for i = 2:settings.nfreqs
   %tmppost = log10(fbands{i}.post); tmppost(find(imag(tmppost)~=0)) = NaN;
   %tmpbl = log10(fbands{i}.post); tmpbl(find(imag(tmpbl)~=0)) = NaN;
   
   %outdata = outputs.data(:,:,i+1,:); outdata(outdata<=0) = NaN; outdata = log10(outdata);
   %effsizetc.(fbandnames{i-1}) = (tmppost-nanmean(tmpbl,3))./...
   %    nanstd(tmp,[],4);
    effsizetc.(fbandnames{i-1}) = nanmean(log10(fbands{i}.post) - nanmean(log10(fbands{i}.bl),3),1)./...
        ((repmat(nanstd(nanmean(log10(fbands{i}.bl),3),[],1),1,1,size(fbands{i}.post,3))+nanstd(log10(fbands{i}.post),[],1))/2);
end
effsizetc = structfun(@squeeze,effsizetc,'UniformOutput',false);

% pseudoeffectsize version
effsizetc.int = nanmean((intdata.post-nanmean(intdata.bl,3))./nanstd(tresdata.grpdata(:,:,1,:),[],4),1);
effsizetc.ple = nanmean((pledata.post-nanmean(pledata.bl,3))./nanstd(tresdata.grpdata(:,:,2,:),[],4),1);
for i = 2:settings.nfreqs
    effsizetc.(fbandnames{i-1}) = nanmean((log10(fbands{i}.post) - nanmean(log10(fbands{i}.bl),3))./...
        nanstd(tresdata.grpdata(:,:,i,:),[],4),1);
end
effsizetc = structfun(@squeeze,effsizetc,'UniformOutput',false);




allstats = [{stats_int stats_ple} stats_fbands(2:end)];

fields = fieldnames(effsizetc);

for i = 1:length(fields)
    effsizesum(i) = nansum(nansum(abs(effsizetc.(fields{i})).*allstats{i}.mask));
end

effsizediff = ones(length(fields));
for i = 1:length(fields)
   if i == 1
       post1 = 10.^intdata.post;
       bl1 = repmat(nanmean(10.^intdata.bl,3),1,1,size(intdata.post,3)); % b/c intdata already logged
   elseif i == 2
       post1 = 10.^pledata.post;
       bl1 = repmat(nanmean(10.^pledata.bl,3),1,1,size(pledata.post,3)); % b/c don't want PLE logged
   else
       post1 = fbands{i-1}.post;
       bl1 = repmat(nanmean(fbands{i-1}.bl,3),1,1,size(fbands{i-1}.post,3));
   end
   mask1 = repmat(permute(allstats{i}.mask,[3 1 2]),size(pledata.post,1),1,1);
   rest1 = repmat(nanstd(squeeze(tresdata.grpdata(:,:,i,:)),[],3),1,1,size(mask1,3));
   effsizesum(i) = effsizefun(post1,bl1,rest1,mask1);
   effsizesum_ci(i,:) = bootci(10000,{effsizefun,post1,bl1,rest1,mask1},'Options',struct('UseParallel',true));
   for ii = 1:(i-1)
       if ii == 1
           post2 = 10.^intdata.post;
           bl2 = repmat(nanmean(10.^intdata.bl,3),1,1,size(intdata.post,3));
       elseif ii == 2
           post2 = 10.^pledata.post;
           bl2 = repmat(nanmean(pledata.bl,3),1,1,size(pledata.post,3));
       else
           post2 = fbands{ii-1}.post;
           bl2 = repmat(nanmean(fbands{ii-1}.bl,3),1,1,size(fbands{ii-1}.post,3));
       end
        rest2 = repmat(nanstd(squeeze(tresdata.grpdata(:,:,ii,:)),[],3),1,1,size(mask1,3));
       mask2 = repmat(permute(allstats{ii}.mask,[3 1 2]),size(pledata.post,1),1,1);
       [~,bootstat,S] = bootci_swt(10000,{@(p1,b1,r1,m1,p2,b2,r2,m2)effsizedifffun(effsizefun,p1,b1,r1,m1,p2,b2,r2,m2),...
            post1,bl1,rest1,mask1,post2,bl2,rest2,mask2},'Options',struct('UseParallel',true));
        effsizediff(i,ii) = ibootp(0,bootstat,S);
   end
end



% indxmat = abs([intindx pleindx fbandindx(:,2:end)]);
% anovap = friedman(indxmat);
% for i = 1:size(indxmat,2)
%     for ii = 1:i
%         mcp(i,ii) = signrank(indxmat(:,i),indxmat(:,ii)); 
%     end
% end


% panel c : effect size
figure

p = panel('no-manage-font');

set(gcf,'units','normalized','position',[0 0.3 1 0.7])


p.pack('h',{60 40})
p(1).select()
b1 = bar([1 2],effsizesum(1:2),'FaceColor',palecol([1 0 0],0.3));
hold on
b2 = bar([4:(length(effsizesum)+1)],effsizesum(3:end),'FaceColor',palecol([0 0 1],0.3));
set(gca,'XTick',[1 2 4:(length(effsizesum)+1)],'XTickLabel',measnames_short)
e = errorbar([1 2 4:(length(effsizesum)+1)],effsizesum,...
    effsizesum'-effsizesum_ci(:,1),effsizesum_ci(:,2)-effsizesum',...
    'LineStyle','none','LineWidth',2,'Color','k');
%ylabel('Summed effect size')
ylabel('Summed absolute percent change')
FixAxes(gca,14)
set(gca,'XLim',[0 8.75])
%sigstar_frommat(1:size(depeffsize_sum,2),pdif_boot(:,:,i));
p(2).select()
[~,~,cbar]=imagesc_stars(vert(effsizesum)-horz(effsizesum),effsizediff,measnames_short)
FixAxes(gca,14)
cbar.Label.String = 'Summed absolute percent change difference';
cbar.Label.FontSize = 14;

p.marginleft = 26;
Normalize_Clim(gcf,1)
p.marginbottom = 20;
p.marginright = 30;

savefig('Fig4c_effsize_sum.fig'); export_fig('Fig4c_effsize_sum.png','-m5')

figure

p = panel('no-manage-font');

set(gcf,'units','normalized','position',[0.15 0.3 0.7 0.7])


p.pack()
p.pack({[0.7 0.7 0.25 0.3]})
p(1).select()
b1 = bar([1 2],mean(indxmat(:,1:2),1),'FaceColor',palecol([1 0 0],0.3));
hold on
b2 = bar([4:(length(effsizesum)+1)],mean(indxmat(:,3:end),1),'FaceColor',palecol([0 0 1],0.3));
set(gca,'XTick',[1 2 4:(length(effsizesum)+1)],'XTickLabel',measnames_short)
%e = errorbar([1 2 4:(length(effsizesum)+1)],effsizesum,...
%    effsizesum'-effsizesum_ci(:,1),effsizesum_ci(:,2)-effsizesum',...
%    'LineStyle','none','LineWidth',2,'Color','k');
ylabel('Summed percent change')
FixAxes(gca,14)
set(gca,'XLim',[0 8.75])
%sigstar_frommat(1:size(depeffsize_sum,2),pdif_boot(:,:,i));
p(2).select()
[~,~,cbar]=imagesc_stars(horz(mean(indxmat,1))-vert(mean(indxmat,1)),mcp,measnames_short)
cbar.Label.String = 'Summed percent change difference difference';
cbar.Label.FontSize = 12;

p.marginleft = 26;
Normalize_Clim(gcf,1)

savefig('Fig1c_prcchange_sum.fig'); export_fig('Fig1c_prcchange_sum.png','-m5')




% p(1,4,2).select()
% stdshade(meanpost.mixd.time,squeeze(mean(mean(meanpost.frac.fourierspctrm,2),3)),'k',0.15,2,'sem')
% Plot_sigmask(gca,stats_bb.mask,'cmapline','LineWidth',3)
% FixAxes(gca,12)
% xlabel('Time (s)')
% ylabel('Fractal broadband power')
% set(gca,'XLim',[min(meanpost.mixd.time) max(meanpost.mixd.time)])




% 
% 
% p(2).pack('h',repmat({1/length(mdl)},1,length(mdl)));
% 
% for i = 1:length(mdl)
%    p(2,i).select()
%    scatter(1:3,mdl{i}.Coefficients.Estimate(2:end),36,[0 0 1],'x','LineWidth',2)
%    hold on
%    %scatter(1:3,partr(i,:),36,[1 0 0],'o','LineWidth',2)
%    er = errorbar(1:3,mdl{i}.Coefficients.Estimate(2:end),mdl{i}.Coefficients.SE(2:end)*1.96,...
%        'LineStyle','none','LineWidth',2,'Color',[0 0 1],'HandleVisibility','off');
%    xl = get(gca,'XLim');
%    line(xl+[-0.1 0.1],[0 0],'LineWidth',1.5,'Color',[0.5 0.5 0.5],'HandleVisibility','off');
%    set(gca,'XLim',xl + [-0.1 0.1],'XTickLabel',{'Oscilatory Power','PLE','Fractal Broadband'})
%    %legend({'Regression Coefficient','Partial R^2'})
%    %ylabel('Regression Coefficient')
%    FixAxes(gca,14)
%    fix_xticklabels(gca,0.1,{'FontSize',14}) 
%    ylabel('Coefficient')
% end

% p(1).de.marginleft = 30;
% p.marginright = 20;
% p.marginleft = 18;
% p.margintop = 8;



