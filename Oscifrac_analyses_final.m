% new script for oscifrac paper

% set up some relevant variables
load lkcmap2
load('atlas_MMP1.0_4k.mat')
%load('mmp_atlas_grouped.mat'); atlas = mmpgroup;
load('Oscifrac_data_abs.mat','commonsubs')

measord = [2 1 3:7];

measnames = {'PLE','Fractal intercept','Delta power','Theta power','Alpha power','Beta power','Gamma power'};
measnames_short = {'PLE','Intercept','Delta','Theta','Alpha','Beta','Gamma'};
measnames_mid = {'PLE','Fractal intercept','Delta','Theta','Alpha','Beta','Gamma'};
measnames = measnames(measord);
measnames_short = measnames_short(measord);
measnames_mid = measnames_mid(measord);

inflated = ft_read_headshape({'/Users/Soren/Documents/MATLAB/megconnectome-3.0/template/Conte69.L.inflated.4k_fs_LR.surf.gii',...
    '/Users/Soren/Documents/MATLAB/megconnectome-3.0/template/Conte69.R.inflated.4k_fs_LR.surf.gii'});
sourcemodel = load('Conte69_fs_LR_8k.mat');

% load in the data

% load in HCP source data
cd newsource
files = dir('*fracpower.mat');
%files = dir('*IRASAmeas_abs.mat');
%files = dir('*group_irasameas.mat');
filesord = [3:5 8:9 6:7 1:2];

for i = 1:length(files)
    load(fullfile(files(filesord(i)).folder,files(filesord(i)).name))
    alloutputs{i} = outputs;
end

for i = 1:length(alloutputs)
    tmp=alloutputs{i}.data(:,:,3:end);
    tmp(find(tmp<0)) = NaN;
    alloutputs{i}.data(:,:,3:end) = log10(tmp);
    alloutputs{i}.data = alloutputs{i}.data(:,:,measord);
    alloutputs{i}.meas = alloutputs{i}.meas(measord);
end

avgrest = alloutputs{1};
avgrest.data = nanmean(cat(4,alloutputs{1}.data,alloutputs{2}.data,alloutputs{3}.data),4);

taskoutputs = alloutputs([4 6 8]); %fix this when you get the other task recordings

taskoutputs2 = alloutputs([5 7 9]);

for i = 1:3
    taskoutputs{i}.data = nanmean(cat(4,taskoutputs{i}.data,taskoutputs2{i}.data),4);
end
clear taskoutputs2

nchans = length(avgrest.chan);
badchans = [find(strcmpi(avgrest.chan,'L_???'));find(strcmpi(avgrest.chan,'R_???'))];
chansind = except(1:nchans,badchans);

% load in time-resolved HCP resting states
files = dir('*_timeres.mat');
load(files(1).name)
tresdata = outputs;
for i = 1:3
    load(files(i).name)
    tresdata.data = cat(4,tresdata.data,outputs.data);
end
 
tmp = tresdata.data(:,:,3:7,:);
tmp(find(tmp < 0)) = NaN;
tresdata.data(:,:,3:7,:) = log10(tmp);
tresdata.data = tresdata.data(:,:,measord,:);
tresdata.meas = tresdata.meas(measord);
cd ..

% alternatively: load in camcan data
% load('rest_continuous_irasameas_fbands_abs.mat'); 
% avgrest = outputs;
% load('task_continuous_irasameas_fbands_abs.mat'); 
% taskoutputs{1} = outputs;
% subids_rest = extractBefore(avgrest.sub,'_MMP');
% subids_task = extractBefore(taskoutputs{1}.sub,'_MMP');
% [i1,i2] = match_str(subids_rest,subids_task);
% avgrest.data = avgrest.data(i1,:,:); avgrest.sub = avgrest.sub(i1);
% taskoutputs{1}.data = taskoutputs{1}.data(i2,:,:); taskoutputs{1}.sub = taskoutputs{1}.sub(i2);
% commonsubs{1,4} = i1; commonsubs{4,1} = i2; 
% 
% tmp=avgrest.data(:,:,3:end);
% tmp(find(tmp<0)) = NaN;
% avgrest.data(:,:,3:end) = log10(tmp);
% avgrest.data = avgrest.data(:,:,measord);
% avgrest.meas = avgrest.meas(measord);
% 
% for i = 1:length(taskoutputs)
%     tmp=taskoutputs{i}.data(:,:,3:end);
%     tmp(find(tmp<0)) = NaN;
%     taskoutputs{i}.data(:,:,3:end) = log10(tmp);
%     taskoutputs{i}.data = taskoutputs{i}.data(:,:,measord);
%     taskoutputs{i}.meas = taskoutputs{i}.meas(measord);
% end

% get statistics for rest-task differences, effect sizes
cfg = []; cfg.multcompare = 'cluster';
cfg.multcompare = 'fdr';
cfg.cluster.datasetinfo.atlas = atlas; cfg.cluster.datasetinfo.atlasname = 'mmp';
cfg.cluster.datasetinfo.label = avgrest.chan;
cfg.cluster.nrand = 10000;
cfg.test = 'signrank'; cfg.effectsize = 'hedgesg';
cfg.channel = chansind;

for i = 1:length(taskoutputs)
    cfg.subs = [commonsubs(i+3,1) commonsubs(1,i+3)];
    stats_resttask(i,:) = ft_measurestatistics(cfg,{taskoutputs{i} avgrest});
end

% get percent change values for rest-task differences
clear prcchange_sub

for i = 1:length(taskoutputs)
    for ii = 1:length(avgrest.meas)
        if ii == 2
            prcchange_sub{i}(:,:,ii) = 100*(taskoutputs{i}.data(commonsubs{i+3,1},:,ii)-...
                avgrest.data(commonsubs{1,i+3},:,ii))./abs(avgrest.data(commonsubs{1,i+3},:,ii));
        else
             prcchange_sub{i}(:,:,ii) = 100*(10.^taskoutputs{i}.data(commonsubs{i+3,1},:,ii)...
                 -10.^avgrest.data(commonsubs{1,i+3},:,ii))./(10.^avgrest.data(commonsubs{1,i+3},:,ii));
        end
    end
    prcchange_sub{i}(:,badchans,:) = [];
end
% for i = 1:7
%     for ii = 1:length(taskoutputs)
%         prcchange_sum(ii,i) = sum(abs(prcchange{ii}(:,i)).*vert(stats_resttask{ii,i}.mask));
%     end
% end

clear pseudoeffsize
for i = 1:length(taskoutputs)
    for ii = 1:length(avgrest.meas)
        r = avgrest.data(commonsubs{1,i+3},chansind,ii); 
        t = taskoutputs{i}.data(commonsubs{i+3,1},chansind,ii);
        tresstd = squeeze(nanstd(tresdata.data(commonsubs{1,i+3},chansind,ii,:),[],4));
        pseudoeffsize_sub{i}(:,:,ii) = (t-r)./tresstd;
    end
end

for i = 1:length(taskoutputs)
    for ii = 1:length(avgrest.meas)
        r = avgrest.data(commonsubs{1,i+3},chansind,ii); 
        t = taskoutputs{i}.data(commonsubs{i+3,1},chansind,ii);
        cohend{i}(1,:,ii) = nanmean(t-r,1)./((nanstd(t,[],1)+nanstd(r,[],1))/2); % keep same dimensionality as other measures
    end
end

% analysis: rest-task diff CI's and differences with pseudo effect size
effsizefunc = @(r,t,tresstd) nansum(abs((nanmean((t-r)./tresstd,1)).*double(fdr(signrank_mat(t,r,1))<0.05)));
%effsizefunc = @(r,t,m,tresstd) nansum(abs((nanmean((t-r)./tresstd,1)).*nanmean(m,1)));

pdif_pseudoeff = cell(1,3);
pseudoeffsize = cell(1,3);
pseudoeffsize_ci = cell(1,3);
%pdif_pseudoeff = ones(length(avgrest.meas),length(avgrest.meas),3);

for i = 1:length(taskoutputs)
    pdif_pseudoeff{i} = ones(length(avgrest.meas),length(avgrest.meas));
    for q = 1:length(avgrest.meas)
        r = avgrest.data(commonsubs{1,i+3},chansind,q);
        t = taskoutputs{i}.data(commonsubs{i+3,1},chansind,q);
        m = repmat(horz(stats_resttask{i,q}.mask),size(r,1),1);
        tresstd = squeeze(nanstd(tresdata.data(commonsubs{1,i+3},chansind,q,:),[],4));
        
        pseudoeffsize{i}(1,q) = effsizefunc(r,t,tresstd);
        pseudoeffsize_ci{i}(1,q,:) = bootci(10000,effsizefunc,r,t,tresstd);
        for qq = 1:(q-1)
            r1 = avgrest.data(commonsubs{1,i+3},chansind,q);
            r2 = avgrest.data(commonsubs{1,i+3},chansind,qq);
            t1 = taskoutputs{i}.data(commonsubs{i+3,1},chansind,q);
            t2 = taskoutputs{i}.data(commonsubs{i+3,1},chansind,qq);
            m1 = repmat(horz(stats_resttask{i,q}.mask),size(r1,1),1);
            m2 = repmat(horz(stats_resttask{i,qq}.mask),size(r1,1),1);
            tresstd1 = squeeze(nanstd(tresdata.data(commonsubs{1,i+3},chansind,q,:),[],4));
            tresstd2 = squeeze(nanstd(tresdata.data(commonsubs{1,i+3},chansind,qq,:),[],4));

            [~,bootstat,S] = bootci_swt(10000,{@(r1,t1,tresstd1,r2,t2,tresstd2) effsizefunc(r1,t1,tresstd1)-effsizefunc(r2,t2,tresstd2),...
                r1,t1,tresstd1,r2,t2,tresstd2},'Options',struct('UseParallel',true));
            pdif_pseudoeff{i}(q,qq) = ibootp(0,bootstat,S);
        end
    end
    pdif_pseudoeff{i}(:,:) = bdiagtomat(bonf_holm(belowDiag(pdif_pseudoeff{i}(:,:)),0.05),7);
end
pdif_pseudoeff = cat(3,pdif_pseudoeff{:});
pseudoeffsize = cat(1,pseudoeffsize{:});
pseudoeffsize_ci = cat(1,pseudoeffsize_ci{:});


for q = 1:length(taskoutputs)
    for i = 1:length(chansind)
        effint = (taskoutputs{q}.data(commonsubs{q+3,1},chansind(i),1)-avgrest.data(commonsubs{1,q+3},chansind(i),1))./...
            squeeze(nanstd(tresdata.data(commonsubs{1,q+3},chansind(i),1,:),[],4));
        effple = (taskoutputs{q}.data(commonsubs{q+3,1},chansind(i),2)-avgrest.data(commonsubs{1,q+3},chansind(i),2))./...
            squeeze(nanstd(tresdata.data(commonsubs{1,q+3},chansind(i),2,:),[],4));
        p_intvsple(q,i) = signrank(effint,effple);
    end
    p_intvsple(q,:) = fdr(p_intvsple(q,:));
end


% analysis: rest-task diff CI's and differences with Cohen's d
effsizefunc = @(r,t) nansum(abs(nanmean(t-r,1)./((nanstd(t,[],1)+nanstd(r,[],1))/2).*double(fdr(signrank_mat(t,r,1))<0.05)));
%effsizefunc = @(r,t,m) nansum(abs(nanmean(t-r,1)./((nanstd(t,[],1)+nanstd(r,[],1))/2).*nanmean(m,1)));

pdif_cohend = cell(1,3);
cohend_sum = cell(1,3);
cohend_ci = cell(1,3);
%pdif_cohend = ones(length(avgrest.meas),length(avgrest.meas),3);

for i = 1:length(taskoutputs)
    pdif_cohend{i} = ones(length(avgrest.meas),length(avgrest.meas));
    for q = 1:length(avgrest.meas)
        r = avgrest.data(commonsubs{1,i+3},chansind,q);
        t = taskoutputs{i}.data(commonsubs{i+3,1},chansind,q);
        m = repmat(horz(stats_resttask{i,q}.mask),size(r,1),1);
        tresstd = squeeze(nanstd(tresdata.data(commonsubs{1,i+3},chansind,q,:),[],4));
        
        cohend_sum{i}(1,q) = effsizefunc(r,t);
        cohend_ci{i}(1,q,:) = bootci(10000,effsizefunc,r,t);
        for qq = 1:(q-1)
            r1 = avgrest.data(commonsubs{1,i+3},chansind,q);
            r2 = avgrest.data(commonsubs{1,i+3},chansind,qq);
            t1 = taskoutputs{i}.data(commonsubs{i+3,1},chansind,q);
            t2 = taskoutputs{i}.data(commonsubs{i+3,1},chansind,qq);
            m1 = repmat(horz(stats_resttask{i,q}.mask),size(r1,1),1);
            m2 = repmat(horz(stats_resttask{i,qq}.mask),size(r1,1),1);
            tresstd1 = squeeze(nanstd(tresdata.data(commonsubs{1,i+3},chansind,q,:),[],4));
            tresstd2 = squeeze(nanstd(tresdata.data(commonsubs{1,i+3},chansind,qq,:),[],4));

            [~,bootstat,S] = bootci_swt(10000,{@(r1,t1,r2,t2) effsizefunc(r1,t1)-effsizefunc(r2,t2),...
                r1,t1,r2,t2},'Options',struct('UseParallel',true));
            pdif_cohend{i}(q,qq) = ibootp(0,bootstat,S);
        end
    end
    pdif_cohend{i}(:,:) = bdiagtomat(bonf_holm(belowDiag(pdif_cohend{i}(:,:)),0.05),7);
end
pdif_cohend = cat(3,pdif_cohend{:});
cohend_sum = cat(1,cohend_sum{:});
cohend_ci = cat(1,cohend_ci{:});


% analysis: rest-task diff CI's and differences with percent change

prcchangefunc = @(r,t) nansum(abs(nanmedian((100*(t-r)./r),1)).*double(fdr(signrank_mat(t,r,1))<0.05));
%prcchangefunc = @(r,t,m) nansum(abs(nanmedian((100*(t-r)./r),1)).*nanmean(m,1));

pdif_prcchange = cell(1,3);
prcchange = cell(1,3);
prcchange_ci = cell(1,3);%pdif_prcchange = ones(length(avgrest.meas),length(avgrest.meas),3);

for i = 1:length(taskoutputs)
    pdif_prcchange{i} = ones(length(avgrest.meas),length(avgrest.meas));
    for q = 1:length(avgrest.meas)
        r = avgrest.data(commonsubs{1,i+3},chansind,q);
        t = taskoutputs{i}.data(commonsubs{i+3,1},chansind,q);
        if q ~= 2
            r = 10.^r; t = 10.^t;
        end
        m = repmat(horz(stats_resttask{i,q}.mask),size(r,1),1);
        
        prcchange{i}(1,q) = prcchangefunc(r,t);
        
        prcchange_ci{i}(1,q,:) = bootci(10000,{prcchangefunc,r,t},'Options',struct('UseParallel',true));
        for qq = 1:(q-1)
            r1 = avgrest.data(commonsubs{1,i+3},chansind,q);
            r2 = avgrest.data(commonsubs{1,i+3},chansind,qq);
            t1 = taskoutputs{i}.data(commonsubs{i+3,1},chansind,q);
            t2 = taskoutputs{i}.data(commonsubs{i+3,1},chansind,qq);
            if q ~= 2
                r1 = 10.^r1; t1 = 10.^t1; 
            end
            if qq ~= 2 
               r2 = 10.^r2; t2 = 10.^t2; 
            end
            m1 = repmat(horz(stats_resttask{i,q}.mask),size(r1,1),1);
            m2 = repmat(horz(stats_resttask{i,qq}.mask),size(r1,1),1);


            [~,bootstat,S] = bootci_swt(10000,{@(r1,t1,r2,t2) prcchangefunc(r1,t1)-prcchangefunc(r2,t2),...
                r1,t1,r2,t2},'Options',struct('UseParallel',true));
            pdif_prcchange{i}(q,qq) = ibootp(0,bootstat,S);
        end
    end
    pdif_prcchange{i}(:,:) = bdiagtomat(bonf_holm(belowDiag(pdif_prcchange{i}(:,:)),0.05),7);
end
pdif_prcchange = cat(3,pdif_prcchange{:});
prcchange = cat(1,prcchange{:});
prcchange_ci = cat(1,prcchange_ci{:});

% analysis: comparing intercept and PLE based on pseudo effect size

for q = 1:length(taskoutputs)
    for i = 1:length(chansind)
        effint = (taskoutputs{q}.data(commonsubs{q+3,1},chansind(i),1)-avgrest.data(commonsubs{1,q+3},chansind(i),1))./...
            squeeze(nanstd(tresdata.data(commonsubs{1,q+3},chansind(i),1,:),[],4));
        effple = (taskoutputs{q}.data(commonsubs{q+3,1},chansind(i),2)-avgrest.data(commonsubs{1,q+3},chansind(i),2))./...
            squeeze(nanstd(tresdata.data(commonsubs{1,q+3},chansind(i),2,:),[],4));
        peff_intvsple(q,i) = signrank(effint,effple);
    end
    peff_intvsple(q,:) = fdr(peff_intvsple(q,:));
end


% analysis: comparing intercept and PLE based on percent change

for q = 1:length(taskoutputs)
    for i = 1:length(chansind)
        prcint = 100*(10.^taskoutputs{q}.data(commonsubs{q+3,1},chansind(i),1)-10.^avgrest.data(commonsubs{1,q+3},chansind(i),1))./...
            squeeze(10.^avgrest.data(commonsubs{1,q+3},chansind(i),1));
        prcple = 100*(taskoutputs{q}.data(commonsubs{q+3,1},chansind(i),2)-avgrest.data(commonsubs{1,q+3},chansind(i),2))./...
            squeeze(avgrest.data(commonsubs{1,q+3},chansind(i),2));
        pprc_intvsple(q,i) = signrank(prcint,prcple);
    end
    pprc_intvsple(q,:) = fdr(pprc_intvsple(q,:));
end

% CV calculation and CV differences in HCP data

for i = 1:size(avgrest.data,3)
    cvinter(i) = squeeze(geocv(nanmean(avgrest.data(:,chansind,i),2),1));
    cvintra(i) = squeeze(nanmedian(geocv(nanmean(tresdata.data(:,chansind,i,:),2),4),1));
end
cvrat = cvinter./cvintra;

% Inter-subject CV
cvinter_pdiff = ones(length(cvinter));
cvdifffunc = @(r1,r2) geocv(r1)-geocv(r2);
for i = 1:length(cvinter)
    rdat = squeeze(nanmean(avgrest.data(:,chansind,i),2));
    cvinter_ci(i,:) = bootci(10000,@geocv,rdat);
    for ii = 1:(i-1)
        rdat2 = squeeze(nanmean(avgrest.data(:,chansind,ii),2));
        [~,bootstat,S] = bootci_swt(10000,{cvdifffunc,rdat,rdat2},'Options',struct('UseParallel',true));
        cvinter_pdiff(i,ii) = ibootp(0,bootstat,S);
    end
end
cvinter_pdiff = bdiagtomat(bonf_holm(belowDiag(cvinter_pdiff),0.05),length(cvinter));

% Intra-subject CV
cvintrafunc = @(rdat) nanmedian(geocv(rdat,2),1);
cvintradifffunc = @(rdat1,rdat2) cvintrafunc(rdat1)-cvintrafunc(rdat2);
cvintra_pdiff = ones(length(cvrat));
for i = 1:length(cvinter)
    rdat = squeeze(nanmean(tresdata.data(:,chansind,i,:),2));
    cvintra_ci(i,:) = bootci(10000,cvintrafunc,rdat);
    for ii = 1:(i-1)
        rdat2 = squeeze(nanmean(tresdata.data(:,chansind,ii,:),2));
        [~,bootstat,S] = bootci_swt(10000,{cvintradifffunc,rdat,rdat2},'Options',struct('UseParallel',true));
        cvintra_pdiff(i,ii) = ibootp(0,bootstat,S);
    end
end
cvintra_pdiff = bdiagtomat(bonf_holm(belowDiag(cvintra_pdiff),0.05),length(cvrat));

%cv ratio
cvratfunc = @(rrep,tresr) squeeze(geocv(nanmean(rrep,2),1))./squeeze(nanmedian(geocv(tresr,2),1));
cvratdifffunc = @(rrep1,tresr1,rrep2,tresr2) cvratfunc(rrep1,tresr1)-cvratfunc(rrep2,tresr2);

cvrat_pdiff = ones(length(cvrat));
for i = 1:length(cvrat)
    rrep = repmat(squeeze(nanmean(avgrest.data(:,chansind,i),2)),1,size(tresdata.data,4));
    tresr = squeeze(nanmean(tresdata.data(:,chansind,i,:),2));
    cvrat_ci(i,:) = bootci(10000,cvratfunc,rrep,tresr);
    for ii = 1:(i-1)
        rrep2 = repmat(squeeze(nanmean(avgrest.data(:,chansind,ii),2)),1,size(tresdata.data,4));
        tresr2 = squeeze(nanmean(tresdata.data(:,chansind,ii,:),2));
        [~,bootstat,S] = bootci_swt(10000,{cvratdifffunc,rrep,tresr,rrep2,tresr2},'Options',struct('UseParallel',true));
        cvrat_pdiff(i,ii) = ibootp(0,bootstat,S);
    end
end
cvrat_pdiff = bdiagtomat(bonf_holm(belowDiag(cvrat_pdiff),0.05),length(cvrat));

% analysis: genetics

twintbl = readtable('RESTRICTED_andreascalabrini_8_16_2018_17_37_5.csv');

subids = cellfun(@str2num,extractBefore(avgrest.sub,'_roidata'),'UniformOutput',true);
[~,subindx] = intersect(twintbl.Subject,subids);
twintbl = twintbl(subindx,:);
uniquefam = unique(twintbl.Family_ID);

mz = [];
dz = [];
nt = [];
for i = 1:length(uniquefam)
    if sum(strcmpi(uniquefam{i},twintbl.Family_ID)) == 2
        indx = find(strcmpi(uniquefam{i},twintbl.Family_ID));
        if any(strcmpi(twintbl.ZygosityGT(indx),'MZ')) || any(strcmpi(twintbl.ZygositySR(indx),'MZ'))
            mz = cat(1,mz,horz(indx));
        elseif any(strcmpi(twintbl.ZygosityGT(indx),'DZ')) || any(strcmpi(twintbl.ZygositySR(indx),'NotMZ'))
            dz = cat(1,dz,horz(indx));
        end
    else
        nt = cat(1,nt,find(strcmpi(uniquefam{i},twintbl.Family_ID)));
    end
end

ntpairs = reshape(nt(randperm(length(nt)-1)),floor(0.5*length(nt)),2);

for q = 1:size(avgrest.data,3)
    mztbl = array2table([nanmean(avgrest.data(mz(:,1),:,q),2) nanmean(avgrest.data(mz(:,2),:,q),2)],...
        'VariableNames',{[measnames_short{q} '_T1'], [measnames_short{q} '_T2']});
    dztbl = array2table([nanmean(avgrest.data(dz(:,1),:,q),2) nanmean(avgrest.data(dz(:,2),:,q),2)],...
        'VariableNames',{[measnames_short{q} '_T1'], [measnames_short{q} '_T2']});
    [ace(q),smry{q},compare{q},reduc(q),reducsmry{q}] = heritability_ACE(mztbl,dztbl);
    herit(q) = ace(q).a^2;
    herit_ci(q,:) = ace(q).ci{2,:}.*abs(ace(q).ci{2,:});
end

heritdiff = ones(size(avgrest.data,3));
heritstat = ones(size(avgrest.data,3));
for i = 1:size(avgrest.data,3)
   for ii = 1:(i-1)
       heritstat(i,ii) = (ace(i).a-ace(ii).a)/sqrt(((ace(i).a-ace(i).ci.l_ci(2))/1.96).^2 +...
           ((ace(ii).a-ace(ii).ci.l_ci(2))/1.96).^2);
       if heritstat(i,ii) > 0
          heritdiff(i,ii) = 2*(1-normcdf(heritstat(i,ii)));
       else
           heritdiff(i,ii) = 2*normcdf(heritstat(i,ii));
       end
   end
end


% analysis: behaviour

behav = readtable('unrestricted_SorenWT_1_24_2020_14_6_21.csv');

subids = cellfun(@str2num,extractBefore(avgrest.sub,'_roidata'),'UniformOutput',true);
[~,subindx] = intersect(behav.Subject,subids);
behav = behav(subindx,:);

% first try with actual values

for i = 1:7
    X = avgrest.data(:,chansind,i);
    X(:,find(~any(~isnan(X),1))) = []; X(:,find(isnan(nanstd(X,[],1)))) = [];
    [vals_pcaweights{i},comps,~,~,explained] = pca(nanzscore(X,[],1),'algorithm','als');
%     Xshuf = X;
%     for q = 1:1000
%         for qq = 1:size(X,2)
%             Xshuf(:,qq) = X(randperm(size(X,1)),qq);
%         end
%         [~,~,~,~,expl_perm(:,q)] = pca(nanzscore(Xshuf,[],1),'rows','pairwise');
%     end
    
    %ncomps = find(1-mean(explained>expl_perm,2)>0.05,1)-1;
    ncomps = find(explained<100/89)-1; % kaiser criterion
    Xcomps = comps(:,1:ncomps);
    Xcomps = [Xcomps Xcomps.^2]; % quadratic model
    [lassoweights,lasstats] = lasso(Xcomps,behav.CardSort_Unadj,'CV',5,'NumLambda',10,'RelTol',1e-3);
    [~,lasstats.IndexMinMSE] = min(lasstats.MSE(any(lassoweights,1)));
    lassoweights = lassoweights(:,lasstats.IndexMinMSE);
    vals_alllassoweights{i} = lassoweights;
    vals_sourcespaceweights{i} = ([vals_pcaweights{i}(:,1:ncomps) vals_pcaweights{i}(:,1:ncomps)])*lassoweights;
    %vals_sourcespaceweights{i} = vals_pcaweights{i}(:,1:ncomps)*lassoweights;
    pred = Xcomps*lassoweights; 
    vals_r_obs(i) = corr(pred,behav.CardSort_Unadj,'rows','pairwise');
    Xcompsshuf = Xcomps;
    for q = 1:1000
        %for qq = 1:size(Xcomps,2)
        %    Xcompsshuf(:,qq) = Xcomps(randperm(size(X,1)),qq);
        %end
        randindx = randperm(size(X,1));
        [lassoweightsshuf,lasstats] = lasso(Xcompsshuf,behav.CardSort_Unadj(randindx),'CV',5,'NumLambda',10,'RelTol',1e-3);
        [~,lasstats.IndexMinMSE] = min(lasstats.MSE(any(lassoweightsshuf,1)));
        lassoweightsshuf = lassoweightsshuf(:,lasstats.IndexMinMSE);
        pred = Xcompsshuf*lassoweightsshuf;
        vals_r_shuf(q,i) = corr(pred,behav.CardSort_Unadj(randindx),'rows','pairwise');
    end
   pbehav_vals_meas(i) = 1-nanmean(vals_r_obs(i)>vals_r_shuf(:,i));
  %  [~,~,~,~,~,~,~,stats_pls_vals(i)] = plsregress_perm(Xcomps,zscore(behav.CardSort_Unadj),1,10000);
end

% now with intra-subject CV

for i = 1:7
    X = geocv(tresdata.data(:,chansind,i,:),4);
    X(:,find(~any(~isnan(X),1))) = []; X(:,find(isnan(nanstd(X,[],1)))) = [];
    [cv_pcaweights{i},comps,~,~,explained] = pca(nanzscore(X,[],1),'algorithm','als');
%     Xshuf = X;
%     for q = 1:1000
%         for qq = 1:size(X,2)
%             Xshuf(:,qq) = X(randperm(size(X,1)),qq);
%         end
%         [~,~,~,~,expl_perm(:,q)] = pca(nanzscore(Xshuf,[],1),'rows','pairwise');
%     end
    %ncomps = find(1-mean(explained>expl_perm,2)>0.05,1)-1;
    ncomps = find(explained<100/89)-1; % kaiser criterion
    Xcomps = comps(:,1:ncomps);
    Xcomps = [Xcomps Xcomps.^2]; % quadratic model
    [lassoweights,lasstats] = lasso(Xcomps,behav.CardSort_Unadj,'CV',5,'NumLambda',10,'RelTol',1e-3);
    [~,lasstats.IndexMinMSE] = min(lasstats.MSE(any(lassoweights,1)));
    lassoweights = lassoweights(:,lasstats.IndexMinMSE);
    cv_alllassoweights{i} = lassoweights;
    cv_sourcespaceweights{i} = ([cv_pcaweights{i}(:,1:ncomps) cv_pcaweights{i}(:,1:ncomps)])*lassoweights;
    %cv_sourcespaceweights{i} = cv_pcaweights{i}(:,1:ncomps)*lassoweights;
    pred = Xcomps*lassoweights; 
    cv_r_obs(i) = corr(pred,behav.CardSort_Unadj,'rows','pairwise');
    Xcompsshuf = Xcomps;
    for q = 1:1000
%         for qq = 1:size(Xcomps,2)
%             Xcompsshuf(:,qq) = Xcomps(randperm(size(X,1)),qq);
%         end
        randindx = randperm(size(X,1));
        [lassoweightsshuf,lasstats] = lasso(Xcompsshuf,behav.CardSort_Unadj(randindx),'CV',5,'NumLambda',10,'RelTol',1e-3);
        [~,lasstats.IndexMinMSE] = min(lasstats.MSE(any(lassoweightsshuf,1)));
        lassoweightsshuf = lassoweightsshuf(:,lasstats.IndexMinMSE);
        pred = Xcompsshuf*lassoweightsshuf;
        cv_r_shuf(q,i) = corr(pred,behav.CardSort_Unadj(randindx),'rows','pairwise');
    end
    pbehav_cv_meas(i) = 1-nanmean(cv_r_obs(i)>cv_r_shuf(:,i));
    [~,~,~,~,~,~,~,stats_pls_cv(i)] = plsregress_perm(Xcomps,zscore(behav.CardSort_Unadj),1,10000);
end

% now with average rest-task percent change 

alltasksubs = intersect(commonsubs{1,4},commonsubs{1,5});
alltasksubs = intersect(alltasksubs,commonsubs{1,6}); % this gets the subject numbers relative to the resting state

for i = 1:3
   taskcomsubs{i} = find(ismember(commonsubs{1,i+3},alltasksubs)); 
end

for i = 1:7
    X = (abs(prcchange_sub{1}(taskcomsubs{1},:,i))+abs(prcchange_sub{2}(taskcomsubs{2},:,i))+abs(prcchange_sub{3}(taskcomsubs{3},:,i)))/3; % mean of absolute task changes
    X(:,find(~any(~isnan(X),1))) = []; X(:,find(isnan(nanstd(X,[],1)))) = [];
    [rt_pcaweights{i},comps,~,~,explained] = pca(nanzscore(X,[],1),'algorithm','als');
%     Xshuf = X;
%     for q = 1:1000
%         for qq = 1:size(X,2)
%             Xshuf(:,qq) = X(randperm(size(X,1)),qq);
%         end
%         [~,~,~,~,expl_perm(:,q)] = pca(nanzscore(Xshuf,[],1),'rows','pairwise');
%     end
    %ncomps = find(1-mean(explained>expl_perm,2)>0.05,1)-1;
    ncomps = find(explained<100/48)-1; % kaiser criterion
    Xcomps = comps(:,1:ncomps);
    Xcomps = [Xcomps Xcomps.^2]; % quadratic model
    [lassoweights,lasstats] = lasso(Xcomps,behav.CardSort_Unadj(alltasksubs),'CV',5,'NumLambda',10,'RelTol',1e-3);
    [~,lasstats.IndexMinMSE] = min(lasstats.MSE(any(lassoweights,1)));
    lassoweights = lassoweights(:,lasstats.IndexMinMSE);
    rt_alllassoweights{i} = lassoweights;
    rt_sourcespaceweights{i} = ([rt_pcaweights{i}(:,1:ncomps) rt_pcaweights{i}(:,1:ncomps)])*lassoweights;
    %rt_sourcespaceweights{i} = rt_pcaweights{i}(:,1:ncomps)*lassoweights;
    pred = Xcomps*lassoweights; 
    rt_r_obs(i) = corr(pred,behav.CardSort_Unadj(alltasksubs),'rows','pairwise');
    Xcompsshuf = Xcomps;
    for q = 1:1000
%         for qq = 1:size(Xcomps,2)
%             Xcompsshuf(:,qq) = Xcomps(randperm(size(X,1)),qq);
%         end
        randindx = randperm(size(X,1));
        [lassoweightsshuf,lasstats] = lasso(Xcompsshuf,behav.CardSort_Unadj(alltasksubs(randindx)),'CV',5,'NumLambda',10,'RelTol',1e-3);
        [~,lasstats.IndexMinMSE] = min(lasstats.MSE(any(lassoweightsshuf,1)));
        lassoweightsshuf = lassoweightsshuf(:,lasstats.IndexMinMSE);
        pred = Xcompsshuf*lassoweightsshuf;
        rt_r_shuf(q,i) = corr(pred,behav.CardSort_Unadj(alltasksubs(randindx)),'rows','pairwise');
    end
    pbehav_rt_meas(i) = 1-nanmean(rt_r_obs(i)>rt_r_shuf(:,i));
%    [~,~,~,~,~,~,~,stats_pls_rt(i)] = plsregress_perm(Xcomps,zscore(behav.CardSort_Unadj(alltasksubs)),1,10000);

end


for i = 1:7
    X = (abs(pseudoeffsize_sub{1}(taskcomsubs{1},:,i))+abs(pseudoeffsize_sub{2}(taskcomsubs{2},:,i))+abs(pseudoeffsize_sub{3}(taskcomsubs{3},:,i)))/3; % mean of absolute task changes
    X(:,find(~any(~isnan(X),1))) = []; X(:,find(isnan(nanstd(X,[],1)))) = [];
    [rt2_pcaweights{i},comps,~,~,explained] = pca(nanzscore(X,[],1),'algorithm','als');
%     Xshuf = X;
%     for q = 1:1000
%         for qq = 1:size(X,2)
%             Xshuf(:,qq) = X(randperm(size(X,1)),qq);
%         end
%         [~,~,~,~,expl_perm(:,q)] = pca(nanzscore(Xshuf,[],1),'rows','pairwise');
%     end
    %ncomps = find(1-mean(explained>expl_perm,2)>0.05,1)-1;
    ncomps = find(explained<100/48)-1; % kaiser criterion
    Xcomps = comps(:,1:ncomps);
    Xcomps = [Xcomps Xcomps.^2]; % quadratic model
    [lassoweights,lasstats] = lasso(Xcomps,behav.CardSort_Unadj(alltasksubs),'CV',5,'NumLambda',10,'RelTol',1e-3);
    [~,lasstats.IndexMinMSE] = min(lasstats.MSE(any(lassoweights,1)));
    lassoweights = lassoweights(:,lasstats.IndexMinMSE);
    rt2_alllassoweights{i} = lassoweights;
    rt2_sourcespaceweights{i} = ([rt2_pcaweights{i}(:,1:ncomps) rt2_pcaweights{i}(:,1:ncomps)])*lassoweights;
    %rt2_sourcespaceweights{i} = rt2_pcaweights{i}(:,1:ncomps)*lassoweights;
    pred = Xcomps*lassoweights; 
    rt2_r_obs(i) = corr(pred,behav.CardSort_Unadj(alltasksubs),'rows','pairwise');
    Xcompsshuf = Xcomps;
    for q = 1:1000
%         for qq = 1:size(Xcomps,2)
%             Xcompsshuf(:,qq) = Xcomps(randperm(size(X,1)),qq);
%         end
        randindx = randperm(size(X,1));
        [lassoweightsshuf,lasstats] = lasso(Xcompsshuf,behav.CardSort_Unadj(alltasksubs(randindx)),'CV',5,'NumLambda',10,'RelTol',1e-3);
        [~,lasstats.IndexMinMSE] = min(lasstats.MSE(any(lassoweightsshuf,1)));
        lassoweightsshuf = lassoweightsshuf(:,lasstats.IndexMinMSE);
        pred = Xcompsshuf*lassoweightsshuf;
        rt2_r_shuf(q,i) = corr(pred,behav.CardSort_Unadj(alltasksubs(randindx)),'rows','pairwise');
    end
    pbehav_rt2_meas(i) = 1-nanmean(rt2_r_obs(i)>rt2_r_shuf(:,i));
%    [~,~,~,~,~,~,~,stats_pls_rt2(i)] = plsregress_perm(Xcomps,zscore(behav.CardSort_Unadj(alltasksubs)),1,10000);

end

% just correlate medians

for i = 1:7
    X = avgrest.data(:,chansind,i);
    [meanvals_r(i),meanvals_p(i)] = corr(nanmedian(X,2),behav.CardSort_Unadj,'Type','Spearman');
    X = geocv(tresdata.data(:,chansind,i,:),4);
    [meancv_r(i),meancv_p(i)] = corr(nanmedian(X,2),behav.CardSort_Unadj,'Type','Spearman');
    X = (abs(prcchange_sub{1}(taskcomsubs{1},:,i))+abs(prcchange_sub{2}(taskcomsubs{2},:,i))+abs(prcchange_sub{3}(taskcomsubs{3},:,i)))/3; % mean of absolute task changes
    [meanrt_r(i),meanrt_p(i)] = corr(nanmedian(X,2),behav.CardSort_Unadj(alltasksubs),'Type','Spearman');
    X = (abs(pseudoeffsize_sub{1}(taskcomsubs{1},:,i))+abs(pseudoeffsize_sub{2}(taskcomsubs{2},:,i))+abs(pseudoeffsize_sub{3}(taskcomsubs{3},:,i)))/3; % mean of absolute task changes
    [meanrt2_r(i),meanrt2_p(i)] = corr(nanmedian(X,2),behav.CardSort_Unadj(alltasksubs),'Type','Spearman');
end