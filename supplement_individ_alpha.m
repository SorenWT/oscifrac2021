% supplementary analysis looking at fractal power rather than intercept

tasks = {'rest-3','rest-4','rest-5','wrkmem-6','wrkmem-7','storym-8','storym-9','motort-10','motort-11'};
basedir = '/group/northoff/share/hcp-meg/source';
cd(basedir)

cfg.measure = {};
cfg.measure{1} = @(EEG)Thetapower_individ_EEG_wrapper(EEG,4,'no','yes');
cfg.measure{2} = @(EEG)Alphapower_individ_EEG_wrapper(EEG,'no','yes');
cfg.measure{3} = @(EEG)Thetapower_individ_EEG_wrapper(EEG,30,'no','yes');
%cfg.measure{2} = @(EEG)IRASAPower_EEG_wrapper(EEG,'frac',[2 85],0);

for i = 1:length(tasks)
    cfg.dir = fullfile(basedir,tasks{i});
    cfg.outfile = fullfile(basedir,tasks{i},[tasks{i} '_irasameas_individalpha.mat']);
    ft_applymeasure(cfg)
end


basedir = '/group/northoff/share/hcp-meg/source';
rests = {'rest-3','rest-4','rest-5'};


for q = 1:3
    cd(fullfile(basedir,rests{q}))
    
    files = dir('*timeres.mat');
    names = extractfield(files,'name');
    files(contains(names,'grpdata')) = [];
            %specs = parload(files(1).name,'timespecs');

    resdata = cell(length(files),1);
    parfor i = 1:length(files)
        specs = parload(files(i).name,'timespecs');

        dat = zeros(1,length(specs),length(cfg.measure),size(specs{1}.mixd,2));
        for ii = 1:length(specs)
           for iii = 1:length(cfg.measure)
              dat(1,ii,iii,:) = cfg.measure{iii}(specs{ii}); 
           end
        end
        resdata{i} = dat;
    end
    numtpoints = cellfun(@(d)size(d,4),resdata,'UniformOutput',true);
    maxtpoints = max(numtpoints);
    for i = 1:length(resdata)
        if size(resdata{i},4) < maxtpoints
            resdata{i}(:,:,:,(end+1):maxtpoints) = NaN;
        end
    end
    resdata = cat(1,resdata{:});
    outputs.data = resdata;
    outputs.sub = extractfield(files,'name');
    save(fullfile(basedir,[rests{q} '_irasameas_individalpha_timeres.mat']),'outputs','-v7.3')
end