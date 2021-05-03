% script for getting IRASA measures in time-resolved specs

basedir = '/group/northoff/share/hcp-meg/source';
rests = {'rest-3','rest-4','rest-5'};


for q = 1:3
    cd(fullfile(basedir,rests{q}))
    
    files = dir('*timeres.mat');
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
    save(fullfile(basedir,[rests{q} '_camcan_IRASAmeas_abs_timeres.mat']),'outputs','-v7.3')
end