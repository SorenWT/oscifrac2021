% script to create sourcemodels for camcan subjects

addpath('/group/northoff/share/fieldtrip-master')
ft_defaults
addpath(genpath('~/Documents/MATLAB'))

cd /group/northoff/share/camcan/mindboggle/mindboggle/release001/freesurfer_subjects/

files = dir('CC*');
names = extractfield(files,'name');

taskfiles = dir('/home/soren/Documents/camcan/Preprocessed/Task/Epoched/*long_1Hz.mat');

tasknames = extractfield(taskfiles,'name');
tasknames = cellstr(extractBefore(tasknames,'_epoched'));

[m1,m2] = match_str(names,tasknames);

files = files(m1);
names = names(m1);

parpool(18)

freesurferdir = '/group/northoff/share/camcan/mindboggle/mindboggle/release001/freesurfer_subjects/';

badsubs = zeros(1,length(files));

parfor i = 1:length(files)
    try
    system(['/group/northoff/share/fieldtrip-master/bin/ft_postfreesurferscript.sh '...
        freesurferdir ' '...
        names{i} ' /home/soren/HCPpipelines/global/templates/standard_mesh_atlases']);
    system(['mkdir ' freesurferdir names{i} '/sourcemodel']);
    system(['cp ' freesurferdir names{i} '/workbench/*midthickness.4k* ' freesurferdir names{i} '/sourcemodel'])
    sourcefiles = dir(fullfile(freesurferdir,names{i},'sourcemodel','*gii'));
    sourcemodel = ft_read_headshape({fullfile(sourcefiles(1).folder,sourcefiles(1).name) fullfile(sourcefiles(2).folder,sourcefiles(2).name)});
    parsave(fullfile(freesurferdir,names{i},'sourcemodel',[names{i} '_sourcemodel_8k.mat']),'sourcemodel',sourcemodel)
    system(['rm -r ' freesurferdir names{i} '/workbench'])
catch
badsubs(i) = 1;
end
end




