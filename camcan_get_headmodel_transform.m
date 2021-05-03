% script for creating headmodels and getting sourcemodel transformation
% matrices

basedir = '/group/northoff/share/camcan/mindboggle/mindboggle/release001/freesurfer_subjects';
cd(basedir)

files = dir('*');
files = files(5:end);

fscommand = 'mri_vol2vol --mov T1.mgz --targ rawavg.mgz --regheader --o T1-in-rawavg.mgz --no-save-reg';

parfor i = 1:length(files)
    try
        if exist(fullfile(basedir,files(i).name,'sourcemodel'),'dir')
        cd(fullfile(basedir,files(i).name,'mri'))
        % convert freesurfer intensity-normalized image to native space 
        system(['export SUBJECTS_DIR=' basedir newline 'export SUBJECTNAME=' files(i).name newline fscommand]);
        
        % read mri, fiducials, and coregister to neuromag coordinates
        mri_native = ft_read_mri('T1-in-rawavg.mgz');
        cfg = []; cfg.method = 'fiducial';
        fid = parload(['/group/northoff/share/camcan/fiducials/fid-native-' files(i).name '.mat'],'fid');
        cfg.fiducial.nas = fid.native.vox.nas;
        cfg.fiducial.lpa = fid.native.vox.lpa;
        cfg.fiducial.rpa = fid.native.vox.rpa;
        cfg.coordsys = 'neuromag';
        mri_nmg = ft_volumerealign(cfg,mri_native);
        transform_vox2nmg = mri_nmg.transform;
        
        % segment mri and create headmodel
        cfg = []; cfg.output = 'brain';
        mri_seg_nmg = ft_volumesegment(cfg,mri_nmg);
        cfg = []; cfg.method = 'singleshell';
        vol = ft_prepare_headmodel(cfg,mri_seg_nmg);
        parsave(fullfile(basedir,files(i).name,'sourcemodel',[files(i).name '_headmodel.mat']),'vol',vol)
        
        % create transform matrix to neuromag
        transform_vox2camcan = mri_native.transform;
        transform_camcan2nmg = transform_vox2nmg/transform_vox2camcan;
        
        % transform sourcemodel to neuromag coordinates
        sourcemodel = parload(fullfile(basedir,files(i).name,'sourcemodel',[files(i).name '_sourcemodel_8k.mat']),'sourcemodel');
        sourcemodel = ft_transform_geometry(transform_camcan2nmg,sourcemodel);
        if isfield(sourcemodel,'atlasroi')
            sourcemodel.inside = sourcemodel.atlasroi>0;
            sourcemodel = rmfield(sourcemodel, 'atlasroi');
        end
        sourcemodel.coordsys = 'neuromag';
        parsave(fullfile(basedir,files(i).name,'sourcemodel',[files(i).name '_sourcemodel_8k_nmg.mat']),'sourcemodel',sourcemodel)
        end
        catch errormsg
        parsave(fullfile(basedir,[files(i).name '_error.mat']),'errormsg',errormsg)
    end
end