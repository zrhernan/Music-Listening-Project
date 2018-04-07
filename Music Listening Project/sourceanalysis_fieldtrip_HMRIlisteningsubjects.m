
%% Preprocessing of the anatomical MRI:
% read in MRI data
cd('Z:\Research_Projects\MUSIC_Karmonik\MUSIC_subject_01_20140410\FreeSurfer\mri')
mri = ft_read_mri('T1.nii');

% impose coordinate system according to MNI convention
cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'spm';
mri_spm    = ft_volumerealign(cfg, mri);

% reslicing
cfg            = [];
cfg.resolution = 1;
cfg.dim        = [256 256 256];
mri_spm_rs     = ft_volumereslice(cfg, mri_spm);
transform_vox2spm = mri_spm_rs.transform;
save('S9007_transform_vox2spm', 'transform_vox2spm');

% save the resliced anatomy in a FreeSurfer compatible format
cfg             = [];
cfg.coordsys    = 'spm';
cfg.filename    = 'S1_FTpreproc';
cfg.filetype    = 'nifti';
cfg.parameter   = 'anatomy';
ft_volumewrite(cfg, mri_spm_rs);

% impose coordinate system according to EEG convention
cfg          = [];
cfg.method   = 'interactive';
cfg.coordsys = 'ctf';
mri_ctf_rs   = ft_volumerealign(cfg, mri_spm_rs);
transform_vox2ctf = mri_ctf_rs.transform;
save('S9007_transform_vox2ctf', 'transform_vox2ctf');

save('S9007_transform_vox2ctf', 'transform_vox2ctf');
%% Source model:
% Volumetric processing in FreeSurfer
cd('Z:\Research_Projects\MUSIC_Karmonik\MUSIC_subject_01_20140410\FreeSurfer\mri')
mri = ft_read_mri('filled.nii');
 
cfg = [];
cfg.funparameter   = 'anatomy';
cfg.maskparameter  = cfg.funparameter;
cfg.interactive = 'yes';
ft_sourceplot(cfg, mri);

% Creation of the mesh using MNE Suite
cd(['Z:\Research_Projects\MUSIC_Karmonik\MUSIC_subject_01_20140410\FreeSurfer\bem'])

sourcespace = ft_read_headshape('S1-oct-6-src.fif', 'format', 'mne_source');

figure
ft_plot_mesh(sourcespace);

% Co-registration of the source space to the sensor-based head coordinate system
T = transform_vox2ctf/transform_vox2spm;

% go to the $SUBJECT/bem directory
sourcespace = ft_read_headshape('S1-oct-6-src.fif', 'format', 'mne_source');
sourcespace = ft_convert_units(sourcespace, 'mm');
sourcespace = ft_transform_geometry(T, sourcespace);

save sourcespace sourcespace

%% Generate volume conduction model
%{
% read in MRI data
cd('\\bmi-nas-01\Contreras-UH\NRI_Project_Nikunj\Experiment_Data\Clinical_study_Data\Subject_S9007\S9007_MRI\mri')
mri = ft_read_mri('S9007_FTpreproc.nii','fileformat','nifti');
load S9007_transform_vox2spm
mri.transform = transform_vox2spm;
mri.coordsys = 'spm';
% tissue segmentation
cfg           = [];
cfg.output    = {'brain','skull','scalp'};
segmentedmri  = ft_volumesegment(cfg, mri);

% create mesh
cfg = [];
cfg.tissue = {'brain','skull','scalp'};
cfg.method = 'projectmesh';
cfg.numvertices = [3000 2000 1000];
bnd = ft_prepare_mesh(cfg,segmentedmri);
%}
% load Freesurfer mesh files
cd('Z:\Research_Projects\MUSIC_Karmonik\MUSIC_subject_01_20140410\FreeSurfer\bem')
brain = ft_read_headshape('S1-inner_skull-5120.surf','fileformat','freesurfer_surf');
bnd(1).pnt = brain.pos;
bnd(1).tri = brain.tri;
bnd(1).unit = brain.unit;
skull = ft_read_headshape('S1-outer_skull-5120.surf','fileformat','freesurfer_surf');
bnd(2).pnt = skull.pos;
bnd(2).tri = skull.tri;
bnd(2).unit = skull.unit;
scalp = ft_read_headshape('S1-outer_skin-5120.surf','fileformat','freesurfer_surf');
bnd(3).pnt = scalp.pos;
bnd(3).tri = scalp.tri;
bnd(3).unit = scalp.unit;

%% create head model
cfg        = [];
cfg.method = 'bemcp'; % You can also specify 'openmeeg', 'bemcp', or another method (dipoli doesn't work).
vol        = ft_prepare_headmodel(cfg, bnd);
for i = 1:length(vol.bnd)
    vol.bnd(i).pnt = vol.bnd(i).pos;
end
save vol vol

%% Averaging and noise-covariance estimation and electrode position initialization
sesslist = {  '1',  '2',    '3',  '4',         '5',  '6',  '7',  '8',  '9',  '10', '11',  '12', '13', '14'};
subjlist = {'S9007','S9009','S9010','S9011','S9012','S9014'}; subj = 1;
cd('\\bmi-nas-01\Contreras-UH\NRI_Project_Nikunj\Experiment_Data\Clinical_study_Data\')
for sess = 1%:length(sesslist)
    %
    filename = ['Subject_',subjlist{subj},'\',subjlist{subj},'_Session',sesslist{sess},...
            '\',subjlist{subj},'_ses',sesslist{sess},'_preprocEEG.mat'];
    load(filename)
    cfg = [];
    cfg.covariance = 'yes';
    cfg.covariancewindow = [-2.5 -1.5]; % it will calculate the covariance matrix 
                                        % on the timepoints that are  
                                        % before the zero-time point in the trials
    datatlck = ft_timelockanalysis(cfg, datafinal);
    save tlck datatlck
    %}
    %
    filename1 = ['Subject_',subjlist{subj},'\',subjlist{subj},'_Session',sesslist{sess},...
            '\',subjlist{subj},'_session',sesslist{sess},'_electrodefile.elp'];
    elec = importfasttrak(filename1);
    emg_chans = {'TP7','TP8','TP9','TP10','FT7','FT8','FT9','FT10'};
    unused_chans = {'PO9','PO10'};
    elec.label = elec.label(~ismember((elec.label),[emg_chans,unused_chans]));
    elec.chanpos = elec.chanpos(~ismember((elec.label),[emg_chans,unused_chans]),:);
    elec.elecpos = elec.elecpos(~ismember((elec.label),[emg_chans,unused_chans]),:);
    save elec elec
    %
end

%% Forward Solution
load tlck;
load sourcespace;
load vol;
load elec;

%%
cfg = [];
cfg.elec = elec;                             % sensor positions
% cfg.grid.pos = sourcespace.pos;              % source points
% cfg.grid.inside = 1:size(sourcespace.pos,1); % all source points are inside of the brain
cfg.headmodel = vol;                               % volume conduction model
leadfield = ft_prepare_leadfield(cfg);

save leadfield leadfield;

