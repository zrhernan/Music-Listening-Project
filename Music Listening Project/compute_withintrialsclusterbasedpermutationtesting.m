%% between-trials cluster-based permutation statistical test
subjList = {'01_20140410'; '02_20140416'; '03_20140429'; '04_20140528';
            '05_20140529'; '06_20140610'; '07_20140620'; '08_20140710';
            '09_20140908'; '10_20140909';'11_20140922';'12_20141017'};

cd(['Z:\Research_Projects\MUSIC_Karmonik\MUSIC_subject_',subjList{subj}])

%% Get channel neighbors
load easycap16_elec
elec.label          = elec.label(~ismember(elec.label,'O1'));
elec.chanpos        = elec.chanpos(~ismember(elec.label,'O1'),:);
cfg_neighb          = [];
cfg_neighb.method   = 'triangulation'; 
cfg_neighb.elec     = elec;
cfg_neighb.feedback = 'yes';
neighbours          = ft_prepare_neighbours(cfg_neighb);
% numneighbours        = cellfun(@length,{neighbours.neighblabel});

%% load data
load('timefreqresults_selfselected.mat')
TFRss = TFRwave;
load('timefreqresults_gagaku.mat')
TFRg = TFRwave;

%% Perform Monte Carlo Significance Testing
cfg = [];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_indepsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;
% prepare_neighbours determines what sensors may form clusters

elec.label          = elec.label(~ismember(elec.label,'O1'));
elec.chanpos        = elec.chanpos(~ismember(elec.label,'O1'),:);
cfg_neighb          = [];
cfg_neighb.method   = 'triangulation'; 
cfg_neighb.elec     = elec;
% cfg_neighb.feedback = 'yes';
cfg.neighbours      = ft_prepare_neighbours(cfg_neighb);


design = zeros(1,size(TFRss.powspctrm,1) + size(TFRg.powspctrm,1));
design(1,1:size(TFRss.powspctrm,1)) = 1;
design(1,(size(TFRss.powspctrm,1)+1):(size(TFRss.powspctrm,1)+...
  size(TFRg.powspctrm,1))) = 2;

cfg.design           = design;
cfg.ivar             = 1;

[stat] = ft_freqstatistics(cfg, TFRss, TFRg);


%% Display Results

% average over trials
cfg = []; 
TFRss = ft_freqdescriptives(cfg, TFRss);
TFRg  = ft_freqdescriptives(cfg, TFRg);

% add raw effect (selfselect - gagaku)
stat.raweffect = TFRss.powspctrm - TFRg.powspctrm;

% Load LAY file of EEG channel 3D locations
load easycap16HMRImusicproj_layout

cfg = [];
cfg.alpha  = 0.025;
cfg.parameter = 'raweffect';
% cfg.zlim   = [-1e-27 1e-27];
cfg.layout = lay;
ft_clusterplot(cfg, stat);
