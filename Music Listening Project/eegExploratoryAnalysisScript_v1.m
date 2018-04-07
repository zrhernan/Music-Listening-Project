% EEG Music Data Processing Script
clc, clear, close all
cd('Z:\Research_Projects\MUSIC_Karmonik\')
subjList = {'01_20140410'; '02_20140416'; '03_20140429'; '04_20140528';
            '05_20140529'; '06_20140610'; '07_20140620'; '08_20140710';
            '09_20140908'; '10_20140909';'11_20140922';'12_20141017'};
subjEEGfileprefixlist = {'MusicFavorite','ep_20140416_MUSICALL',...
    'ep_20140429_MUSICAL','ep_20140528_MUSICAL','ep_20140529_MUSICAL',...
    'ep_20140610_MUSICAL','ep_20140620_MUSICAL','ep_20140710_MUSICAL',...
    'ep_20140908_MUSICAL','ep_20140909_MUSICAL','ep_20140922_MUSICAL',...
    'ep_20141017_MUSICAL'};
SelfSelectNumList   = {'03','04','03','03','03','03','04','04','03','03','03','03'}; %04'};
BachNumList         = {'04','06','04','0405','04','04','05','06','04','04','04','0508'};
% subj = 11;
datafft = cell(12,1);
freqband = [1 4; 4 7; 7 12; 12 15; 15 24; 24 30; 30 50];
avg_pow = zeros(12,size(freqband,1));

%% Read in EEG Data
cfg                     = [];
cfg.dataset             = ['MUSIC_subject_',subjList{subj},'\EEG\',...
                            subjEEGfileprefixlist{subj},...
                            'Edit Markers_rem_sync_SSs_R-edf_run',...
                            SelfSelectNumList{subj},'.edf'];
cfg.continuous          = 'yes';
% cfg.detrend             = 'yes';
% cfg.bpfilter            = 'yes';
% cfg.bpfilttype          = 'but';
% cfg.bpfreq              = [0.1 100];
% cfg.bpfiltdir           = 'twopass';
% cfg.bpfiltord           = 4;
% cfg.bpinstabilityfix    = 'reduce';

data           = ft_preprocessing(cfg);

%% Re-reference to O1 channel and Remove ECG channel
newdata  = zeros(size(data.trial{1}));
newlabel = cell(length(data.label),1);
for cp = 1:length(data.label)
    newdata(cp,:) = sum(data.trial{1}(1:cp,:),1);
    newlabel{cp}  = regexprep(data.label{cp},'\w{2,3}_','');
end

data.trial{1}   = newdata;
data.label      = newlabel;
data.hdr.label  = newlabel;

% Remove ECG channel
cfg         = [];
cfg.channel = data.label(1:15);
data_noecg  = ft_selectdata(cfg,data);
%{
%% Define event trials
% music was played 24 seconds on and off
cfg           = [];
cfg.length    = 23.5;
datatrl       = ft_redefinetrial(cfg,data_noecg);
%{
%% De-noise data using ICA to remove artifactual ICs
[dComp, dataICclean] = uh_ica(datatrl,'reject');
%}
%% Use ICA in different way
% convert to EEGLAB structure
EEG=fieldtrip2eeglab(data_noecg,0);

% add EEG channel locations
load('BrainVision_1020_16ChanLocs.mat')
EEG.chanlocs = chanLocs;

% run ICA from EEGLAB function
EEGica = pop_runica(EEG,'icatype','runica');
EEGica.chanlocs = chanLocs;

% select components to remove
EEG_ICclean = pop_selectcomps(EEGica,(1:15));

% display just one component
% pop_prop(EEG,0,10)

%return 
dataICclean = eeglab2fieldtrip(EEG_ICclean,'componentanalysis','none');
%
%{
%% De-noise data by removal of EOG/EMG artifacts using FORCe method
[dataDenoise] = uh_force_denoisebyTrial(datatrl);
%
%% De-noise data w/Fieldtrip's zscore detection method
% Channel jumps
Data_NoJumps  = uh_removejumps(datatrl);

% EMG activity
% Data_NoEMG    = uh_removeemg(Data_NoJumps);
    
% EOG activity
Data_NoEOG  = uh_removeeog(Data_NoJumps);
%}

%% Visually inspect data
cfg          = [];
cfg.viewmode = 'vertical';
ft_databrowser(cfg, datatrl);
%}

%{
%% Perform Time-Frequency Analysis
cfg           = [];
cfg.output    = 'pow';
cfg.pad       = 'nextpow2';
cfg.channel   = 'all';
cfg.method    = 'wavelet';
cfg.width     = 14;
cfg.foi       = (1:40);
% cfg.t_ftimwin = 5./cfg.foi;
% cfg.tapsoffrq = 0.4 * cfg.foi;
cfg.toi       = data_noecg.time{1};
TFRwave       = ft_freqanalysis(cfg, data_noecg);
%}
%% Perform Fourier Spectral Analysis
cfg                 = [];
cfg.method          = 'mtmfft';
cfg.output          = 'pow';
cfg.pad             = 'nextpow2';
cfg.foi             = (1:0.1:100);
cfg.taper           = 'hanning';
datafft{subj}       = ft_freqanalysis(cfg, data_noecg);

for fq = 1:size(freqband,1)
    avg_pow(subj,fq) = bandpower(mean(datafft{subj}.powspctrm),datafft{subj}.freq, freqband(fq,:), 'psd');
end

%% Load LAY file of EEG channel 3D locations
load easycap16HMRImusicproj_layout

%% plot results
cfgp            = [];
cfgp.xlim       = [data_fft.freq(1) data_fft.freq(end)];
cfgp.ylim       = [fix(min(data_fft.powspctrm(:))) 1];%max(data_fft.powspctrm(:))];
cfgp.layout     = lay;
cfgp.showlabels = 'yes';
cfgp.interactive= 'yes';
ft_multiplotER(cfgp, data_fft);
%% Display Results
cfg = [];
cfg.baseline = [0 23.5];
cfg.baselinetype = 'absolute';
cfg.showlabels = 'yes';
cfg.layout = lay;
% cfg.zlim = 'minzero';
figure
ft_multiplotTFR(cfg, TFRwave)

%% Display on topoplot by 24 second epochs
figure
for i =2:13;
    cfg = [];
    cfg.baseline = [0.004 23.5];
    cfg.baselinetype = 'absolute';
    cfg.xlim = [datatrl.sampleinfo(i,1) datatrl.sampleinfo(i,2)]./datatrl.fsample;
    cfg.ylim = [0.1 3];
    cfg.zlim = [-200 200];
    cfg.marker = 'labels';
    cfg.layout = lay;
    cfg.comment = 'xlim';
    cfg.commentpos = 'middlebottom';
    subplot(6,2,i-1)
    ft_topoplotTFR(cfg, TFRwave)
end