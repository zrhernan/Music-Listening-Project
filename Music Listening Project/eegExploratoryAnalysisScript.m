% EEG Music Data Processing Script
clc, clear, close all
cd('Z:\Research_Projects\MUSIC_Karmonik\')
subjfile = 'MUSIC_subject_01_20140410\EEG\';
%% Read in EEG Data
cfg            = [];
cfg.dataset    = [subjfile,'MusicFavoriteEdit Markers_rem_sync_SSs_R-edf_run03.edf'];
cfg.continuous = 'yes';
cfg.channel    = 'all';
data           = ft_preprocessing(cfg);

%% Take out O2-ECG channel and replace with O2-O1 channel
o1_ecg = sum(data.trial{1}(2:16,:),1);
o2_ecg = data.trial{1}(16,:);
o2_o1 = o2_ecg - o1_ecg;
data.trial{1}(16,:) = o2_o1;
data.label{16} = 'O2_O1';
data.hdr.label{16} = 'O2_O1';

%% Check if EEG labels have underscore and correct if needed
if ~isempty(cell2mat(cellfun(@strfind,data.label,repmat({'_'},16,1),'uniformoutput',0)))
    data.label = cellfun(@strrep,data.label,repmat({'_'},16,1),repmat({'-'},16,1),'uniformoutput',0);
end

if ~isempty(cell2mat(cellfun(@strfind,data.hdr.label,repmat({'_'},16,1),'uniformoutput',0)))
    data.hdr.label = cellfun(@strrep,data.hdr.label,repmat({'_'},16,1),repmat({'-'},16,1),'uniformoutput',0);
end

%% Visually inspect data
cfg          = [];
cfg.viewmode = 'vertical';
ft_databrowser(cfg, datafinal);

%% Define event trials
% music was played 30 seconds on and off

%% Read in Channel Location Data
cfg          = [];
cfg.elecfile = 'MusicFavorite_Edit Markers_rem_sync_SSs_R_locs.sfp';
lay          = ft_prepare_layout(cfg);

%% Perform Time-Frequency Analysis
cfg           = [];
cfg.output    = 'pow';
cfg.channel   = 'all';
cfg.method    = 'wavelet';
cfg.width     = 7;
cfg.foi       = (1:2:50);
% cfg.t_ftimwin = 5./cfg.foi;
% cfg.tapsoffrq = 0.4 * cfg.foi;
cfg.toi       = data.time{1};
TFRwave       = ft_freqanalysis(cfg, data);

%% Display Results
cfg = [];
% cfg.baseline = [0 4];
% cfg.baselinetype = 'absolute';
% cfg.zlim = [-3e25 3e25];
cfg.showlabels = 'yes';
cfg.layout = lay;
figure
ft_multiplotTFR(cfg, TFRwave)

