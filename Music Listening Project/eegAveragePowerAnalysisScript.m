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
SelfSelectList      = { '03','04','03','03', '03','03','04','04', '03', '03','03', '03'}; %option '04' for subj 12
BachList            = { '04','06','04','05', '04','04','05','06', '04', '04','04', '05'}; %option '08 for subj 12
GagakuList          = { '07','08','08','07', '07','07','06','07', '07', '06','07', '09'};
ClickList           = { '09','09','09','08','NaN','08','07','08', '08', '07','08', '10'};
CronkiteList        = {'NaN','10','10','09','NaN','09','08','09','NaN','NaN','09','NaN'};
ChaplinList         = {'NaN','11','11','10','NaN','10','09','10','NaN','NaN','10','NaN'};

CondList = {SelfSelectList,BachList,GagakuList,ClickList,CronkiteList,ChaplinList};

datafft = cell(12,1);
freqband = [1 4; 4 7; 7 12; 12 15; 15 24; 24 30; 30 50];
avg_pow = zeros(12,size(freqband,1),length(CondList));
%%
for subj = 1:12
    for cond = 1:length(CondList)
        cond_idxlist = CondList{cond};
        if isnan(str2double(cond_idxlist{subj})), continue, end
        filename = ['MUSIC_subject_',subjList{subj},'\EEG\',subjEEGfileprefixlist{subj},...
                    'Edit Markers_rem_sync_SSs_R-edf_run',cond_idxlist{subj},'.edf'];
                
        %% Read in EEG Data
        cfg                     = [];
        cfg.dataset             = filename;
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
        %
        %% De-noise data using ICA to remove artifactual ICs
        [dComp, dataICclean] = uh_ica(data_noecg,'reject');
        %
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
        %
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
        %
        
        %% De-noise data using ASR
        eegplot_flag = 0;
        [EEG] = fieldtrip2eeglab(data_noecg,eegplot_flag);
        
        % run ASR
        ASRdata = uh_runASR_NF(EEG);
%         vis_artifacts(ASRdata,EEG)

        %return to fieldtrip
        data_asr = eeglab2fieldtrip(ASRdata,'preprocessing','none');
        
        if isfield(ASRdata.etc,'clean_channel_mask')
            chanlocs_labels = {EEG.chanlocs.labels};
            missingchans = chanlocs_labels(ASRdata.etc.clean_channel_mask == 0);
            % interpolate missing channels
            data_preproc = uh_removechannel(data_asr,'interpolate',missingchans',0);            
        end
        %{
        %% Visually inspect data
        cfg          = [];
        cfg.viewmode = 'vertical';
        ft_databrowser(cfg, data_interp);
        %}
        
        %% Calculate some Time-Series Statistical Features
        extractedfeats.max(subj,cond,:)  = max(data_preproc.trial{1},[],2);         % maximum values per channel
        extractedfeats.min(subj,cond,:)  = min(data_preproc.trial{1},[],2);         % minimum values per channel
        extractedfeats.var(subj,cond,:)  = var(data_preproc.trial{1},0,2);          % variance values per channel
        extractedfeats.skew(subj,cond,:) = skewness(data_preproc.trial{1},0,2);     % skewness values per channel
        extractedfeats.kurt(subj,cond,:) = kurt(data_preproc.trial{1}');            % kurtosis values per channel
        extractedfeats.iqr(subj,cond,:)  = iqr(data_preproc.trial{1},2);            % interquartile range per channel
        extractedfeats.mad(subj,cond,:)  = mad(data_preproc.trial{1},0,2);          % median absolute deviation values per channel

        %% Adding entropy measure as well
        
        [EL,BIN] = hist(data_preproc.trial{1}',100);
        pdf_allchns = EL./(numel(EL)*diff(BIN(1:2)));
        H = nan(length(data_preproc.label),1);
        for ch = 1:length(data_preproc.label)
            chnPr = pdf_allchns(:,ch);
            H(ch) = sum(-(chnPr(chnPr>0).*(log2(chnPr(chnPr>0)))));
        end
        extractedfeats.entpy(subj,cond,:)  = H;          % Shannon entropy values per channel
        %{
        %% and percentiles
        extractedfeats.prctile5(subj,cond,:) = prctile(zscore(data_preproc.trial{1}'),5); %5th percentile
        extractedfeats.prctile25(subj,cond,:) = prctile(zscore(data_preproc.trial{1}'),25); %25th percentile
        extractedfeats.prctile50(subj,cond,:) = prctile(zscore(data_preproc.trial{1}'),50); %50th percentile
        extractedfeats.prctile75(subj,cond,:) = prctile(zscore(data_preproc.trial{1}'),75); %75th percentile
        extractedfeats.prctile95(subj,cond,:) = prctile(zscore(data_preproc.trial{1}'),95); %95th percentile
        %}  
        % not yet, need to check Jesus' feature extraction code on descretization
        %% Perform Fourier Spectral Analysis
        cfg                 = [];
        cfg.method          = 'mtmfft';
        cfg.output          = 'pow';
        cfg.pad             = 'nextpow2';
        cfg.foi             = (1:0.1:100);
        cfg.taper           = 'hanning';
        cfg.tapsmofrq       = 0.4 * cfg.foi;
        datafft{subj}       = ft_freqanalysis(cfg, data_preproc);

        for fq = 1:size(freqband,1)
            avg_pow(subj,fq,cond) = bandpower(mean(datafft{subj}.powspctrm),...
                datafft{subj}.freq, freqband(fq,:), 'psd');
            for chn = 1:length(datafft{subj}.label) % to save as extracted features
                extractedfeats.delt_mupow(subj,cond,chn) = bandpower(datafft{subj}.powspctrm(chn,:),...
                datafft{subj}.freq, freqband(1,:), 'psd');
                extractedfeats.thet_mupow(subj,cond,chn) = bandpower(datafft{subj}.powspctrm(chn,:),...
                datafft{subj}.freq, freqband(2,:), 'psd');
                extractedfeats.loalph_mupow(subj,cond,chn) = bandpower(datafft{subj}.powspctrm(chn,:),...
                datafft{subj}.freq, freqband(3,:), 'psd');        
                extractedfeats.hialph_mupow(subj,cond,chn) = bandpower(datafft{subj}.powspctrm(chn,:),...
                datafft{subj}.freq, freqband(4,:), 'psd');
                extractedfeats.lobeta_mupow(subj,cond,chn) = bandpower(datafft{subj}.powspctrm(chn,:),...
                datafft{subj}.freq, freqband(5,:), 'psd');        
                extractedfeats.hibeta_mupow(subj,cond,chn) = bandpower(datafft{subj}.powspctrm(chn,:),...
                datafft{subj}.freq, freqband(6,:), 'psd');
                extractedfeats.logamm_mupow(subj,cond,chn) = bandpower(datafft{subj}.powspctrm(chn,:),...
                datafft{subj}.freq, freqband(7,:), 'psd');
                extractedfeats.broadband_mupow(subj,cond,chn) = bandpower(datafft{subj}.powspctrm(chn,:),...
                datafft{subj}.freq, [1 100], 'psd');            
            end
        end
         
    end
end
%{
%% Load LAY file of EEG channel 3D locations
load easycap16HMRImusicproj_layout

%% plot Fourier spectra as PSDs
cfgp            = [];
cfgp.xlim       = [datafft{subj}.freq(1) datafft{subj}.freq(end)];
cfgp.ylim       = [fix(min(datafft{subj}.powspctrm(:))) 1];%max(datafft{subj}.powspctrm(:))];
cfgp.layout     = lay;
cfgp.showlabels = 'yes';
cfgp.interactive= 'yes';
ft_multiplotER(cfgp, datafft{subj});
%}

%% Save Extracted Features and Grand-Averaged Mean Power
save('grand_avg_pow_ASRdenoisedEEG.mat',avg_pow)
save('extractfeats_ASRdenoisedEEG.mat',extractedfeats)

%% Plot Grand-Averaged Mean Power


%% Plot Similarity Matrix

