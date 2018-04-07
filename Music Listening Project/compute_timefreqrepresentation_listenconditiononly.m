clearvars,clc
load 'ASRdenoisedEEG.mat'
subjList = {'01_20140410'; '02_20140416'; '03_20140429'; '04_20140528';
            '05_20140529'; '06_20140610'; '07_20140620'; '08_20140710';
            '09_20140908'; '10_20140909';'11_20140922';'12_20141017'};
freqband = [1 4; 4 7; 7 12; 12 15; 15 24; 24 30; 30 50];
glob_rel_pow = zeros(12,size(freqband,1),6);
Conds = {'selfselected','bach','gagaku','clicking','cronkite','chaplin'};

%%
for subj = 1%:12
    for cond = [1 3]%:6
        if isempty(preprocdata_allsubjANDcond{subj,cond}), continue, end
        data_preproc = preprocdata_allsubjANDcond{subj,cond};
       
        %% Define event trials
        % focusing on 12 seconds before and after the listening onset (no onsets)
        r2lonset = (1:6000:length(data_preproc.time{1}))'; %rest-to-listen onsets
        trl = [r2lonset, r2lonset + 6000, 0*ones(length(r2lonset),1)];
        trl(1:2:length(r2lonset),4) = ones(length(1:2:length(r2lonset)),1);
        trl(2:2:length(r2lonset),4) = 2*ones(length(2:2:length(r2lonset)),1);
        if trl(end,2) > length(data_preproc.time{1}), trl(end,2) = length(data_preproc.time{1}); end
        cfg           = [];
        cfg.trl       = trl;
        datatrl       = ft_redefinetrial(cfg,data_preproc);
        
        %% Perform Time-Frequency Analysis (using multitapers)
        cfg = [];
        cfg.output     = 'pow';
        cfg.pad        = 'nextpow2';
        cfg.trials     = find(datatrl.trialinfo==2);
        cfg.method     = 'mtmconvol';%'wavelet';%
        cfg.foi        = 4;%(1:0.1:4);
        cfg.toi        = datatrl.time{1};
        cfg.keeptrials = 'yes';
        cfg.t_ftimwin  = 5./cfg.foi;
        cfg.tapsmofrq  = 0.4*cfg.foi;  
        TFRwave = ft_freqanalysis(cfg, datatrl);
%{        
        %% Display the results
        % Load LAY file of EEG channel 3D locations
        load easycap16HMRImusicproj_layout
        cfg = [];	
        cfg.showlabels   = 'yes';	        
        cfg.layout       = lay;
        figure
        ft_multiplotTFR(cfg, TFRwave)        
%}        
        %% Save Extracted Features and Grand-Averaged Mean Power
        cd(['Z:\Research_Projects\MUSIC_Karmonik\MUSIC_subject_',subjList{subj}])
        save(['timefreqresults_',Conds{cond},'.mat'],'TFRwave','-v7.3') 
    end
end
