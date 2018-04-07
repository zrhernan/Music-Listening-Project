clearvars,clc
load 'ASRdenoisedEEG.mat'
freqband = [1 4; 4 7; 7 12; 12 15; 15 24; 24 30; 30 50];
glob_rel_pow = zeros(12,size(freqband,1),6);
%
for subj = 1:12
    for cond = 1:6
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
        
%% Perform Fourier Spectral Analysis
        cfg                 = [];
        cfg.trials          = find(datatrl.trialinfo==2);
        cfg.method          = 'mtmfft';
        cfg.output          = 'pow';
        cfg.pad             = 'nextpow2';
        cfg.foi             = (0.1:0.1:100);
        cfg.taper           = 'hanning';
        cfg.tapsmofrq       = 2;
        datafft             = ft_freqanalysis(cfg, datatrl);
        
        for fbndoi = 1:size(freqband,1)
            glob_rel_pow(subj,fbndoi,cond) = bandpower(mean(datafft.powspctrm),...
                datafft.freq, freqband(fbndoi,:), 'psd') ./ ...
                bandpower(mean(datafft.powspctrm),datafft.freq, [1 50], 'psd');
            for chn = 1:length(datafft.label) % to save as extracted features
                relpowbands.delta(subj,cond,chn) = bandpower(datafft.powspctrm(chn,:),...
                datafft.freq, freqband(1,:), 'psd') ./ bandpower(datafft.powspctrm(chn,:),...
                datafft.freq, [1 50], 'psd');
                relpowbands.theta(subj,cond,chn) = bandpower(datafft.powspctrm(chn,:),...
                datafft.freq, freqband(2,:), 'psd') ./ bandpower(datafft.powspctrm(chn,:),...
                datafft.freq, [1 50], 'psd');
                relpowbands.loalpha(subj,cond,chn) = bandpower(datafft.powspctrm(chn,:),...
                datafft.freq, freqband(3,:), 'psd') ./ bandpower(datafft.powspctrm(chn,:),...
                datafft.freq, [1 50], 'psd');        
                relpowbands.hialpha(subj,cond,chn) = bandpower(datafft.powspctrm(chn,:),...
                datafft.freq, freqband(4,:), 'psd') ./ bandpower(datafft.powspctrm(chn,:),...
                datafft.freq, [1 50], 'psd');
                relpowbands.lobeta(subj,cond,chn) = bandpower(datafft.powspctrm(chn,:),...
                datafft.freq, freqband(5,:), 'psd') ./ bandpower(datafft.powspctrm(chn,:),...
                datafft.freq, [1 50], 'psd');        
                relpowbands.hibeta(subj,cond,chn) = bandpower(datafft.powspctrm(chn,:),...
                datafft.freq, freqband(6,:), 'psd') ./ bandpower(datafft.powspctrm(chn,:),...
                datafft.freq, [1 50], 'psd');
                relpowbands.logamma(subj,cond,chn) = bandpower(datafft.powspctrm(chn,:),...
                datafft.freq, freqband(7,:), 'psd') ./ bandpower(datafft.powspctrm(chn,:),...
                datafft.freq, [1 50], 'psd');
                relpowbands.broadband(subj,cond,chn) = bandpower(datafft.powspctrm(chn,:),...
                datafft.freq, [1 50], 'psd');            
            end
        end
    end
end

%% Save Extracted Features and Grand-Averaged Mean Power
cd('Z:\Research_Projects\MUSIC_Karmonik\')
save('global_rel_pow_ASRdenoisedEEG_onlylistenphase.mat','glob_rel_pow')
save('relpowbands_ASRdenoisedEEG_onlylistenphase.mat','relpowbands')

%% Plot Grand-Averaged Mean Power
plot_lines_grandavgpower_allsubj_allcond();

%% Plot Similarity Matrix
% plot_similaritymatrix_extractedfeats_allsubj_allmusicconds();
plot_similaritymatrix_extractedfeats_allsubj_allmusicconds_arrangedbysubj();