load 'ASRdenoisedEEG.mat'
subjcondlist = {};
psdgrandarray = [];
psdminusRESTgrandarray = [];
for subj = 1:12
    for cond = 1:6
        if isempty(preprocdata_allsubjANDcond{subj,cond}), continue, end
        data_preproc = preprocdata_allsubjANDcond{subj,cond};
        
%% Define event trials
        % focusing on 12 seconds before and after the listening onset (no onsets)
        r2lonset = (1:6000:length(data_preproc.time{1}))'; %rest-to-listen onsets
        trl = [r2lonset, r2lonset + 5999, 0*ones(length(r2lonset),1)];
        trl(1:2:length(r2lonset),4) = ones(length(1:2:length(r2lonset)),1);
        trl(2:2:length(r2lonset),4) = 2*ones(length(2:2:length(r2lonset)),1);
        if trl(end,2) > length(data_preproc.time{1}), trl(end,2) = length(data_preproc.time{1}); end
        cfg           = [];
        cfg.trl       = trl;
        datatrl       = ft_redefinetrial(cfg,data_preproc);

%% Perform Fourier Spectral Analysis
% for listening
        cfg                 = [];
        cfg.trials          = find(datatrl.trialinfo==2);
        cfg.method          = 'mtmfft';
        cfg.output          = 'pow';
        cfg.pad             = 'nextpow2';
        cfg.foi             = (0.1:0.1:50);
        cfg.taper           = 'hanning';
        cfg.tapsmofrq       = 2;
        datafftLISTEN       = ft_freqanalysis(cfg, datatrl);
        
        
% for resting
        cfg                 = [];
        cfg.trials          = find(datatrl.trialinfo==1);
        cfg.method          = 'mtmfft';
        cfg.output          = 'pow';
        cfg.pad             = 'nextpow2';
        cfg.foi             = (0.1:0.1:50);
        cfg.taper           = 'hanning';
        cfg.tapsmofrq       = 2;
        datafftREST         = ft_freqanalysis(cfg, datatrl);
      
% Then take the difference of each PSD using ft_math
        cfg                 = [];
        cfg.operation       = 'subtract';
        cfg.parameter       = 'powspctrm';
        datafftLISTENvsREST = ft_math(cfg,datafftLISTEN,datafftREST); 
        
% save for each subject and condition
        psdgrandarray          = [psdgrandarray; datafftLISTEN.powspctrm];
        psdminusRESTgrandarray = [psdminusRESTgrandarray; datafftLISTENvsREST.powspctrm];

        for chn = 1:length(datafftLISTEN.label) % to save as extracted features
            subjcondlist = [subjcondlist; ['s',num2str(subj),'c',num2str(cond),datafftLISTEN.label{chn}]];
 
        end
    end
end
%% Open large file of PSDs
load allPSDs_ASRdenoisedEEG.mat

%% extract only one channel
% chnidx = find(ismember(datafft.label,'C4'));
close all
PSDs = psdminusRESTgrandarray(:,121:149); %psdgrandarray;  %
for chn = 1:15
    chnidx = chn;
    chnidcs = (chnidx:15:size(PSDs,1));
    chnpsds = (PSDs(chnidcs,:));
    kclust = 2;
    [label, model, energy] = knKmeans(chnpsds',kclust);
    plot_cluster(kclust,label,chnpsds,freqbins(121:149))
    condWithLabel = plot_piepercluster(kclust,label,subjcondlist,chnidcs);
end
    %% Get classification rate to condition
    finalpred = zeros(6,2);
    for actC = 1:6
        prctpred = zeros(6,1);
        for predC = 1:kclust
            prctpred(predC) = sum(predconds(condnums==actC)==predC) ./ numel(predconds(condnums==actC));
        end
        [M, i] = max(prctpred);
        finalpred(actC,:) = [i, M];
    end

    numOfActConds(chn) = numel(unique(finalpred(:,1)))./kclust;
% end
%% Display results
figure(12)
for i = 1:length(chnidcs)
    scatter(freqbins,log(chnpsds(i,:)),[],predconds)
    hold on
end

%% 
figure(12)
for i = 1:length(freqbins)
    scatter((1:length(subjcondlist(chnidcs))),log(chnpsds(:,i)),[],predidcs)
    hold on
end
set(gca,'xtick',(1:length(subjcondlist(T8idcs))),'xticklabel',subjcondlist(T8idcs))
rotateXLabels(gca,45)