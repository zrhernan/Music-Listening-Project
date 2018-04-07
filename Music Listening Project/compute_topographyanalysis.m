%% Load pre-processed EEG data
clearvars
load('ASRdenoisedEEG.mat')

%%
diss_rest_listen = cell(12,6);
for subj = 1:12
    for cond = 1:6
        if isempty(preprocdata_allsubjANDcond{subj,cond}), continue, end
        data_preproc = preprocdata_allsubjANDcond{subj,cond};
        
        %% Add the reference O1 electrode back into set of electrodes
        data_preproc.label{13} = 'O1'; 
        data_preproc.label(14:16) = {'O2';'P4';'P8'};
        dattemp = data_preproc.trial{1}(13:15,:);
        data_preproc.trial{1}(13,:) = zeros(1,length(data_preproc.time{1}));
        data_preproc.trial{1}(14:16,:) = dattemp; clear dattemp
        load easycap16_elec
        data_preproc.elec = elec;
        
        %% Define event trials
        % focusing on 12 seconds before and after the listening onset
        r2lonset = (1:6000:length(data_preproc.time{1}))'; %rest-to-listen onsets
        trl = [r2lonset, r2lonset + 6000, 0*ones(length(r2lonset),1)];
        trl(1:2:length(r2lonset),4) = ones(length(1:2:length(r2lonset)),1);
        trl(2:2:length(r2lonset),4) = 2*ones(length(2:2:length(r2lonset)),1);
        if trl(end,2) > length(data_preproc.time{1}), trl(end,2) = length(data_preproc.time{1}); end
        cfg           = [];
        cfg.trl       = trl;
        datatrl       = ft_redefinetrial(cfg,data_preproc);
        
        %% Compute Common Average Reference
        cfg             = [];
        cfg.reref       = 'yes';
        cfg.refchannel  = 'all';
        cfg.refmethod   = 'avg';       
        dataCAR         = ft_preprocessing(cfg,datatrl);
        
        %% Take average across rest trials
        cfg              = [];
        cfg.vartrllength = 2;
        cfg.keeptrials   = 'yes';
        cfg.trials       = find(dataCAR.trialinfo==1);
        erp_rest         = ft_timelockanalysis(cfg,dataCAR);
        
        %% Take average across listen trials
        cfg              = [];
        cfg.vartrllength = 2;
        cfg.keeptrials   = 'yes';
        cfg.trials       = find(dataCAR.trialinfo==2);
        erp_listen       = ft_timelockanalysis(cfg,dataCAR);       
        
        
        %% Compute Global Field Power (GFP)
        cfg         = [];
        gfp_rest    = ft_globalmeanfield(cfg, erp_rest);
        gfp_listen  = ft_globalmeanfield(cfg, erp_listen);
        
        %% Compute GFP-normalized waveforms per channel
        % rest task
        gfpnorm_rest = erp_rest;
        gfprep_rest = repmat(gfp_rest.avg,length(erp_rest.label),1);
        gfpnorm_rest.avg = erp_rest.avg ./ gfprep_rest;
        
        % listening task
        gfpnorm_listen = erp_listen;
        gfprep_listen = repmat(gfp_listen.avg,length(erp_listen.label),1);
        gfpnorm_listen.avg = erp_listen.avg ./ gfprep_listen;
        
        %% Compute Global Dissimilarity (DISS)
        cfg                = [];
        cfg.operation      = '(x1 - x2)';
        cfg.parameter      = 'avg';
        gfpnormdiff        = ft_math(cfg, gfpnorm_rest, gfpnorm_listen);
        cfg                = [];
        diss_rest_listen{subj,cond} = ft_globalmeanfield(cfg, gfpnormdiff);
        
        %% Get channel neighbors
        cfg_neighb          = [];
        cfg_neighb.method   = 'triangulation'; 
        cfg_neighb.elec     = dataCAR.elec;
        cfg_neighb.feedback = 'yes';
        neighbours          = ft_prepare_neighbours(cfg_neighb, dataCAR);
        numneighbours        = cellfun(@length,{neighbours.neighblabel});
        %% Monte Carlo simulation       
        cfg = [];
        cfg.method = 'montecarlo';       % use the Monte Carlo Method to calculate the significance probability
        cfg.statistic = 'ft_statfun_indepsamplesT'; % use the independent samples T-statistic as a measure to 
                                         % evaluate the effect at the sample level
        cfg.correctm = 'cluster';
        cfg.clusteralpha = 0.05;         % alpha level of the sample-specific test statistic that 
                                         % will be used for thresholding
        cfg.clusterstatistic = 'maxsum'; % test statistic that will be evaluated under the 
                                         % permutation distribution. 
        cfg.minnbchan = 2;               % minimum number of neighborhood channels that is 
                                             % required for a selected sample to be included 
                                         % in the clustering algorithm (default=0).  or  ceil(0.4*numneighbours);%
        cfg.neighbours = neighbours;   % see below
        cfg.tail = 0;                    % -1, 1 or 0 (default = 0); one-sided or two-sided test
        cfg.clustertail = 0;
        cfg.alpha = 0.025;               % alpha level of the permutation test
        cfg.numrandomization = 100;      % number of draws from the permutation distribution

        design = zeros(1,size(gfpnorm_rest.trial,1) + size(gfpnorm_listen.trial,1));
        design(1,1:size(gfpnorm_rest.trial,1)) = 1;
        design(1,(size(gfpnorm_rest.trial,1)+1):(size(gfpnorm_rest.trial,1) + size(gfpnorm_listen.trial,1)))= 2;

        cfg.design = design;             % design matrix
        cfg.ivar  = 1;                   % number or list with indices indicating the independent variable(s)
        
        tanova = ft_timelockstatistics(cfg,gfpnorm_rest,gfpnorm_listen);
        
    end
end