%% Load pre-processed EEG data
clearvars
load('ASRdenoisedEEG.mat')

%%
diss_rest_listen = cell(12,6);
for subj = 2:12
    condcnt=0;
    data_preprocall = cell(1,length([preprocdata_allsubjANDcond{subj,:}]));
    for cond = 1:6
        if isempty(preprocdata_allsubjANDcond{subj,cond}), continue, end
        data_preproc = preprocdata_allsubjANDcond{subj,cond};
        condcnt = condcnt + 1;
        %% Add the reference O1 electrode back into set of electrodes
        data_preproc.label{13} = 'O1'; 
        data_preproc.label(14:16) = {'O2';'P4';'P8'};
        dattemp = data_preproc.trial{1}(13:15,:);
        data_preproc.trial{1}(13,:) = zeros(1,length(data_preproc.time{1}));
        data_preproc.trial{1}(14:16,:) = dattemp; clear dattemp
        load easycap16_elec
        data_preproc.elec = elec;
        data_preprocall{cond} = data_preproc;
    end
    %% 
    cfg = [];
    dataallcond = ft_appenddata(cfg,data_preprocall{:});
    dataallcond.trialinfo = [1;2;3;4;5;6];
        
        %% Compute Common Average Reference
        cfg             = [];
        cfg.reref       = 'yes';
        cfg.refchannel  = 'all';
        cfg.refmethod   = 'avg';       
        dataCAR         = ft_preprocessing(cfg,dataallcond);
        
        %% Compute Global Field Power (GFP) for each condition
        gfp = cell(1,length(dataallcond.trialinfo));
        for trl = dataallcond.trialinfo'
            data1trl    = ft_selectdata(dataCAR, 'rpt', trl);
            data1trl.avg = data1trl.trial;
            gfp{trl}    = ft_globalmeanfield([], data1trl);
        end
        %% Compute GFP-normalized waveforms per channel per condition
        % rest task
        for trl = dataallcond.trialinfo'
            
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