%% load all subject data
clearvars,clc
load 'ASRdenoisedEEG.mat'
freqband = [1 4; 4 7; 7 12; 12 15; 15 24; 24 30; 30 50];
Conds = {'selfselected','bach','gagaku','clicking','cronkite','chaplin'};
%% Note: need to work on spoken language next (06/15/2017)
timeavg_CondGroup1 = cell(11,1);     timeavg_CondGroup2 = cell(11,1);
subjcnt = 0;
for subj = [1:4,6:12]%[2:4,6:8,11]%
    subjcnt = subjcnt + 1;
    condspersubj = cell(6,1);
    noavggrp1flag = 0;      noavggrp2flag = 0;
    for cond = 1:6
        if isempty(preprocdata_allsubjANDcond{subj,cond}), continue, end
        data_preproc = preprocdata_allsubjANDcond{subj,cond};
        cfg                     = [];
        cfg.continuous          = 'yes';
        cfg.bpfilter            = 'yes';
        cfg.bpfilttype          = 'fir';
        cfg.bpfreq              = [4 7];
        cfg.hilbert             = 'abs';
        data_AM                 = ft_preprocessing(cfg,data_preproc);

        
    %{
    %% Visually inspect data
        %filter data first
        cfg                     = [];
        cfg.continuous          = 'yes';
        cfg.bpfilter            = 'yes';
        cfg.bpfilttype          = 'fir';
        cfg.bpfreq              = [4 7];
        cfg.dftfilter           = 'yes';
        % cfg.bpfiltdir           = 'twopass';
        % cfg.bpfiltord           = 4;
        % cfg.bpinstabilityfix    = 'reduce';
        cfg.hilbert             = 'abs';
        data_preproc_filt        = ft_preprocessing(cfg,data_preproc);
        cfg          = [];
        cfg.viewmode = 'vertical';
        ft_databrowser(cfg, data_preproc_filt);
    %}    
    %% Define event trials
    % focusing on 12 seconds before and after the listening onset (no onsets)
        r2lonset = (6000:12000:length(data_preproc.time{1}))'; %rest-to-listen onsets
        trl = [r2lonset - 500, r2lonset + 500, -500*ones(length(r2lonset),1)];
        trl(1:2:length(r2lonset),4) = ones(length(1:2:length(r2lonset)),1);
        trl(2:2:length(r2lonset),4) = 2*ones(length(2:2:length(r2lonset)),1);
        if trl(end,2) > length(data_preproc.time{1}), trl(end,2) = length(data_preproc.time{1}); end
        cfg           = [];
        cfg.trl       = trl;
%         datatrl       = ft_redefinetrial(cfg,data_preproc);
        datatrl       = ft_redefinetrial(cfg,data_AM);        
        condspersubj{cond} = datatrl;
    end
%% append each condition group per subject
    grp1_idcs = [1];%[6];
    grp2_idcs = [2 3];%[4];
    grp1_idcskept = grp1_idcs(~cellfun('isempty',{condspersubj{grp1_idcs}})); % music conditions
    grp2_idcskept = grp2_idcs(~cellfun('isempty',{condspersubj{grp2_idcs}})); % language conditions
    cfg          = [];
    numGrp1Condskept = length({condspersubj{grp1_idcskept}});
    numGrp2Condskept = length({condspersubj{grp2_idcskept}});
    if numGrp1Condskept > 1
        CondGroup1perSubj = ft_appenddata(cfg,condspersubj{grp1_idcskept});
    elseif numGrp1Condskept == 1 
        CondGroup1perSubj = condspersubj{grp1_idcskept};
    elseif numGrp1Condskept == 0
        noavggrp1flag = 1;
    end
    if numGrp2Condskept > 1
        CondGroup2perSubj = ft_appenddata(cfg,condspersubj{grp2_idcskept});
    elseif numGrp2Condskept == 1 
        CondGroup2perSubj = condspersubj{grp2_idcskept};
    elseif numGrp2Condskept == 0
        noavggrp2flag = 1;
    end
       
%% calculate average for each condition (listen block only)
    cfg                  = [];
    cfg.vartrllength     = 2;
    cfg.keeptrials       = 'yes';
    cfg.trials           = find(datatrl.trialinfo==2);
    if noavggrp1flag ~= 1
        timeavg_CondGroup1{subjcnt} = ft_timelockanalysis(cfg,CondGroup1perSubj);
    end
    if noavggrp2flag ~= 1
        timeavg_CondGroup2{subjcnt} = ft_timelockanalysis(cfg,CondGroup2perSubj);
    end
end

%% calculate grand average across subjects for each condition group   
cfg                 = [];
cfg.channel         = 'all';
cfg.latency         = 'all';
cfg.parameter       = 'avg';
cfg.keepindividual  = 'yes';
GACondGroup1        = ft_timelockgrandaverage(cfg,timeavg_CondGroup1{~cellfun('isempty',{timeavg_CondGroup1{:}})}); 
GACondGroup2        = ft_timelockgrandaverage(cfg,timeavg_CondGroup2{~cellfun('isempty',{timeavg_CondGroup2{:}})});
% "{:}" means to use data from all elements of the variable

%% visualize
%{
close all
% GACondGroup2_subjrmvd = GACondGroup2;
% GACondGroup2_subjrmvd.individual = GACondGroup2.individual([2:7 10],:,:);

%filter data first
cfg                     = [];
cfg.continuous          = 'yes';
% cfg.lpfilter            = 'yes';
% cfg.lpfilttype          = 'fir';
% cfg.lpfreq              = 11;
cfg.bpfilter            = 'yes';
cfg.bpfilttype          = 'fir';
cfg.bpfreq              = [4 7];
cfg.hilbert             = 'abs';
% GACondFamiliar_musiconly_filt        = ft_preprocessing(cfg,GACondGroup1);
% GACondUnfamiliar_musiconly_filt     = ft_preprocessing(cfg,GACondGroup2);
GACondFamiliar_filt        = ft_preprocessing(cfg,GACondGroup1);
GACondUnfamiliar_filt     = ft_preprocessing(cfg,GACondGroup2);
%}
%% Load LAY file of EEG channel 3D locations
load easycap16HMRImusicproj_layout
cfg              = [];
% cfg.xlim         = [0 1];
cfg.showlabels   = 'yes';
cfg.layout    	 = lay;
myfig=figure;
% ft_multiplotER(cfg,GACondFamiliar_filt, GACondUnfamiliar_filt);
ft_multiplotER(cfg,GACondGroup1, GACondGroup2);
legend({'Familiar Conditions (Self-Selected Only)','Unfamiliar Conditions (Bach and Gagaku Only)'})

%% Save the files!
cd('C:\Users\TMHZXH5\Documents\ERP Analysis Figures')
% filename = 'InstantaneousAmplitudePerChannelFrequency_FamiliarvsUnfamiliar_lowgammaband.tif';
% filename = 'MonteCarloSigTest_InstAmpPerChanFreq_FamiliarSelfSelectvsUnfamiliarBachGagaku_thetaband.tif';
filename = 'DifferenceInstAmpThetaBand_FamiliarSelfSelectvsUnfamiliarBachGagaku.tif';
set(gcf,'paperpositionmode','auto')
print(myfig,'-dpng','-r300',filename)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STATISTICAL TESTING: Time Only %%%%%%%
cfg = [];
numsubjs = intersect(numel(timeavg_CondGroup1),numel(timeavg_CondGroup1));
% numsubjs = intersect(size(timelock1.trial,1),size(timelock2.trial,1));
% cfg.design = [ 1*ones(1,numsubjs) 2*ones(1,numsubjs) ];
cfg.design(1,:) = [ 1*ones(1,numsubjs) 2*ones(1,numsubjs) ];
cfg.design(2,:) = [ 1:numsubjs 1:numsubjs ];
cfg.parameter = 'avg';
cfg.avgovertime = 'yes';
cfg.uvar = 2;
cfg.ivar = 1;
cfg.latency = [0 1];
cfg.method = 'montecarlo';
% cfg.statistic = 'indepsamplesT';
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';
% cfg.clusteralpha = 0.05;
% cfg.alpha = 0.1;
cfg.numrandomization = 2000;
cfg.neighbours = []; % only cluster over time, not over channels
% stat = ft_timelockstatistics(cfg, timelock1, timelock2);
TA_CG1_fullcells = {timeavg_CondGroup1{~cellfun('isempty',{timeavg_CondGroup1{:}})}};
TA_CG2_fullcells = {timeavg_CondGroup2{~cellfun('isempty',{timeavg_CondGroup2{:}})}};
stat = ft_timelockstatistics(cfg, TA_CG1_fullcells{:}, TA_CG2_fullcells{:});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STATISTICAL TESTING: Channels and Time %%%%%%
%% Get channel neighbors
cfg_neighb          = [];
cfg_neighb.method   = 'triangulation'; 
cfg_neighb.elec     = timeavg_CondGroup1{1,1}.elec;
cfg_neighb.feedback = 'yes';
neighbours          = ft_prepare_neighbours(cfg_neighb, timeavg_CondGroup1{1,1});
numneighbours        = cellfun(@length,{neighbours.neighblabel});
%% Monte Carlo simulation
cfg = [];
numsubjs = intersect(numel(timeavg_CondGroup1),numel(timeavg_CondGroup1));
cfg.design(1,:) = [ 1*ones(1,numsubjs) 2*ones(1,numsubjs) ];
cfg.design(2,:) = [ 1:numsubjs 1:numsubjs ];
cfg.parameter = 'avg';
% cfg.avgovertime = 'yes';
cfg.ivar = 1;
cfg.uvar = 2;
cfg.latency = [0.25 0.55];
cfg.ivar = 1;                     % number or list with indices indicating the independent variable(s)
cfg.method = 'montecarlo';        % Monte Carlo Method to calculate the significance probability
cfg.statistic = 'depsamplesT';  % dependent samples T-statistic at sample level
cfg.correctm = 'cluster';
% cfg.clusteralpha = 0.05;          %  sample-specific test statistic alpha level
% cfg.clusterstatistic = 'maxsum';  % test statistic evaluated via permutation distribution
cfg.minnbchan = 2;                % minimum number of neighboring channels required to be included in clustering algorithm
cfg.neighbours = neighbours;      % channel neighbour information
% cfg.tail = 0;                     % -1, 1 or 0 (default = 0); one-sided or two-sided test
% cfg.clustertail = 0;
cfg.correcttail = 'prob';
cfg.alpha = 0.05;                % alpha level of the permutation test
cfg.numrandomization = 2000;      % number of draws from the permutation distribution
TA_CG1_fullcells = {timeavg_CondGroup1{~cellfun('isempty',{timeavg_CondGroup1{:}})}};
TA_CG2_fullcells = {timeavg_CondGroup2{~cellfun('isempty',{timeavg_CondGroup2{:}})}};
stat = ft_timelockstatistics(cfg, TA_CG1_fullcells{:}, TA_CG2_fullcells{:});

% calculate the difference between conditions
% timelock1 = GACondFamiliar_filt;   timelock2 = GACondUnfamiliar_filt;
% timelock1.trial = timelock1.trial{1:11};
% timelock2.trial = timelock2.trial([2:7 10],:,:);
cfg = [];
cfg.operation = 'x1-x2';
% cfg.parameter = 'trial';
cfg.parameter = 'individual';
% difference = ft_math(cfg, timelock1, timelock2);
difference = ft_math(cfg, GACondGroup1, GACondGroup2);
%{
figure;
load easycap16HMRImusicproj_layout
cfg = [];
% cfg.xlim         = [0 1];
cfg.showlabels   = 'yes';
cfg.layout    	 = lay;
ft_multiplotER(cfg, difference);
legend({'difference'})
%}

% Plot t-values, condition differences, cluster probabilities, and significance binary values across time 
figure;
Tdiffintrsect = ismember(difference.time,stat.time);
xmin = min(difference.time(Tdiffintrsect));
xmax = max(difference.time(Tdiffintrsect));
subplot(4,1,1); plot(stat.time, stat.stat); ylabel('t-value'); xlim([xmin xmax])
subplot(4,1,2); plot(stat.time, squeeze(mean(difference.individual(:,:,Tdiffintrsect),1))); ylabel('Familiar-Unfamiliar (uV)'); xlim([xmin xmax])
subplot(4,1,3); semilogy(stat.time, stat.prob); ylabel('cluster probabilities'); axis([xmin xmax 0.01 1])
subplot(4,1,4); plot(stat.time, stat.mask,'.'); ylabel('binary significance value'); axis([xmin xmax -0.1 1.1])
%% Get info on probabilities
stat.label(sum(stat.mask,2) > 1)

minprob = min(stat.prob,[],2);    bb=[];
for ch=1:15, bb(ch,:)=ismember(stat.prob(ch,:),minprob(ch));end
bb = logical(bb);    statdetails = cell(16,4);
statdetails{1,1} = 'Label';     statdetails{1,2} = 'Time Start';
statdetails{1,3} = 'Time End';  statdetails{1,4} = 'Prob';
for ch=1:15
    statdetails{ch+1,1} = stat.label(ch);
    statdetails{ch+1,2} = min(stat.time(bb(ch,:)));
    statdetails{ch+1,3} = max(stat.time(bb(ch,:)));
    statdetails{ch+1,4} = min(stat.prob(ch,:));
end
%% plot difference w/regions of statisical significance
mask = zeros(length(stat.label),length(Tdiffintrsect)); mcntr=1;
for i = 1:length(Tdiffintrsect),
    if Tdiffintrsect(i)
mask(:,i) = stat.mask(:,mcntr); mcntr=mcntr+1;  end
end
difference.mask = mask;
myfig = figure;
cfg               = [];
cfg.showlabels    = 'yes';
cfg.xlim          = [0 1];
cfg.layout    	  = lay;
cfg.maskparameter = 'mask';
% cfg.maskstyle     = 'thickness';
cfg.renderer      = 'opengl';
ft_multiplotER(cfg, difference);
legend({'difference'})

%% per channel too
close all
for chn = 1:length(difference.label)
    close all

    cfg = [];
    cfg.maskparameter = 'mask';
    cfg.channel = difference.label(chn);
    cfg.showlabels   = 'yes';
    cfg.layout    	 = lay;
    cfg.xlim         = [0 1];
    ft_singleplotER(cfg, difference);
    set(gcf,'position',[200 300 416 634])
    
    % Save the files!
    cd(['C:\Users\TMHZXH5\Documents\ERP Analysis Figures\',...
        'InstantaneousAmplitudeThetaBand_FamiliarSelfSelectMinus',...
        'UnfamiliarBachDifference Plots by Channel_CORRECTEDclusterstats'])
    myfig = figure(1);
%     filename = ['DifferenceInstAmpThetaBand_FamiliarSelfSelectvsUnfamiliarBachGagaku_Chan',cfg.channel{:},'.tif'];
%     filename = ['DifferenceInstAmpThetaBand_FamiliarSelfSelectvsUnfamiliarBach_Chan',cfg.channel{:},'.tif'];
    filename = ['Chan',cfg.channel{:},'_correctedclustertesting.tif'];
%     filename = ['DifferenceInstAmpThetaBand_FamiliarvsUnfamiliar_Chan',cfg.channel{:},'.tif'];
    set(gcf,'paperpositionmode','auto')
    print(myfig,'-dpng','-r300',filename)
end

