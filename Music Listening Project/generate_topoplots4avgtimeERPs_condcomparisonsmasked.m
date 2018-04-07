%% load all subject data
clearvars,clc
addpath('C:\Users\zrhernan\Documents\MATLAB\Listening Project Matlab Code\Imported Matlab Data')
load 'ASRdenoisedEEG.mat'
freqband = [1 4; 4 7; 7 12; 12 15; 15 24; 24 30; 30 50];
AudSampList = {'Self-Selected','Bach','Gagaku','Clicking','Cronkite','Chaplin'};
AudSampListofSubjs = cell(6,1);
timeavg_AudSampList = cell(6,1);

%%
tic;
for audsamp = 1:6
    subjsperaudsamp = cell(11,1);
    subjcnt = 0;
    for subj = [1:4,6:12]%[2:4,6:8,11]%
        subjcnt = subjcnt + 1;
    
        if isempty(preprocdata_allsubjANDcond{subj,audsamp}), continue, end
        data_preproc = preprocdata_allsubjANDcond{subj,audsamp};
        cfg                     = [];
        cfg.continuous          = 'yes';
        cfg.bpfilter            = 'yes';
        cfg.bpfilttype          = 'fir';
        cfg.bpfreq              = [4 7];
        cfg.hilbert             = 'abs';
        data_AM                 = ft_preprocessing(cfg,data_preproc);
    
    %% Define event trials
    % focusing on 2 seconds before and after the listening onset (no onsets)
        r2lonset = (6000:12000:length(data_preproc.time{1}))'; %rest-to-listen onsets
        trl = [r2lonset - 500, r2lonset + 500, -500*ones(length(r2lonset),1)];
        trl(1:2:length(r2lonset),4) = ones(length(1:2:length(r2lonset)),1);
        trl(2:2:length(r2lonset),4) = 2*ones(length(2:2:length(r2lonset)),1);
        if trl(end,2) > length(data_preproc.time{1}), trl(end,2) = length(data_preproc.time{1}); end
        cfg           = [];
        cfg.trl       = trl;
%         datatrl       = ft_redefinetrial(cfg,data_preproc);
        datatrl       = ft_redefinetrial(cfg,data_AM);        
        subjsperaudsamp{subjcnt} = datatrl;
    end
%% append all subjects per condition
    subjsperaudsamp_noempties = subjsperaudsamp(~cellfun(@isempty,{subjsperaudsamp{:}}));
    AudSampListofSubjs{audsamp} = ft_appenddata(cfg,subjsperaudsamp_noempties{:});
      
%% calculate average for each condition epoch
    cfg                  = [];
    cfg.vartrllength     = 2;
    cfg.keeptrials       = 'yes';
    cfg.trials           = find(datatrl.trialinfo==2);
    timeavg_AudSampList{audsamp} = ft_timelockanalysis(cfg,AudSampListofSubjs{audsamp});

end
timecompute(toc);

%% Display topoplots as subplots
close all
cd('C:\Users\zrhernan\Documents\MATLAB\Listening Project Matlab Code')
load('EEG Montage Templates\BrainVision_1020_64ChanLocs.mat')
load('Imported Matlab Data\ERPdata_GroupedbyAuditorySample.mat')
chanLocs = chanLocs((ismember({chanLocs.labels},strrep(data_AM.label,'FP','Fp'))));
statlist{1} = load('Imported Matlab Data\ERPdata_FamiliarSelfSelectedvsUnfamiliarBach_correctedcluststats_v2.mat', 'stat');
statlist{2} = load('Imported Matlab Data\ERPdata_FamiliarSelfSelectedvsUnfamiliarGagaku_correctedcluststats_v2.mat', 'stat');
statlist{1}.condlabel = {'Self-Selected','Bach'};
statlist{1}.condidx = [1 2];
statlist{2}.condlabel = {'Self-Selected','Gagaku'};
statlist{2}.condidx = [1 3];

% Extract indices to match channel order in 'chanLocs' variable
chnIDCS = zeros(length(chanLocs),1);
for chn = 1:length(chanLocs)
    labelsREVSD = strrep(data_AM.label,'FP','Fp');
    chnIDCS(chn) = find(ismember(labelsREVSD,{chanLocs(chn).labels}));
end

% Get global and local (to each auditory sample) max and min
localmax = nan(length(timeavg_AudSampList),1);      localmin = nan(size(localmax));
for audsamp = 1:length(timeavg_AudSampList)
    localmax(audsamp) = max(max(timeavg_AudSampList{audsamp}.avg(:)));
    localmin(audsamp) = min(min(timeavg_AudSampList{audsamp}.avg(:)));
end
globalmax = max(localmax);
globalmin = min(localmin);

for st = 1: length(statlist)
    stat = statlist{st}.stat;
    % Expand mask array to cover same time as amplitude envelope-filtered EEG
    Tdiffintrsect = ismember(timeavg_AudSampList{1}.time,stat.time);
    mask = zeros(length(stat.label),length(Tdiffintrsect)); mcntr=1;
    for i = 1:length(Tdiffintrsect),
        if Tdiffintrsect(i), mask(:,i) = stat.mask(:,mcntr); mcntr=mcntr+1;  end
    end

    fsample = 1/mode(diff(timeavg_AudSampList{1, 1}.time)); %Hz
    tslength = length(timeavg_AudSampList{1}.time);
    downsampfactor = 20;
    Tstart = find(round(timeavg_AudSampList{1}.time,3)==0.252); 
    Tend   = find(round(timeavg_AudSampList{1}.time,3)==0.752);
    TSindices = round((Tstart:(fsample/downsampfactor):Tend)); 
    time_segment = timeavg_AudSampList{1}.time(TSindices);



    % generate the subplots
    myfig=figure;
    hndl = tight_subplot(length(statlist{st}.condidx),length(TSindices),[.001 .02],[.05 .02],[.01 .08]);
    SPcntr = 0;
    for audsamp = statlist{st}.condidx
        for ts =1:length(TSindices)
            SPcntr = SPcntr + 1;
            axes(hndl(SPcntr));

            datavector = timeavg_AudSampList{audsamp}.avg(chnIDCS,TSindices(ts));
            datamask = mask(chnIDCS,TSindices(ts));
            topoplot(datavector, chanLocs,'maplimits',[globalmin globalmax],'conv','on',...
                    'plotrad',0.65,'pmask',datamask);
            colormap((brewermap(100,'YlGnBu')));    
            caxis([globalmin globalmax])

            % specific labels and colorbars on places on the figure
            if ts == length(TSindices) && audsamp == statlist{st}.condidx(end),
                cbh = colorbar('location','east');
                set(cbh,'position',[0.947 0.147 0.016 0.747],'fontweight','bold','fontsize',10)
                ylabel(cbh,'EEG Amplitude (\muV)','FontWeight','bold','FontSize',12,'Position',[-1.173 149 0]) 
            end
            if audsamp == 1, title(round(time_segment(ts),2),'fontsize',14,'fontweight','bold'); end

        end
    end
    set(gcf,'paperpositionmode','auto','position',[1 1 1780 323])
    %% Save the figures
    filenamesave = ['Imported Matlab Data\ThetaSignalEnvelope_',strrep(strjoin(statlist{st}.condlabel,'vs'),'-',''),'_TopoplotFigure.png'];
    print(myfig,'-dpng', filenamesave, '-r600')
end

%%
    cfg = [];
    cfg.xlim = [-2:0.5:2];
%     cfg.zlim = [min(diff_TL.avg(:)) max(diff_TL.avg(:))];  % across whole epoch
    cfg.zlim = [min(min(timeavg_AudSampList{1}.avg(:))) max(max(timeavg_AudSampList{1}.avg(:)))]; %for only specified segment 
    cfg.marker = 'labels';
    cfg.layout = lay;
    cfg.comment = 'xlim';
    cfg.commentpos = 'title';
    cfg.colorbar = 'yes';
    cfg.colormap = flipud(brewermap(100,'RdBu'));
    cfg.fontsize = 14;
    cfgtopo = ft_topoplotER(cfg, timeavg_AudSampList{1});
