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
cd('C:\Users\zrhernan\Documents\MATLAB\Listening Project Matlab Code')
load('EEG Montage Templates\BrainVision_1020_16ChanLocs.mat')
load('Imported Matlab Data\ERPdata_GroupedbyAuditorySample.mat')
fsample = 250; %Hz
tslength = length(timeavg_AudSampList{1}.time);
TSindices = ((tslength-1)*0.5+1 : fsample/5 : (tslength-1)*0.75+1); 
time_segment = timeavg_AudSampList{1}.time(TSindices);
myfig=figure;
hndl = tight_subplot(length(timeavg_AudSampList),length(TSindices),[.001 .02],[.05 .02],[.01 .08]);
SPcntr = 0;

% get global and local (to each auditory sample) max and min
localmax = nan(length(timeavg_AudSampList),1);      localmin = nan(size(localmax));
for audsamp = 1:length(timeavg_AudSampList)
    localmax(audsamp) = max(max(timeavg_AudSampList{audsamp}.avg(:)));
    localmin(audsamp) = min(min(timeavg_AudSampList{audsamp}.avg(:)));
end
globalmax = max(localmax);
globalmin = min(localmin);

% generate the subplots
for audsamp = 1:length(timeavg_AudSampList)
    for ts =1:length(TSindices)
        SPcntr = SPcntr + 1;
        axes(hndl(SPcntr));
        
        datavector = timeavg_AudSampList{audsamp}.avg(:,TSindices(ts));
        topoplot(datavector, chanLocs,'maplimits',[globalmin globalmax],'conv','on',...
                'plotrad',0.65);
        colormap((brewermap(100,'YlGnBu')));    
        caxis([globalmin globalmax])
             
        % specific labels and colorbars on places on the figure
        if ts == length(TSindices) && audsamp == length(timeavg_AudSampList),
            cbh = colorbar('location','east');
            set(cbh,'position',[0.947 0.324 0.025 0.356],'fontweight','bold','fontsize',10)
            ylabel(cbh,'EEG Amplitude (\muV)','FontWeight','bold','FontSize',12,'Position',[-1.173 149 0]) 
        end
        if audsamp == 1, title(time_segment(ts),'fontsize',14,'fontweight','bold'); end
            
    end
end

set(gcf,'position',[1 1 827 810])
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
