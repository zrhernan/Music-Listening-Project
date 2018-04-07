clearvars,clc
load 'ASRdenoisedEEG.mat'
freqband = [1 4; 4 7; 7 12; 12 15; 15 24; 24 30; 30 50];
glob_rel_pow = zeros(12,size(freqband,1),6);
Conds = {'selfselected','bach','gagaku','clicking','cronkite','chaplin'};
%%
datafft_CondGroup1 = cell(11,1);     subjcnt = 0;
subjectselection = (1:12);%[1:4,6:12];
for subj = subjectselection
    subjcnt = subjcnt + 1;
    condspersubj = cell(6,1);
    nofftgrp1flag = 0;      %nofftgrp2flag = 0;
    for cond = 2%1:6
        if isempty(preprocdata_allsubjANDcond{subj,cond}), continue, end
        data_preproc = preprocdata_allsubjANDcond{subj,cond};
       
        %% Define event trials
        % focusing on 12 seconds before and after the listening onset (no onsets)
        r2lonset = (1:6000:length(data_preproc.time{1}))'; %rest-to-listen onsets
        trl = [r2lonset, r2lonset+5999, 0*ones(length(r2lonset),1)];
        trl(1:2:length(r2lonset),4) = ones(length(1:2:length(r2lonset)),1);
        trl(2:2:length(r2lonset),4) = 2*ones(length(2:2:length(r2lonset)),1);
        if trl(end,2) > length(data_preproc.time{1}), trl(end,2) = length(data_preproc.time{1}); end
        cfg           = [];
        cfg.trl       = trl;
        datatrl       = ft_redefinetrial(cfg,data_preproc);
        condspersubj{cond} = datatrl;
    end
%% append each condition group per subject (all conditions)
    grp1_idcs = [1:6];%[6];
%     grp2_idcs = [2];%[4];
    grp1_idcskept = grp1_idcs(~cellfun('isempty',{condspersubj{grp1_idcs}})); % music conditions
%     grp2_idcskept = grp2_idcs(~cellfun('isempty',{condspersubj{grp2_idcs}})); % language conditions
    cfg          = [];
    numGrp1Condskept = length({condspersubj{grp1_idcskept}});
%     numGrp2Condskept = length({condspersubj{grp2_idcskept}});
    if numGrp1Condskept > 1
        CondGroup1perSubj = ft_appenddata(cfg,condspersubj{grp1_idcskept});
    elseif numGrp1Condskept == 1 
        CondGroup1perSubj = condspersubj{grp1_idcskept};
    elseif numGrp1Condskept == 0
        nofftgrp1flag = 1;
    end
%     if numGrp2Condskept > 1
%         CondGroup2perSubj = ft_appenddata(cfg,condspersubj{grp2_idcskept});
%     elseif numGrp2Condskept == 1 
%         CondGroup2perSubj = condspersubj{grp2_idcskept};
%     elseif numGrp2Condskept == 0
%         nofftgrp2flag = 1;
%     end

%% calculate Fourier spectra for each condition (listen block only)
    cfg                 = [];
    cfg.trials          = find(datatrl.trialinfo==2);
    cfg.method          = 'mtmfft';
    cfg.output          = 'pow';
    cfg.pad             = 'nextpow2';
    cfg.foi             = (0.1:0.1:100);
    cfg.taper           = 'hanning';
    cfg.tapsmofrq       = 2;

    if nofftgrp1flag ~= 1
        datafft_CondGroup1{subjcnt} = ft_freqanalysis(cfg, CondGroup1perSubj);
    end

end

%% calculate grand average across subjects for each condition group   
cfg                 = [];
cfg.channel         = 'all';
cfg.foilim          = 'all';
cfg.parameter       = 'powspctrm';
cfg.keepindividual  = 'yes';
FGACondGroup1        = ft_freqgrandaverage(cfg,datafft_CondGroup1{~cellfun('isempty',{datafft_CondGroup1{:}})}); 
% "{:}" means to use data from all elements of the variable




%% plot results (per subject)
close all
load easycap16HMRImusicproj_layout
for s=1:length(subjectselection)
close all
myfig = figure;
cfgp            = [];
cfgp.xlim       = [datafft_CondGroup1{s}.freq(1) datafft_CondGroup1{s}.freq(end)];
cfgp.ylim       = [fix(min(datafft_CondGroup1{s}.powspctrm(:))) 3];%max(datafft_CondGroup1{s}.powspctrm(:))]; %
cfgp.layout     = lay;
cfgp.showlabels = 'yes';
cfgp.interactive= 'yes';
ft_multiplotER(cfgp, datafft_CondGroup1{s});


% Save the files!
cd('C:\Users\TMHZXH5\Documents\Fourier Spectra Figures')
pause(0.00001);
warning('off','verbose')
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);

filename = ['FourierSpectra_AllConditions_Subject',num2str(subjectselection(s)),'.tif'];
set(gcf,'paperpositionmode','auto')
print(myfig,'-dpng','-r300',filename)
end

%% plot results (all subjects)
load easycap16HMRImusicproj_layout
myfig = figure;
cfgp            = [];
cfgp.xlim       = [FGACondGroup1.freq(1) FGACondGroup1.freq(end)];
cfgp.ylim       = [fix(min(FGACondGroup1.powspctrm(:)))  3];%max(FGACondGroup1.powspctrm(:))]; %
cfgp.layout     = lay;
cfgp.showlabels = 'yes';
cfgp.interactive= 'yes';
ft_multiplotER(cfgp, FGACondGroup1);

%% Save the figure
cd('C:\Users\TMHZXH5\Documents\Fourier Spectra Figures')
filename = 'FourierSpectra_AllConditions_AllSubjects.tif';
set(gcf,'paperpositionmode','auto')
print(myfig,'-dpng','-r300',filename)
