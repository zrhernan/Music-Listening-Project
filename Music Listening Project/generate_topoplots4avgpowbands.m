clearvars; close all; clc
load('relpowbands_ASRdenoisedEEG_onlylistenphase.mat')
Conds = {'selfselected','bach','gagaku','clicking','cronkite','chaplin'};
%%
powband_names = fieldnames(relpowbands);
FB_size = length(powband_names);
numchns = 15;
powbands_CondGroup1 = cell(FB_size,1);
powbands_CondGroup2 = cell(FB_size,1);
avgpow_CondGroup1 = nan(FB_size,numchns);
avgpow_CondGroup2 = nan(FB_size,numchns);
for frqbnd = 1:FB_size
    for subj = 1:12
        for cond = [1 5]  % group familiar listening conditions
            if eval(['sum(relpowbands.',powband_names{frqbnd},'(subj,cond,:),3) == 0']), continue, end
            eval(['powbands_CondGroup1{frqbnd} = [powbands_CondGroup1{frqbnd}; transpose(squeeze(relpowbands.',...
                powband_names{frqbnd},'(subj,cond,:)))];']) 
        end
        
        for cond = [3 4];  % group unfamiliar listening conditions
            if eval(['sum(relpowbands.',powband_names{frqbnd},'(subj,cond,:),3) == 0']), continue, end
            eval(['powbands_CondGroup2{frqbnd} = [powbands_CondGroup2{frqbnd}; transpose(squeeze(relpowbands.',...
                powband_names{frqbnd},'(subj,cond,:)))];']) 
        end        
        
    end
    avgpow_CondGroup1(frqbnd,:) = nanmean(powbands_CondGroup1{frqbnd},1);
    avgpow_CondGroup2(frqbnd,:) = nanmean(powbands_CondGroup2{frqbnd},1);
end

%% display familiar vs. unfamiliar average power bands on topoplots
load('BrainVision_1020_16ChanLocs.mat')
chanLocs(ismember({chanLocs.labels},'O1')) = []; %remove O1 channel
globmax = max(max(avgpow_CondGroup1(1:end-1,:)));
globmin = min(min(avgpow_CondGroup1(1:end-1,:)));
fbcntr = 0;
displayNames = {'Delta Band (1-4 Hz)','Theta Band (4-7 Hz)',...
    'Low Alpha Band (7-12 Hz)','High Alpha Band (12-15 Hz)',...
    'Low Beta Band (15-24 Hz)','High Beta Band (24-30 Hz)',...
    'Low Gamma Band (30-50 Hz)',};
freqband = [1 4; 4 7; 7 12; 12 15; 15 24; 24 30; 30 50];
ha = tight_subplot(2,FB_size-1,[.001 .02],[.1 .01],[.01 .01]);

for plotnum = 1:(FB_size-1)
    fbcntr = fbcntr + 1;
    axes(ha(plotnum));
%     subplot(2,FB_size-1,plotnum)
    title(displayNames{fbcntr})
%     if plotnum == 1, ylabel(gca,'Familiar Listening'), end
    datavector = avgpow_CondGroup1(fbcntr,:);
    localmax = max(max(avgpow_CondGroup1(fbcntr,:)), max(avgpow_CondGroup2(fbcntr,:)));
    localmin = min(min(avgpow_CondGroup1(fbcntr,:)), min(avgpow_CondGroup2(fbcntr,:)));

    topoplot(datavector, chanLocs,'maplimits',[localmin localmax],'conv','on',...
        'plotrad',0.65);
    colormap(gray(100))
    colorbar('Location','SouthOutside')    
    axes(ha(plotnum+7));
%     subplot(2,FB_size-1,plotnum+7)
%     if plotnum == 1, ylabel('Unfamiliar Listening'), end
    datavector2 = avgpow_CondGroup2(fbcntr,:);
    
    topoplot(datavector2, chanLocs,'maplimits',[localmin localmax],'conv','on',...
        'plotrad',0.65);
    colormap(gray(100))
    cb = colorbar('Location','SouthOutside');
    if plotnum == 4, ylabel(cb,'Average Relative Power (%)','fontsize',12), end
    
 
end
myfig=gcf;
set(myfig,'position',[1 1 1899 728])

%% Save it!
cd('C:\Users\TMHZXH5\Documents')
filename = 'TopographicPlots_AverageRelativePower_FamiliarvsUnfamiliar.tif';%MusicvsLangauge.tif';
set(gcf,'paperpositionmode','auto')
print(myfig,'-dpng','-r300',filename)
