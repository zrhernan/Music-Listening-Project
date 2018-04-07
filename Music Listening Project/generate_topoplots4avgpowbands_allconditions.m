clearvars; close all; clc
load('relpowbands_ASRdenoisedEEG_onlylistenphase.mat')
Conds = {'selfselected','bach','gagaku','clicking','cronkite','chaplin'};
%%
powband_names = fieldnames(relpowbands);
FB_size = length(powband_names);
numchns = 15;       numconds = length(Conds);
powbands_perCond = cell(FB_size,numconds);
avgpow_perCond = cell(FB_size,numconds);
for frqbnd = 1:FB_size
    for cond = [1 5 6 2 3 4];%1:numconds
        for subj = [1:4,6:12]%1:12
            if eval(['sum(relpowbands.',powband_names{frqbnd},'(subj,cond,:),3) == 0']), continue, end
            
            eval(['powbands_perCond{frqbnd,cond} = [powbands_perCond{frqbnd,cond}; transpose(squeeze(relpowbands.',...
                powband_names{frqbnd},'(subj,cond,:)))];'])
        end        
        avgpow_perCond{frqbnd,cond} = nanmean(powbands_perCond{frqbnd,cond},1);
    end
    
end

%% display familiar vs. unfamiliar average power bands on topoplots
load('BrainVision_1020_16ChanLocs.mat')
chanLocs(ismember({chanLocs.labels},'O1')) = []; %remove O1 channel
globmax = max([avgpow_perCond{1:end-1,:}]);
globmin = min([avgpow_perCond{1:end-1,:}]);
fbcntr = 0;
displayNames = {'Delta Band (1-4 Hz)','Theta Band (4-7 Hz)',...
    'Low Alpha Band (7-12 Hz)','High Alpha Band (12-15 Hz)',...
    'Low Beta Band (15-24 Hz)','High Beta Band (24-30 Hz)',...
    'Low Gamma Band (30-50 Hz)',};
freqband = [1 4; 4 7; 7 12; 12 15; 15 24; 24 30; 30 50];
hndl = tight_subplot(numconds,FB_size-1,[.001 .02],[.05 .03],[.01 .01]);
CBxpos = [0.045 0.188 0.331 0.474 0.617 0.76 0.902];

for fb = 1:(FB_size-1)
    fbcntr = fbcntr + 1;
    axes(hndl(fb));

    title(displayNames{fbcntr})
    
    localmax = max([avgpow_perCond{fbcntr,:}]);
    localmin = min([avgpow_perCond{fbcntr,:}]);

    datavector = avgpow_perCond{fbcntr,1};
    topoplot(datavector, chanLocs,'maplimits',[localmin localmax],'conv','on',...
        'plotrad',0.65);
    colormap(gray(100))
    caxis([localmin localmax])
    
    for cond = 2:numconds
        axes(hndl(fb+7*(cond-1)));   
        datavector2 = avgpow_perCond{fbcntr,cond};

        topoplot(datavector2, chanLocs,'maplimits',[localmin localmax],'conv','on',...
            'plotrad',0.65);
        caxis([localmin localmax])
        colormap(jet(100))
        if cond == 6
            cb = colorbar('Location','SouthOutside');
            set(cb,'position',[CBxpos(fb) 0.043 0.053 0.012])
            set(cb,'ticks',[ceil(100*localmin)/100;floor(100*localmax)/100])
        end
        if fb == 4 && cond == 6, ylabel(cb,'Average Relative Power (%)','fontsize',12), end
    end
 
end
myfig=gcf;
set(myfig,'position',[1 1 1364 909])

%% Save it!
cd('C:\Users\TMHZXH5\Documents')
filename = 'TopographicPlots_AverageRelativePower_AllConditions.tif';
set(gcf,'paperpositionmode','auto')
print(myfig,'-dpng','-r300',filename)
