%% Load extracted features list
cd('Z:\Research_Projects\MUSIC_Karmonik\')
% load('extractfeats_ASRdenoisedEEG.mat')
load('relpowbands_ASRdenoisedEEG_onlylistenphase.mat')
extractedfeats = relpowbands;
%% Concatenate subject-condition dimensions together for all features
feat_names = fields(extractedfeats);
feat_size = length(feat_names);
MusicCondLabels = {'Self-Selected','Bach','Gagaku','Clicking Langauge','Cronkite','Chaplin'};
subjNo = [(1:12),(1:12),(1:12),(1:4),(6:12),(2:4),(6:8),11,(2:4),(6:8),11];


% for 
    ftr = 1;%9:15%:feat_size
    close all
    allfeat_cat = [];
    eval(['allfeat_cat = [allfeat_cat, squeeze(cat(1,extractedfeats.',feat_names{ftr},...
        '(:,1,:),extractedfeats.',feat_names{ftr},...
        '(:,2,:),extractedfeats.',feat_names{ftr},...
        '(:,3,:),extractedfeats.',feat_names{ftr},...
        '(:,4,:),extractedfeats.',feat_names{ftr},...
        '(:,5,:),extractedfeats.',feat_names{ftr},'(:,6,:)))];']);


allfeat_cat(sum(allfeat_cat,2)==0,:)=[];
smlr_arry = corrcoef(allfeat_cat);

%% Generate Similarity Matrix based on Features
% close all
myfig = figure(555);
imagesc(abs(smlr_arry)); axis image

set(gca,'xtick',[6.5 18.5 30.5 42 51 58],'xticklabel',MusicCondLabels,'ytick',...
    [6.5 18.5 30.5 42 51 58],'yticklabel',MusicCondLabels)
set(gca,'fontsize',12,'fontweight','bold')
xlabel('Music Listening Condition','fontweight','bold')
ylabel('Music Listening Condition','fontweight','bold')


caxis([0 1]), cbh = colorbar;
ylabel(cbh,'Absolute Correlation Coefficient','rotation',-90,...
    'verticalalignment','bottom','fontsize',14,'fontweight','bold')

% horizontal lines
linlength = (0.5:0.5:61.5);
line(linlength,12.5*ones(length(linlength)),'linewidth',2,'color','k')
line(linlength,24.5*ones(length(linlength)),'linewidth',2,'color','k')
line(linlength,36.5*ones(length(linlength)),'linewidth',2,'color','k')
line(linlength,47.5*ones(length(linlength)),'linewidth',2,'color','k')
line(linlength,54.5*ones(length(linlength)),'linewidth',2,'color','k')

% vertical lines
line(12.5*ones(length(linlength)),linlength,'linewidth',2,'color','k')
line(24.5*ones(length(linlength)),linlength,'linewidth',2,'color','k')
line(36.5*ones(length(linlength)),linlength,'linewidth',2,'color','k')
line(47.5*ones(length(linlength)),linlength,'linewidth',2,'color','k')
line(54.5*ones(length(linlength)),linlength,'linewidth',2,'color','k')

set(gcf,'position',[1 1 1727 1727])

%% Save it!
cd('C:\Users\TMHZXH5\Documents')
filename = ['SimilarityMatrix_arrangedXcond_allsubj_allMusicCond_',feat_names{ftr},'bandfeats_ASRdenoised.tif'];
set(gcf,'paperpositionmode','auto')
print(myfig,'-dpng','-r300',filename)

% end