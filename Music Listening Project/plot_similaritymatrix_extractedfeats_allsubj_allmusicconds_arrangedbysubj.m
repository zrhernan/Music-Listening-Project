%% Load extracted features list
cd('Z:\Research_Projects\MUSIC_Karmonik\')
% load('extractfeats_ASRdenoisedEEG.mat')
load('relpowbands_ASRdenoisedEEG_onlylistenphase.mat')
extractedfeats = relpowbands;
%% Generate Similarity Matrix based on Features
close all
myfig = figure(555);
feat_names = fields(extractedfeats);
feat_size = length(feat_names);
subjNo = [(1:12),(1:12),(1:12),(1:4),(6:12),(2:4),(6:8),11,(2:4),(6:8),11];
subjNoidcs = {find(subjNo==1),find(subjNo==2),find(subjNo==3),...
    find(subjNo==4),find(subjNo==5),find(subjNo==6),find(subjNo==7),...
    find(subjNo==8),find(subjNo==9),find(subjNo==10),find(subjNo==11),find(subjNo==12)};
allfeat_cat = [];

for ftr = 1%1:feat_size
    eval(['allfeat_cat = [allfeat_cat, squeeze(cat(1,extractedfeats.',feat_names{ftr},...
        '(:,1,:),extractedfeats.',feat_names{ftr},...
        '(:,2,:),extractedfeats.',feat_names{ftr},...
        '(:,3,:),extractedfeats.',feat_names{ftr},...
        '(:,4,:),extractedfeats.',feat_names{ftr},...
        '(:,5,:),extractedfeats.',feat_names{ftr},'(:,6,:)))];']);
end

allfeat_cat(sum(allfeat_cat,2)==0,:)=[];
allfeat_arrangeXsubj = [];
for sbjt = 1:length(subjNoidcs)
    allfeat_arrangeXsubj = [allfeat_arrangeXsubj; allfeat_cat(subjNoidcs{sbjt},:)];
end
smlr_arry = corrcoef(allfeat_arrangeXsubj');
imagesc(smlr_arry); axis image

set(gca,'ytick',[2.5 7.5 13.5 19.5 24 28.5 34.5 40.5 45.5 49.5 54.5 59.5],'yticklabel',(1:12),...
    'xtick',[2.5 7.5 13.5 19.5 24 28.5 34.5 40.5 45.5 49.5 54.5 59.5],'xticklabel',(1:12))
set(gca,'fontsize',12,'fontweight','bold')
xlabel('Subjects','fontweight','bold')
ylabel('Subjects','fontweight','bold')


caxis([0 1]), cbh = colorbar;
ylabel(cbh,'Correlation Coefficient','rotation',-90,...
    'verticalalignment','bottom','fontsize',14,'fontweight','bold')

% horizontal lines
linlength = (0.5:0.5:61.5);
line(linlength,4.5*ones(length(linlength)),'linewidth',2,'color','k')
line(linlength,10.5*ones(length(linlength)),'linewidth',2,'color','k')
line(linlength,16.5*ones(length(linlength)),'linewidth',2,'color','k')
line(linlength,22.5*ones(length(linlength)),'linewidth',2,'color','k')
line(linlength,25.5*ones(length(linlength)),'linewidth',2,'color','k')
line(linlength,31.5*ones(length(linlength)),'linewidth',2,'color','k')
line(linlength,37.5*ones(length(linlength)),'linewidth',2,'color','k')
line(linlength,43.5*ones(length(linlength)),'linewidth',2,'color','k')
line(linlength,47.5*ones(length(linlength)),'linewidth',2,'color','k')
line(linlength,51.5*ones(length(linlength)),'linewidth',2,'color','k')
line(linlength,57.5*ones(length(linlength)),'linewidth',2,'color','k')

% vertical lines
line(4.5*ones(length(linlength)),linlength,'linewidth',2,'color','k')
line(10.5*ones(length(linlength)),linlength,'linewidth',2,'color','k')
line(16.5*ones(length(linlength)),linlength,'linewidth',2,'color','k')
line(22.5*ones(length(linlength)),linlength,'linewidth',2,'color','k')
line(25.5*ones(length(linlength)),linlength,'linewidth',2,'color','k')
line(31.5*ones(length(linlength)),linlength,'linewidth',2,'color','k')
line(37.5*ones(length(linlength)),linlength,'linewidth',2,'color','k')
line(43.5*ones(length(linlength)),linlength,'linewidth',2,'color','k')
line(47.5*ones(length(linlength)),linlength,'linewidth',2,'color','k')
line(51.5*ones(length(linlength)),linlength,'linewidth',2,'color','k')
line(57.5*ones(length(linlength)),linlength,'linewidth',2,'color','k')

set(gcf,'position',[1 1 1727 1727])

%% Save it!
cd('C:\Users\TMHZXH5\Documents')
filename = 'SimilarityMatrix_arrangedXsubj_allsubj_allMusicCond_relpowlogammabandfeats_ASRdenoised.tif';
set(gcf,'paperpositionmode','auto')
print(myfig,'-dpng','-r300',filename)