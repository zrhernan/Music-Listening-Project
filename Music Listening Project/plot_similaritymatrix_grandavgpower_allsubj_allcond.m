%% Load grand-averaged power list
load('grand_avg_pow_ASRdenoisedEEG.mat')

%% Generate Similarity Matrix
myfig = figure(6);
AP_cat = cat(1,avg_pow(:,:,1),avg_pow(:,:,2),avg_pow(:,:,3));
smlr_vect = pdist(AP_cat,'correlation');
smlr_arry = squareform(smlr_vect);
subjNo = [(1:12),(1:12),(1:12),(1:4),(6:12),(2:4),(6:8),11,(2:4),(6:8),11];

imagesc(smlr_arry); axis image
set(gca,'ytick',(1:length(subjNo)),'yticklabel',subjNo,'xtick', (1:length(subjNo)),'xticklabel',...
    subjNo)
set(gca,'fontsize',12,'fontweight','bold')
xlabel(['Self-Selected                                Bach                                   Gagaku     ',...
10,'Music Listening Condition'],...
'fontweight','bold')
ylabel(['Music Listening Condition',10,...
'  Gagaku                                     Bach                                 Self-Selected'],...
'fontweight','bold')


caxis([0 1]), cbh = colorbar;
ylabel(cbh,'Correlation Coefficient','rotation',-90,...
    'verticalalignment','bottom','fontsize',14,'fontweight','bold')

% horizontal lines
linlength = (0.5:0.5:36.5);
line(linlength,12.5*ones(length(linlength)),'linewidth',2,'color','w')
line(linlength,24.5*ones(length(linlength)),'linewidth',2,'color','w')

% vertical lines
line(12.5*ones(length(linlength)),linlength,'linewidth',2,'color','w')
line(24.5*ones(length(linlength)),linlength,'linewidth',2,'color','w')

set(gcf,'position',[1 1 1126 865])%1126 865

%% Save it!
cd('C:\Users\TMHZXH5\Documents')
filename = 'SimilarityMatrix_allsubj_allMusicCond_7avgpows_ASRdenoised.tif';
set(gcf,'paperpositionmode','auto')
print(myfig,'-dpng','-r600',filename)
