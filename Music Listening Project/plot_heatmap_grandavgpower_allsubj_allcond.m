%% Load grand-averaged power list
load('grand_avg_pow_ASRdenoisedEEG.mat')

%% Plot the average power over multiple freq bands, subjects, and conditions
myfig = figure(5);
frqbnd_labels = {'1-4';'4-7';'7-12'; '12-15'; '15-24'; '24 30'; '30-50'};
avg_pow_cat=cat(1,avg_pow(:,:,1),avg_pow(:,:,2),avg_pow(:,:,3),avg_pow(:,:,4),avg_pow(:,:,5),avg_pow(:,:,6));
avg_pow_cat(sum(avg_pow_cat,2)==0,:)=[];

subjNo = [(1:12),(1:12),(1:12),(1:4),(6:12),(2:4),(6:8),11,(2:4),(6:8),11];

imagesc(avg_pow_cat')
set(gca,'xtick',(1:length(subjNo)),'xticklabel',subjNo,'ytick', (1:7),'yticklabel',...
    frqbnd_labels)
set(gca,'fontsize',12,'fontweight','bold')


xlabel(['                     Self-Selected                               Bach                                        Gagaku                             Click Language                 Cronkite                 Chaplin              ',...
10,'Music Listening Condition'],'fontweight','bold')
ylabel('Frequency Bands (Hz)','fontweight','bold')


caxis([0 0.5]), cbh = colorbar;


% vertical lines
linlength = (0.5:0.5:7.5);
line(12.5*ones(length(linlength)),linlength,'linewidth',2,'color','w')
line(24.5*ones(length(linlength)),linlength,'linewidth',2,'color','w')
line(36.5*ones(length(linlength)),linlength,'linewidth',2,'color','w')
line(47.5*ones(length(linlength)),linlength,'linewidth',2,'color','w')
line(54.5*ones(length(linlength)),linlength,'linewidth',2,'color','w')

set(gcf,'position',[1 1 1727 449])
% ylabel(cbh,'Average Power Amplitude (A.U.)','FontWeight','bold',...
%     'verticalalignment','bottom','rotation',-90)
cbh.Label.String = 'Spectral Amplitude Power  (A.U.)';
cbh.Label.FontSize = 14;
cbh.Label.Rotation = -90;
cbh.Label.VerticalAlignment = 'bottom';

%% save it
cd('C:\Users\TMHZXH5\Documents')
filename = 'FreqBandAveragePower_allsubj_3MusicCond_ASRdenoised.tif';
set(gcf,'paperpositionmode','auto')
print(myfig,'-dpng','-r600',filename)
