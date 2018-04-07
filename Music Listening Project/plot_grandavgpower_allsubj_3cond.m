%% Load grand-averaged power list
load('grand_avg_pow_ASRdenoisedEEG.mat')

%% Plot the average power over multiple freq bands, subjects, and conditions
myfig = figure(5);
frqbnd_labels = {'1-4';'4-7';'7-12'; '12-15'; '15-24'; '24 30'; '30-50'};
avg_pow_cat=cat(1,avg_pow(:,:,1),avg_pow(:,:,2),avg_pow(:,:,3));
imagesc(avg_pow_cat)
set(gca,'ytick',(1:36),'yticklabel',(1:12),'xtick', (1:7),'xticklabel',...
    frqbnd_labels)
set(gca,'fontsize',12,'fontweight','bold')

xlabel('Frequency Bands (Hz)','fontweight','bold')
ylabel(['Music Listening Condition',10,...
'  Gagaku                                     Bach                                 Self-Selected'],...
'fontweight','bold')

caxis([0 1]), cbh = colorbar;


% horizontal lines
linlength = (0.5:0.5:7.5);
line(linlength,12.5*ones(length(linlength)),'linewidth',2,'color','k')
line(linlength,24.5*ones(length(linlength)),'linewidth',2,'color','k')

set(gcf,'position',[1 1 578 838])
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
