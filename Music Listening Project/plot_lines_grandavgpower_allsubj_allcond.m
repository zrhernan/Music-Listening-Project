%% Load grand-averaged power list
load('global_rel_pow_ASRdenoisedEEG.mat')

%% Concatenate Subject and Condition Dimensions together
frqbnd_labels = {'1-4';'4-7';'7-12'; '12-15'; '15-24'; '24-30'; '30-50'};
MusicCondLabels = {'Self-Selected','Bach','Gagaku','Clicking Langauge','Cronkite','Chaplin'};
rel_pow_cat=cat(1,glob_rel_pow(:,:,1),glob_rel_pow(:,:,2),glob_rel_pow(:,:,3),...
    glob_rel_pow(:,:,4),glob_rel_pow(:,:,5),glob_rel_pow(:,:,6));
rel_pow_cat(sum(rel_pow_cat,2)==0,:)=[];

grand_avg_rel_pow = [mean(rel_pow_cat(1:12,:));mean(rel_pow_cat(13:24,:));...
                     mean(rel_pow_cat(25:36,:));mean(rel_pow_cat(37:47,:));...
                     mean(rel_pow_cat(48:54,:));mean(rel_pow_cat(55:61,:))];
grand_std_rel_pow = [std(rel_pow_cat(1:12,:));std(rel_pow_cat(13:24,:));...
                     std(rel_pow_cat(25:36,:));std(rel_pow_cat(37:47,:));...
                     std(rel_pow_cat(48:54,:));std(rel_pow_cat(55:61,:))];
subjidcs = {(1:12),(13:24),(25:36),(37:47),(48:54),(55:61)};
%% Make boxplots
for c=1:6
    subplot(3,2,c)
    
    boxplot(rel_pow_cat(subjidcs{c},:))
    title(MusicCondLabels{c})
    ylim([0 0.7])
    set(gca,'xtick',(1:7),'xticklabel',frqbnd_labels,'ytick',(0:0.1:0.7))
    set(gca,'fontsize',12,'fontweight','bold','yminortick','on','yminorgrid','on')
    if c==6, xlabel('Frequency Bands (Hz)','fontweight','bold'), end
    if c==3, ylabel('Global Relative Power Across Subjects per Condition','fontweight','bold'), end

end

%% Plot the average power over multiple freq bands, subjects, and conditions
myfig = figure(5);
% % display as points on a graph
% colorlist = {'r*','b*','k*','m*','y*','c*',};
% for fq = 1:length(MusicCondLabels)
%     plot((1:7),subjs_grand_avg_pow(fq,:),colorlist{fq})
%     hold on
% end

% display as bar plot
% bar(grand_avg_rel_pow','grouped')
cmap = parula(10);
errorbar_groups(grand_avg_rel_pow,grand_std_rel_pow,...
    'bar_colors',cmap(4:end,:),'bar_names',frqbnd_labels,'FigID',5);
grid on
% set(gca,'xtick', (1:7),'xticklabel',frqbnd_labels,'ytick',(0:0.1:0.5))
set(gca,'fontsize',12,'fontweight','bold','yminortick','on','yminorgrid','on')

xlabel('Frequency Bands (Hz)','fontweight','bold')
ylabel(['Global Relative Power Averaged',10,'Across Subjects per Condition'],'fontweight','bold')


legend(MusicCondLabels)

% set(gcf,'position',[1 1 1727 449])


%% save it
cd('C:\Users\TMHZXH5\Documents')
filename = 'FreqBandGrandAverageGlobalPower_allsubj_allMusicCond_ASRdenoised_onlylistenphase.tif';
set(gcf,'paperpositionmode','auto')
print(myfig,'-dpng','-r300',filename)
