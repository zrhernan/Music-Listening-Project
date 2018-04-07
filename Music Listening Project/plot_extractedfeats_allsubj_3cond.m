%% Load extracted features list
cd('Z:\Research_Projects\MUSIC_Karmonik\')
load('extractfeats_ASRdenoisedEEG.mat')

%% Plot the statistical features over channels, subjects, and conditions
feat_names = fields(extractedfeats);
feat_size = length(feat_names);
featNameLabels = {'Maximum';'Minimum';'Variance';'Skewness';...
    'Kurtosis';'Interquartile Range';'Median Absolute Deviation';'Shannon Entropy';'Mean Power (1-4 Hz)';...
    'Mean Power (4-7 Hz)';'Mean Power (7-12 Hz)';'Mean Power (12-15 Hz)';...
    'Mean Power (15-24 Hz)';'Mean Power (24-30 Hz)';'Mean Power (30-50 Hz)';'Mean Power (1-100 Hz)'};
for ftr = 1:feat_size
    
    myfig = figure(ftr+100);
    eval(['feat_cat=squeeze(cat(1,extractedfeats.',feat_names{ftr},...
        '(:,1,:),extractedfeats.',feat_names{ftr},...
        '(:,2,:),extractedfeats.',feat_names{ftr},'(:,3,:)));']);
    imagesc(feat_cat)
    set(gca,'ytick',(1:36),'yticklabel',(1:12),'xtick', (1:15),'xticklabel',...
        data_preproc.label)
    set(gca,'fontsize',12,'fontweight','bold')

    xlabel('EEG Channels','fontweight','bold')
    ylabel(['Music Listening Condition',10,...
    '  Gagaku                                     Bach                                 Self-Selected'],...
    'fontweight','bold')

    caxis([min(feat_cat(:)) max(feat_cat(:))]), cbh = colorbar;


    % horizontal lines
    linlength = (0.5:0.5:15.5);
    line(linlength,12.5*ones(length(linlength)),'linewidth',2,'color','k')
    line(linlength,24.5*ones(length(linlength)),'linewidth',2,'color','k')

    set(gcf,'position',[1 1 678 838])
    % ylabel(cbh,'Average Power Amplitude (A.U.)','FontWeight','bold',...
    %     'verticalalignment','bottom','rotation',-90)
    cbh.Label.String = featNameLabels{ftr};
    cbh.Label.FontSize = 14;
    cbh.Label.Rotation = -90;
    cbh.Label.VerticalAlignment = 'bottom';
    
%% Save it!
cd('C:\Users\TMHZXH5\Documents')
filename = ['ExtractedFeatures_',feat_names{ftr},'_allsubj_3MusicCond_ASRdenoised.tif'];
set(gcf,'paperpositionmode','auto')
print(myfig,'-dpng','-r600',filename)
end