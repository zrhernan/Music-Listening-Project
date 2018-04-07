%% Load grand-averaged power list
load('global_rel_pow_ASRdenoisedEEG.mat')

%% Initializations
frqbnd_labels = {'1-4';'4-7';'7-12'; '12-15'; '15-24'; '24-30'; '30-50'};
MusicCondLabels = {'Self-Selected','Bach','Gagaku','Clicking Langauge','Cronkite','Chaplin'};
subjNo = [(1:12),(1:12),(1:12),(1:4),(6:12),(2:4),(6:8),11,(2:4),(6:8),11];
condNo = [ones(1,12),2*ones(1,12),3*ones(1,12),4*ones(1,11),5*ones(1,7),6*ones(1,7)];
frqbnd_names = {'Delta','Theta','LowAlpha','HighAlpha','LowBeta','HighBeta','LowGamma'};

%% Concatenate Subject and Condition Dimensions together
rel_pow_cat=cat(1,glob_rel_pow(:,:,1),glob_rel_pow(:,:,2),glob_rel_pow(:,:,3),...
    glob_rel_pow(:,:,4),glob_rel_pow(:,:,5),glob_rel_pow(:,:,6));
rel_pow_cat(sum(rel_pow_cat,2)==0,:)=[];

%% Linear regression and anova
for ftr = 1:length(frqbnd_names)
    eval([frqbnd_names{ftr},' = rel_pow_cat(:,ftr).*100;']);
end


eval(['d = dataset(',frqbnd_names{1},',',frqbnd_names{2},',',frqbnd_names{3},',',...
    frqbnd_names{4},',',frqbnd_names{5},',',frqbnd_names{6},',',frqbnd_names{7},');']);
d.Condition = (condNo');
d.Subject = (subjNo');
d.thetalobetaratio = log(theta./lobeta);
d.thetahibetaratio = log(theta./hibeta);
d.thetabetaratio = log(theta./(lobeta + hibeta));
d.thetaloalpharatio = log(theta./loalpha);
d.thetahialpharatio = log(theta./hialpha);
d.thetaalpharatio = log(theta./(loalpha + hialpha));


% eval(['lm = LinearModel.fit(d,''Condition ~ ',frqbnd_names{1},'+',frqbnd_names{2},'+',frqbnd_names{3},'+',...
%    frqbnd_names{4},'+',frqbnd_names{5},'+',frqbnd_names{6},'+',frqbnd_names{7},' + thetalobetaratio + ',...
%    'thetahibetaratio + thetabetaratio + thetaloalphaaratio + thetahialphaaratio + thetaalphaaratio'')']);

%% Plot Interactions of each Power Band to the Conditions
myfig = figure(345);
for fb = 1:length(frqbnd_names)
    subplot(2,4,fb)
    eval(['lm = LinearModel.fit(d,''',frqbnd_names{fb},' ~ Condition'');']);
    plot(lm)
    if fb ~= 1, legend('off'), end
    set(gca,'xlim',[0.5 6.5],'xtick',(1:6),'xticklabel',MusicCondLabels)
    ylabel(['Relative Power at ',frqbnd_labels{fb},' Hz (%)'])
    rotateXLabels(gca,45)
    
    % save table of statistics
    %                 t-test statistic, p-value, R^2, R^2 adjusted, B1, B0, RMSE
    statsumm(fb,:) = [lm.Coefficients{2,3:4},lm.Rsquared.Ordinary,lm.Rsquared.Adjusted,lm.Coefficients{2,1},lm.Coefficients{1,1},lm.RMSE];
end
h = subplot(2,4,8);
set(h,'color',[.94 .94 .94],'xcolor',[.94 .94 .94],'ycolor',[.94 .94 .94])
cnames = {'t','p-Value','R^2','R^2 Adj.','B1','B0','RMSE'};
t = uitable('Parent',myfig,'Data',statsumm,'ColumnName',cnames,... 
            'RowName',frqbnd_names,'Position',[1227 198 451 148],'ColumnWidth',{50});
set(gcf,'position',[1 1 1722 839])

%% save it
cd('C:\Users\TMHZXH5\Documents')
filename = 'LinearRegressionBetweenGlobalPowerBands2AllSubjects_ASRdenoised.tif';
set(gcf,'paperpositionmode','auto')
print(myfig,'-dpng','-r300',filename)