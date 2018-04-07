%% Load extracted features list
clc, clear
cd('Z:\Research_Projects\MUSIC_Karmonik\')
load('extractfeats_ASRdenoisedEEG.mat')

%% Concatenate subject-condition dimensions together for all features
feat_names = fields(extractedfeats);
feat_size = length(feat_names);
allpowbands_array = [];
%%
for ftr = 9:15%1:feat_size
    %% Take only subjects with all conditions
    subjWOallconds = [1,5,9,10,12];
    eval(['extractedfeats.',feat_names{ftr},'(subjWOallconds,:,:)=[];']);
    %%
    allcond_array = [];
    for cndtn = 1:6
        eval(['condarray = squeeze(extractedfeats.',feat_names{ftr},'(:,cndtn,:));']);
        condvect = reshape(condarray,numel(condarray),1); %reshape s.t. each condition contains groups of subjects per channel
        allcond_array = [allcond_array, condvect];
    end
    allpowbands_array = [allpowbands_array; allcond_array];
end

dat = allpowbands_array;


%% Load grand-averaged power list
load('global_rel_pow_ASRdenoisedEEG.mat')

%% Initializations
frqbnd_labels = {'1-4';'4-7';'7-12'; '12-15'; '15-24'; '24-30'; '30-50'};
MusicCondLabels = {'Self-Selected','Bach','Gagaku','Clicking Langauge','Cronkite','Chaplin'};
subjNo = [(1:12),(1:12),(1:12),(1:4),(6:12),(2:4),(6:8),11,(2:4),(6:8),11];
condNo = [ones(1,12),2*ones(1,12),3*ones(1,12),4*ones(1,11),5*ones(1,7),6*ones(1,7)];
frqbnd_names = {'Delta','Theta','LowAlpha','HighAlpha','LowBeta','HighBeta','LowGamma'};

%% Take only subjects with all conditions
subjWOallconds = [1,5,9,10,12];
glob_rel_pow(subjWOallconds,:,:)=[];

%% Concatenate Subject and Freq Band Dimensions together
dat = cat(1,squeeze(glob_rel_pow(1,:,:)),squeeze(glob_rel_pow(2,:,:)),squeeze(glob_rel_pow(3,:,:)),...
squeeze(glob_rel_pow(4,:,:)),squeeze(glob_rel_pow(5,:,:)),squeeze(glob_rel_pow(6,:,:)),squeeze(glob_rel_pow(7,:,:)));

%% Remove freq. bands
idcs_rmv = sort([1:7:49,2:7:49,3:7:49,4:7:49,5:7:49,7:7:49]);
dat(idcs_rmv,:)=[];

%% Compute Two-Way ANOVA with Conditions as columns
[~,~,stats] = anova2(dat,1);

%% Compute Multiple Comparisons Test by Condition
c1 = multcompare(stats);

%% Compute Multiple Comparisons Test by Subject
c2 = multcompare(stats,'Estimate','row');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Concatenate Condition and Freq Band Dimensions together
dat2 = cat(1,glob_rel_pow(:,:,1),glob_rel_pow(:,:,2),glob_rel_pow(:,:,3),...
glob_rel_pow(:,:,4),glob_rel_pow(:,:,5),glob_rel_pow(:,:,6));
% rel_pow_cat(sum(rel_pow_cat,2)==0,:)=[];

%% Compute Two-Way ANOVA with Subjects as columns
[~,~,stats2] = anova2(dat2,7);

%% Compute Multiple Comparisons Test by Subject
c12 = multcompare(stats2);

%% Compute Multiple Comparisons Test by Condition
c22 = multcompare(stats2,'Estimate','row');