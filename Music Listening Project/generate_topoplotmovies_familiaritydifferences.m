%% Load ERP Analysis data
clc, clear
% load('ERPdata_FamiliarSelfSelectvsUnfamiliarBach_correctedcluststats.mat')
% load('ERPdata_FamiliarSelfSelectvsUnfamiliarGagaku_correctedcluststats.mat')
load('ERPdata_FamiliarvsUnfamiliar_correctedcluststats.mat')

%% calculate the difference between conditions
cfg = [];
cfg.operation = 'x1-x2';
cfg.parameter = 'individual';
difference = ft_math(cfg, GACondGroup1, GACondGroup2);
cfg = [];
diff_timelock = ft_timelockanalysis(cfg, difference);

%% Add mask to difference dataset
Tdiffintrsect = ismember(difference.time,stat.time);
mask = zeros(length(stat.label),length(Tdiffintrsect)); mcntr=1;
for ts = 1:length(Tdiffintrsect),
    if Tdiffintrsect(ts)
mask(:,ts) = stat.mask(:,mcntr); mcntr=mcntr+1;  end
end
diff_timelock.mask = mask;

%% Display topoplot movie (by saving and compiling frames) 
mkdir('C:\Users\TMHZXH5\Documents\','TopoplotFrames')
cd('C:\Users\TMHZXH5\Documents\TopoplotFrames')
load easycap16HMRImusicproj_layout

time_segment = stat.time;
for ts =1:length(time_segment)
    myfig=figure;  
    cfg = [];
    cfg.xlim = [stat.time(ts) stat.time(ts)];
%     cfg.zlim = [min(diff_TL.avg(:)) max(diff_TL.avg(:))];  % across whole epoch
    cfg.zlim = [min(min(diff_timelock.avg(:,Tdiffintrsect))) max(max(diff_timelock.avg(:,Tdiffintrsect)))]; %for only specified segment 
    cfg.marker = 'labels';
    cfg.layout = lay;
    cfg.comment = 'xlim';
    cfg.commentpos = 'middlebottom';
    cfg.colorbar = 'yes';
    cfg.colormap = flipud(brewermap(64,'RdBu'));
    ft_topoplotER(cfg, diff_timelock)
    print(myfig,'-dpng', ['FamiliarDiffTopoPlot_frame',num2str(ts),'.png'], '-r300')
    close all
end

%% Create a video writer object
writerObj = VideoWriter('FamiliarSelfSelectvsGagakuVideo.avi');

writerObj.FrameRate = 15;   % Set frame rate

% Open video writer object and write frames sequentially
open(writerObj)
cd('C:\Users\TMHZXH5\Documents\TopoplotFrames')
for ts = 1:length(time_segment)
     % Read frame
     frame = sprintf('FamiliarDiffTopoPlot_frame%d.png', ts);
     input = imread(frame);

     % Write frame now
     writeVideo(writerObj, input);
end

% Close the video writer object
close(writerObj);
%{
%% Display topoplot movie (using FieldTrip command)
figure
load easycap16HMRImusicproj_layout
cfg = [];
cfg.xlim = [stat.time(1) stat.time(end)];
% cfg.xlim = [diff_timelock.time(1) diff_timelock.time(end)];
cfg.layout=lay;
cfg.parameter = 'avg';
cfg.colorbar     = 'yes';
ft_movieplotER(cfg, diff_timelock);
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
%}