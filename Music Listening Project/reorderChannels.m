function [ ChnOrder_Grouped, idcsCOG, COGlabels, COGlab_idcs ] = reorderChannels(allchns)
%% Group Channels
if ~exist('allchns','var')
    load('BrainVision_1020_16ChanLocs.mat')
    allchns = {chanLocs.labels}'; clear chanLocs
    allchns(ismember(allchns,'O1')) = [];
end
PF = {'FP2','FP1'};
FL = {'F7','F3'};
FR = {'F4','F8'};
CR = {'T8','C4'};
CL = {'C3','T7'};
PL = {'P7','P3'};
O  = {'O2'};
PR = {'P4','P8'};

ChnOrder_Grouped = [PF,FL,FR,CR,CL,PL,O,PR];
COGlabels = {'PF','FL','F','FR','CR','C','CL','PL','O','PR'};
COGlab_idcs = [1.5,3.5,5.5,7.5,9.5,11.5,13,14.5];

idcsCOG = zeros(length(allchns),1);
for ch = 1:length(allchns)
    idcsCOG(ch) = find(ismember(allchns,ChnOrder_Grouped(ch)));
end

end  %EOF