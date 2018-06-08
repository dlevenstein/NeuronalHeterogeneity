function [ classIDs,classmeans ] = ManualMapClassification( maps )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% DEV
maps = ISIstats.ISIhist.NREMstate.return;

%%
nummaps = size(maps,3);
classIDs = zeros(nummaps,1);
%% Manual Classifier
histcolors = flipud(gray);
figure
colormap(histcolors)
for mm = 1:nummaps
%mm = 1;

%This map
subplot(2,2,4)
    imagesc(maps(:,:,mm))
    set(gca,'ytick',[]);set(gca,'xtick',[]);
    axis xy
    title([num2str(mm),' of ',num2str(nummaps)])

%Mean of other map classes
mapclasses = unique(classIDs);
numclasses = length(mapclasses);
for cc = 1:numclasses
    inclassmaps = classIDs==mapclasses(cc);
    classmean = nanmean(maps(:,:,inclassmaps),3);
    
    subplot(4,4,cc)
        imagesc(classmean)
        title(num2str(mapclasses(cc)))
        set(gca,'ytick',[]);set(gca,'xtick',[]);
        axis xy

end

%User: which class is this (press 0 to leave uncategorized, come back to 0
%later?)
classIDs(mm) = input('Which Class?');

end

end

