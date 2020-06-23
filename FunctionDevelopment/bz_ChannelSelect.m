function [selectedchannel] = bz_ChannelSelect(channelheatmap,chanIDs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%load('YMV09_171204.AnalysisResults.ISILFPMap.mat')
%%
%channelheatmap = MIMap.NA.NREMstate.pE;
%sorting = MIMap.NA.SGorder;
%chanIDs = MIMap.NA.ChanID;

%%
figure
imagesc(channelheatmap')
hold on
title('Select A Channel')
[~,chanclick] = ginput(1);
chanclick = round(chanclick);
plot(xlim(gca),chanclick+[0.5 0.5],'r')
plot(xlim(gca),chanclick-[0.5 0.5],'r')

%selectedchannel = MIMap.NA.SGorder(chanclick);
selectedchannel = chanIDs(chanclick);
end

