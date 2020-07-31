reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/HeadDirectionTuningAnalysis'];


[HDALL,baseNames] = GetMatResults(figfolder,'HeadDirectionTuningAnalysis');
HDALL = bz_CollapseStruct(HDALL);

%%
ISIbyHD_align = bz_CollapseStruct(HDALL.ISIbyHD_align,3,'justcat',true);
%ISIbyHD_alignGam = bz_CollapseStruct(HDALL.ISIbyHD_alignGam,3,'justcat',true);
MutInfo = bz_CollapseStruct(HDALL.MutInfo,'match','justcat',true);
cellISIStats = bz_CollapseStruct(HDALL.cellISIStats,3,'justcat',true);
%%
%MutInfo.Skaggs = MutInfo.SkaggsInf;
MIkinds = {'Skaggs','Rate','ISI'};

MIthresh.Rate = 0.02;
MIthresh.ISI = 0.02;
MIthresh.Skaggs = 1;

for kk = 1:3
    tunedcells.(MIkinds{kk}) = MutInfo.(MIkinds{kk})>MIthresh.(MIkinds{kk});
    MeanPlaceField.(MIkinds{kk}).pISI = nanmean(ISIbyHD_align.Dist.pYX(:,:,tunedcells.(MIkinds{kk})),3);
    MeanPlaceField.(MIkinds{kk}).Rate = nanmean(ISIbyHD_align.Dist.SpikeRate(:,:,tunedcells.(MIkinds{kk})),3);
%     MeanPlaceField.(MIkinds{kk}).pGS = nanmean(ISIbyHD_alignGam.GammaModes.GSweights(:,:,tunedcells.(MIkinds{kk})),3);
%     MeanPlaceField.(MIkinds{kk}).GSrate = nanmean(ISIbyHD_alignGam.GammaModes.GSlogrates(:,:,tunedcells.(MIkinds{kk})),3);

    [~,sortMutInfo.(MIkinds{kk})] = sort(MutInfo.(MIkinds{kk}));
    [~,sortAR.(MIkinds{kk})] = sort(MutInfo.GSweight);
    sortAR.(MIkinds{kk}) = sortAR.(MIkinds{kk})(ismember(sortAR.(MIkinds{kk}),find(tunedcells.(MIkinds{kk}))));
end

numcells = length(MutInfo.GSrate);

%%
cellISIStats.ARmod = squeeze((1-cellISIStats.GammaModes.GSweights(1,4,:))-(1-cellISIStats.GammaModes.GSweights(1,3,:)));


%%
%%
figure
subplot(3,3,1)
plot(cellISIStats.ARmod,log10(squeeze(MutInfo.Skaggs)),'.')
xlabel('AR Modulation');ylabel('MI Skaggs')
subplot(3,3,2)
plot(cellISIStats.ARmod,log10(squeeze(MutInfo.Rate)),'.')
xlabel('AR Modulation');ylabel('MI Rate')
subplot(3,3,3)
plot(cellISIStats.ARmod,log10(squeeze(MutInfo.ISI)),'.')
xlabel('AR Modulation');ylabel('MI ISI')
subplot(3,3,4)
plot(cellISIStats.ARmod,(squeeze(MutInfo.GSrate)),'.')
xlabel('AR Modulation');ylabel('GS Rate')

%%
clear diffAR
diffAR = squeeze((1-cellISIStats.GammaModes.GSweights(1,4,:))-(1-cellISIStats.GammaModes.GSweights(1,3,:)));

%for kk = 1:3
%meanISIhist = bz_CollapseStruct(cellISIStats.allISIhist(tunedcells.(MIkinds{kk})),3,'mean',true);

%end
%% In/Out Field
histcolors = flipud(gray);
NREMhistcolors = makeColorMap([1 1 1],[0 0 0.8]);

figure
colormap(gcf,histcolors)
for kk = 1:3
subplot(3,3,6+kk)
hold on
for ss = 3:4
    plot(meanISIhist.logbins,mean(cellISIStats.allISIhist.(cellISIStats.statenames{ss}).log(:,:,tunedcells.(MIkinds{kk})),3))
end
axis tight
legend(cellISIStats.statenames{3:4},'location','southoutside')



for ss = 3:4
subplot(3,3,(ss-2)*3+kk-3)
hold on
    imagesc(meanISIhist.logbins,[0 1],squeeze(cellISIStats.allISIhist.(cellISIStats.statenames{ss}).log(:,:,tunedcells.(MIkinds{kk})))')
    axis tight
end
%legend(cellISIStats.statenames{1:3})

end

%For each cell, find the intervals in the field, calculate in-field runing,
%out-field runining, non-running, ISI distributions
NiceSave('InOutField',figfolder,[])

%% In/Out Field: REturn Maps


figure
colormap(gcf,histcolors)
for kk = 1:3
subplot(3,3,6+kk)
hold on
for ss = 3:4
    plot(meanISIhist.logbins,mean(cellISIStats.allISIhist.(cellISIStats.statenames{ss}).log(:,:,tunedcells.(MIkinds{kk})),3))
end
axis tight
legend(cellISIStats.statenames{3:4},'location','southoutside')


for ss = 3:4
subplot(3,3,(ss-2)*3+kk-3)
    imagesc(meanISIhist.logbins,meanISIhist.logbins,mean(cellISIStats.allISIhist.(cellISIStats.statenames{ss}).return(:,:,tunedcells.(MIkinds{kk})),3))
    axis xy
    axis tight
    LogScale('xy',10,'exp',true)
end
%legend(cellISIStats.statenames{1:3})

end

NiceSave('InOutField_Return',figfolder,[])
%%
figure
for kk = 1:3
subplot(3,3,kk)
plot(squeeze(1-cellISIStats.GammaModes.GSweights(1,3,tunedcells.(MIkinds{kk}))),squeeze(1-cellISIStats.GammaModes.GSweights(1,4,tunedcells.(MIkinds{kk}))),'k.','markersize',8)
hold on
UnityLine
ylabel('In-Field AR');xlabel('Out-Field AR')
title(MIkinds{kk})

subplot(3,3,3+kk)
plot(squeeze(cellISIStats.GammaModes.GSlogrates(1,3,tunedcells.(MIkinds{kk}))),squeeze(cellISIStats.GammaModes.GSlogrates(1,4,tunedcells.(MIkinds{kk}))),'k.','markersize',8)
hold on
axis tight
UnityLine
box off
ylabel('In-Field GS Rate');xlabel('Out-Field GS Rate')
LogScale('xy',10,'exp',true,'nohalf',true)

diffASweight = (cellISIStats.GammaModes.ASweights(4,:,tunedcells.(MIkinds{kk}))-cellISIStats.GammaModes.ASweights(3,:,tunedcells.(MIkinds{kk})));
subplot(3,3,6+kk)
hold on
for aa = 1:5
scatter(-cellISIStats.GammaModes.ASlogrates(1,aa,tunedcells.(MIkinds{kk})),...
    log10(cellISIStats.GammaModes.ASCVs(1,aa,tunedcells.(MIkinds{kk}))),...
    30*cellISIStats.GammaModes.ASweights(3,aa,tunedcells.(MIkinds{kk}))+eps,...
    squeeze(diffASweight(1,aa,:)),'filled')
end
axis tight
plot(xlim(gca),[0 0],'k--')
ColorbarWithAxis([-0.15 0.15],'Delta AR','inclusive',{'<','>'})
crameri('berlin','pivot',0)
LogScale('x',10,'exp',true,'nohalf',true)
LogScale('y',10,'exp',false,'nohalf',true)
xlabel('AS Mean ISI (s)');ylabel('AS Mode CV')


end
NiceSave('GSASField',figfolder,[])
%% Figure Information Metrics
figure
for kk = 1:3
subplot(3,2,1+(kk-1)*2)
    imagesc(ISIbyHD_align.Dist.Xbins(1,:,1),[0 numcells],squeeze(log10(ISIbyHD_align.Dist.SpikeRate(:,:,sortMutInfo.(MIkinds{kk}))))')
    hold on
    plot(xlim(gca),numcells-sum(tunedcells.(MIkinds{kk})).*[1 1],'r--','linewidth',1)
    ylabel(['Sort by I ',(MIkinds{kk})])
    bz_piTickLabel('x')
end

subplot(3,2,2)
scatter(log10(MutInfo.Skaggs),log10(MutInfo.Rate),3,diffAR,'filled')
axis tight
hold on
plot(log10(MIthresh.Skaggs).*[1 1],ylim(gca),'r--')
plot(xlim(gca),log10(MIthresh.Rate).*[1 1],'r--')
xlabel('I Skaggs');ylabel('MI Rate')
LogScale('xy',10)
crameri('berlin','pivot',0)
ColorbarWithAxis([-0.3 0.3],'AR Modulation')

subplot(3,2,4)
scatter(log10(MutInfo.Skaggs),log10(MutInfo.ISI),3,diffAR,'filled')
axis tight
hold on
plot(log10(MIthresh.Skaggs).*[1 1],ylim(gca),'r--')
plot(xlim(gca),log10(MIthresh.ISI).*[1 1],'r--')
xlabel('I Skaggs');ylabel('MI ISI')
LogScale('xy',10)
crameri('berlin','pivot',0)
ColorbarWithAxis([-0.3 0.3],'AR Modulation')

subplot(3,2,6)
scatter(log10(MutInfo.Rate),log10(MutInfo.ISI),3,diffAR,'filled')
axis tight
hold on
plot(xlim(gca),log10(MIthresh.ISI).*[1 1],'r--')
plot(log10(MIthresh.Rate).*[1 1],ylim(gca),'r--')
ylabel('MI ISI');xlabel('MI Rate') 
LogScale('xy',10)
crameri('berlin','pivot',0)
ColorbarWithAxis([-0.3 0.3],'AR Modulation')

NiceSave('InfoNetrics',figfolder,[])
%% Figure: Information Metrics and GS/AS
figure
for kk = 1:3
subplot(3,3,1+(kk-1)*3)
%scatter(log10(MutInfo.(MIkinds{kk})),MutInfo.GSrate,5,1-MutInfo.GSweight,'filled')
plot(log10(MutInfo.(MIkinds{kk})),MutInfo.GSrate,'k.')
axis tight
hold on
plot(log10(MIthresh.(MIkinds{kk})).*[1 1],ylim(gca),'r--')
box off
xlabel(['I ',(MIkinds{kk})]);ylabel('GS Rate (HZ)')
LogScale('xy',10,'nohalf',true,'exp',true)



subplot(3,3,2+(kk-1)*3)
scatter(log10(MutInfo.(MIkinds{kk})),1-MutInfo.GSweight,5,MutInfo.GSrate,'filled')
%plot(log10(MutInfo.(MIkinds{kk})),1-MutInfo.GSweight,'.k')
hold on
plot(log10(MIthresh.(MIkinds{kk})).*[1 1],ylim(gca),'r--')
box off
axis tight
if kk==3
    plot(xlim(gca).*[0 1]+log10(MIthresh.(MIkinds{kk})).*[1 0],0.5.*[1 1],'k--')
end
caxis([-0.5 1.25])
%colorbar
LogScale('x',10,'nohalf',true,'exp',true)
LogScale('c',10)
xlabel(['I ',(MIkinds{kk})]);ylabel('Act. Ratio')


subplot(3,3,3+(kk-1)*3)
%scatter(log10(MutInfo.(MIkinds{kk})),MutInfo.GSrate,5,1-MutInfo.GSweight,'filled')
plot(log10(MutInfo.(MIkinds{kk})),cellISIStats.ARmod,'k.')
axis tight
hold on
plot(xlim(gca),[0 0],'k')
plot(log10(MIthresh.(MIkinds{kk})).*[1 1],ylim(gca),'r--')
box off
xlabel(['I ',(MIkinds{kk})]);ylabel('AR Modulation')
LogScale('x',10,'nohalf',true,'exp',true)

%LogScale('xy',10)


end
NiceSave('HDGSAS',figfolder,[])

%%

figure
for kk = 1:3
subplot(3,3,kk+3)
    imagesc(ISIbyHD_align.Dist.Xbins(1,:,1),ISIbyHD_align.Dist.Ybins(1,:,1),MeanPlaceField.(MIkinds{kk}).pISI')
    hold on
    imagesc(ISIbyHD_align.Dist.Xbins(1,:,1)+2*pi,ISIbyHD_align.Dist.Ybins(1,:,1),MeanPlaceField.(MIkinds{kk}).pISI')
    imagesc(ISIbyHD_align.Dist.Xbins(1,:,1)-2*pi,ISIbyHD_align.Dist.Ybins(1,:,1),MeanPlaceField.(MIkinds{kk}).pISI')
    plot(ISIbyHD_align.Dist.Xbins(1,:,1),-log10(MeanPlaceField.(MIkinds{kk}).Rate),'r')
    plot(ISIbyHD_align.Dist.Xbins(1,:,1)+2*pi,-log10(MeanPlaceField.(MIkinds{kk}).Rate),'r')
    plot(ISIbyHD_align.Dist.Xbins(1,:,1)-2*pi,-log10(MeanPlaceField.(MIkinds{kk}).Rate),'r')
    LogScale('y',10,'nohalf',true,'exp',true)
    ylabel('ISI (s)')
    bz_AddRightRateAxis
    xlabel('Position relative to HD Peak (m)')
    xlim([-1.5*pi 1.5.*pi])
    bz_piTickLabel('x')
    
    
% subplot(3,3,6+(kk))
%     plot(ISIbyHD_alignGam.Dist.Xbins(1,:,1),MeanPlaceField.(MIkinds{kk}).pGS,'k')
%     hold on
%     plot(ISIbyHD_alignGam.Dist.Xbins(1,:,1)+2*pi,MeanPlaceField.(MIkinds{kk}).pGS,'k')
%     plot(ISIbyHD_alignGam.Dist.Xbins(1,:,1)-2*pi,MeanPlaceField.(MIkinds{kk}).pGS,'k')
%     box off
%     xlim([-1.5*pi 1.5.*pi])
%     ylim([0.1 0.5])
%     bz_piTickLabel('x')
%     ylabel('pGS')


subplot(3,3,kk)
    imagesc(ISIbyHD_align.Dist.Xbins(1,:,1),[0 length(sortAR.(MIkinds{kk}))],...
        squeeze(log10(ISIbyHD_align.Dist.SpikeRate(:,:,sortAR.(MIkinds{kk}))))')
    hold on
    %ylabel(['Sort by I ',(MIkinds{kk})])
    bz_piTickLabel('x')
    title([(MIkinds{kk}),'-tuned cells (',num2str(length(sortAR.(MIkinds{kk}))),')'])
    ylabel('Sorted by AR')

    
end  
    
NiceSave('HDCoding',figfolder,[])

%% Polar plots
figure
for kk = 1:3

% subplot(3,3,3+(kk))
%     polarplot(ISIbyHD_alignGam.Dist.Xbins(1,:,1),1-MeanPlaceField.(MIkinds{kk}).pGS,...
%         'k','linewidth',2)
% ax = gca;
% ax.ThetaAxisUnits = 'radians';
% %ax.ThetaDir = 'clockwise';
% ax.ThetaZeroLocation = 'top';
% rlim([0 1])


% subplot(3,3,6+(kk))
%     polarplot(ISIbyHD_alignGam.Dist.Xbins(1,:,1),MeanPlaceField.(MIkinds{kk}).GSrate,...
%         'k','linewidth',2)
% ax = gca;
% ax.ThetaAxisUnits = 'radians';
% %ax.ThetaDir = 'clockwise';
% ax.ThetaZeroLocation = 'top';
% rlim([-0.3 1])

end

%% Groups

hilowAR = 0.5;
groups = {'NonHiGS','TunedHiAR','NonLoGS','TunedLoAR'};
tunedcells.NonHiGS = MutInfo.GSweight>0.95 & MutInfo.GSrate>-0.5;
tunedcells.NonLoGS = MutInfo.GSweight>0.95 & MutInfo.GSrate<-0.5;
tunedcells.TunedHiAR = tunedcells.ISI & MutInfo.GSweight<hilowAR;
tunedcells.TunedLoAR = tunedcells.ISI & MutInfo.GSweight>hilowAR;
%tunedcells.TunedLoAR = ~tunedcells.ISI & MutInfo.GSweight<0.2;


for kk = 1:4
    %tunedcells.(groups{kk}) = MutInfo.(groups{kk})>MIthresh.(groups{kk});
    MeanPlaceField.(groups{kk}).pISI = nanmean(ISIbyHD_align.Dist.pYX(:,:,tunedcells.(groups{kk})),3);
    MeanPlaceField.(groups{kk}).Rate = nanmean(ISIbyHD_align.Dist.SpikeRate(:,:,tunedcells.(groups{kk})),3);
%     MeanPlaceField.(groups{kk}).pGS = nanmean(ISIbyHD_alignGam.GammaModes.GSweights(:,:,tunedcells.(groups{kk})),3);

    %[~,sortMutInfo.(groups{kk})] = sort(MutInfo.(groups{kk}));
    %[~,sortAR.(groups{kk})] = sort(MutInfo.GSweight);
    %sortAR.(groups{kk}) = sortAR.(groups{kk})(ismember(sortAR.(groups{kk}),find(tunedcells.(groups{kk}))));
end
%Non-activating (pGS>0.95)
%Divide into high and low GS rate

%ISI-tuned
%Divide into high AR, low AR, non-activating?
%%
figure
ScatterWithLinFit(MutInfo.GSweight(tunedcells.ISI),MutInfo.peakwidth(tunedcells.ISI))
box off
xlabel('GS Weight');ylabel('Peak Width')
%%
figure
for kk = 1:4
    subplot(3,2,kk)
        imagesc(ISIbyHD_align.Dist.Xbins(1,:,1),ISIbyHD_align.Dist.Ybins(1,:,1),MeanPlaceField.(groups{kk}).pISI')
        hold on
        imagesc(ISIbyHD_align.Dist.Xbins(1,:,1)+2*pi,ISIbyHD_align.Dist.Ybins(1,:,1),MeanPlaceField.(groups{kk}).pISI')
        imagesc(ISIbyHD_align.Dist.Xbins(1,:,1)-2*pi,ISIbyHD_align.Dist.Ybins(1,:,1),MeanPlaceField.(groups{kk}).pISI')
        plot(ISIbyHD_align.Dist.Xbins(1,:,1),-log10(MeanPlaceField.(groups{kk}).Rate),'r')
        plot(ISIbyHD_align.Dist.Xbins(1,:,1)+2*pi,-log10(MeanPlaceField.(groups{kk}).Rate),'r')
        plot(ISIbyHD_align.Dist.Xbins(1,:,1)-2*pi,-log10(MeanPlaceField.(groups{kk}).Rate),'r')
        LogScale('y',10,'nohalf',true)
        ylabel('ISI (s)')
        bz_AddRightRateAxis
        xlabel('Position relative to HD Peak (m)')
        xlim([-1.5*pi 1.5.*pi])
        bz_piTickLabel('x')


%     subplot(3,2,mod(kk+1,2)+5)
%         plot(ISIbyHD_alignGam.Dist.Xbins(1,:,1),MeanPlaceField.(groups{kk}).pGS,'k')
%         hold on
%         plot(ISIbyHD_alignGam.Dist.Xbins(1,:,1)+2*pi,MeanPlaceField.(groups{kk}).pGS,'k')
%         plot(ISIbyHD_alignGam.Dist.Xbins(1,:,1)-2*pi,MeanPlaceField.(groups{kk}).pGS,'k')
%         box off
%         xlim([-1.5*pi 1.5.*pi])
%         bz_piTickLabel('x')
    
end  
NiceSave('ARGroups',figfolder,[])