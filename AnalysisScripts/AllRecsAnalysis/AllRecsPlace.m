reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/PlaceTuningAnalysis'];


[PlaceALL,baseNames] = GetMatResults(figfolder,'PlaceTuningAnalysis');
PlaceALL = bz_CollapseStruct(PlaceALL);

%%
ISIbyPOS_norm = bz_CollapseStruct(PlaceALL.ISIbyPOS_norm,'match','justcat',true);
MutInfo = bz_CollapseStruct(PlaceALL.MutInfo,'match','justcat',true);
cellISIStats = bz_CollapseStruct(PlaceALL.cellISIStats,3,'justcat',true);
%%
%MutInfo.Skaggs = MutInfo.SkaggsInf;
MIkinds = {'Skaggs','Rate','ISI','hasfield'};

MIthresh.Rate = 0.02;
MIthresh.ISI = 0.02;
MIthresh.ISI = 0.05;
MIthresh.Skaggs = 0.25;
MIthresh.hasfield = 0.5;

MIthresh.numspks = 500;

MutInfo.goodcells = MutInfo.numspks>MIthresh.numspks & MutInfo.cellclass.pE';

for kk = 1:4

        tunedcells.(MIkinds{kk}) = MutInfo.goodcells & MutInfo.(MIkinds{kk})>MIthresh.(MIkinds{kk});

    MeanPlaceField.(MIkinds{kk}).pISI = nanmean(ISIbyPOS_norm.Dist.pYX(:,:,tunedcells.(MIkinds{kk})),3);
    MeanPlaceField.(MIkinds{kk}).Rate = nanmean(ISIbyPOS_norm.Dist.SpikeRate(:,:,tunedcells.(MIkinds{kk})),3);

    [~,sortMutInfo.(MIkinds{kk})] = sort(MutInfo.(MIkinds{kk}));
    sortMutInfo.(MIkinds{kk}) = sortMutInfo.(MIkinds{kk})(ismember(sortMutInfo.(MIkinds{kk}),find(MutInfo.goodcells)));

    [~,sortAR.(MIkinds{kk})] = sort(MutInfo.GSweight);
    sortAR.(MIkinds{kk}) = sortAR.(MIkinds{kk})(ismember(sortAR.(MIkinds{kk}),find(tunedcells.(MIkinds{kk}))));
end

numcells = length(MutInfo.GSrate);


%%
meanISIhist = bz_CollapseStruct(cellISIStats.allISIhist,3,'mean',true);



%% In/Out Field
histcolors = flipud(gray);
NREMhistcolors = makeColorMap([1 1 1],[0 0 0.8]);

figure
colormap(gcf,histcolors)
subplot(3,3,7)
hold on
for ss = 2:3
    plot(meanISIhist.logbins,meanISIhist.(cellISIStats.statenames{ss}).log)
end
axis tight
legend(cellISIStats.statenames{2:3},'location','southoutside')

subplot(3,3,8)
hold on
for ss = 4:5
    plot(meanISIhist.logbins,meanISIhist.(cellISIStats.statenames{ss}).log)
end
axis tight
legend(cellISIStats.statenames{4:5},'location','southoutside')

subplot(3,3,9)
hold on
for ss = 6:7
    plot(meanISIhist.logbins,meanISIhist.(cellISIStats.statenames{ss}).log)
end
axis tight
legend(cellISIStats.statenames{6:7},'location','southoutside')

for ss = 2:3
subplot(3,3,(ss-2)*3+1)
hold on
    imagesc(meanISIhist.logbins,[0 1],squeeze(cellISIStats.allISIhist.(cellISIStats.statenames{ss}).log)')
    axis tight
end
%legend(cellISIStats.statenames{1:3})

for ss = 4:5
subplot(3,3,(ss-4)*3+2)
hold on
    imagesc(meanISIhist.logbins,[0 1],squeeze(cellISIStats.allISIhist.(cellISIStats.statenames{ss}).log)')
    axis tight
end


for ss = 6:7
subplot(3,3,(ss-6)*3+3)
hold on
    imagesc(meanISIhist.logbins,[0 1],squeeze(cellISIStats.allISIhist.(cellISIStats.statenames{ss}).log)')
    colormap(gca,NREMhistcolors)
    axis tight
end

%For each cell, find the intervals in the field, calculate in-field runing,
%out-field runining, non-running, ISI distributions
NiceSave('InOutField',figfolder,[])

%% In/Out Field: REturn Maps


figure
colormap(gcf,histcolors)
subplot(3,3,7)
hold on
for ss = 2:3
    plot(meanISIhist.logbins,meanISIhist.(cellISIStats.statenames{ss}).log)
end
    axis tight
    box off
legend(cellISIStats.statenames{2:3},'location','southoutside')

subplot(3,3,8)
hold on
for ss = 4:5
    plot(meanISIhist.logbins,meanISIhist.(cellISIStats.statenames{ss}).log)

end
    axis tight
    box off
legend(cellISIStats.statenames{4:5},'location','southoutside')

subplot(3,3,9)
hold on
for ss = 6:7
    plot(meanISIhist.logbins,meanISIhist.(cellISIStats.statenames{ss}).log)

end
    axis tight
    box off
legend(cellISIStats.statenames{6:7},'location','southoutside')

for ss = 2:3
subplot(3,3,(ss-2)*3+1)
    imagesc(meanISIhist.logbins,meanISIhist.logbins,meanISIhist.(cellISIStats.statenames{ss}).return)
    axis xy
    axis tight
    LogScale('xy',10,'exp',true)
end
%legend(cellISIStats.statenames{1:3})

for ss = 4:5
subplot(3,3,(ss-4)*3+2)
    imagesc(meanISIhist.logbins,meanISIhist.logbins,meanISIhist.(cellISIStats.statenames{ss}).return)
    axis xy
    axis tight
    LogScale('xy',10,'exp',true)
end

for ss = 6:7
subplot(3,3,(ss-6)*3+3)
    imagesc(meanISIhist.logbins,meanISIhist.logbins,meanISIhist.(cellISIStats.statenames{ss}).return)
    axis xy
    colormap(gca,NREMhistcolors)
    axis tight
    LogScale('xy',10,'exp',true)
end

NiceSave('InOutField_Return',figfolder,[])


%%
figure
subplot(3,3,1)
plot(squeeze(1-cellISIStats.GammaModes.GSweights(1,2,:)),squeeze(1-cellISIStats.GammaModes.GSweights(1,3,:)),'k.','markersize',8)
hold on
UnityLine
ylabel('In-Field AR');xlabel('Out-Field AR')

subplot(3,3,2)
plot(squeeze(cellISIStats.GammaModes.GSlogrates(1,2,:)),squeeze(cellISIStats.GammaModes.GSlogrates(1,3,:)),'k.','markersize',8)
hold on
axis tight
UnityLine
box off
ylabel('In-Field GS Rate');xlabel('Out-Field GS Rate')
LogScale('xy',10,'exp',true,'nohalf',true)

diffASweight = (cellISIStats.GammaModes.ASweights(3,:,:)-cellISIStats.GammaModes.ASweights(2,:,:));
subplot(3,3,3)
hold on
for aa = 1:5
scatter(-cellISIStats.GammaModes.ASlogrates(1,aa,:),log10(cellISIStats.GammaModes.ASCVs(1,aa,:)),30*cellISIStats.GammaModes.ASweights(3,aa,:)+eps,...
    squeeze(diffASweight(1,aa,:)),'filled')
end
axis tight
plot(xlim(gca),[0 0],'k--')
ColorbarWithAxis([-0.2 0.2],'Delta AR','inclusive',{'<','>'})
crameri('berlin','pivot',0)
LogScale('x',10,'exp',true,'nohalf',true)
LogScale('y',10,'exp',false,'nohalf',true)
xlabel('AS Mean ISI (s)');ylabel('AS Mode CV')

diffAR = (1-cellISIStats.GammaModes.GSweights(1,3,:))-(1-cellISIStats.GammaModes.GSweights(1,2,:));

% subplot(2,2,4)
% plot(squeeze(cellISIStats.MIskaggs),squeeze(diffAR),'.')
NiceSave('GSASField',figfolder,[])
%%
figure
subplot(2,2,1)
plot(log10(MutInfo.numspks(MutInfo.goodcells)),log10(MutInfo.ISI(MutInfo.goodcells)),'.')
hold on
plot(log10(MutInfo.numspks(MutInfo.hasfield)),log10(MutInfo.ISI(MutInfo.hasfield)),'ko','markersize',4)
xlabel('numspks');ylabel('MI ISI')

subplot(2,2,2)
plot(log10(MutInfo.numspks(MutInfo.goodcells)),log10(MutInfo.Skaggs(MutInfo.goodcells)),'.')
hold on
plot(log10(MutInfo.numspks(MutInfo.hasfield)),log10(MutInfo.Skaggs(MutInfo.hasfield)),'ko','markersize',4)
xlabel('# Spks on Track');ylabel('Skaggs Info')

subplot(2,2,3)
scatter(log10(MutInfo.Skaggs(MutInfo.goodcells)),log10(MutInfo.ISI(MutInfo.goodcells)),5,log10(MutInfo.numspks(MutInfo.goodcells)))
hold on
plot(log10(MutInfo.Skaggs(MutInfo.hasfield)),log10(MutInfo.ISI(MutInfo.hasfield)),'ko','markersize',4)
%%
%MIthresh = 0.03;
%MeanISIPlaceField = nanmean(ISIbyPOS_norm.Dist.pYX(:,:,MutInfo.ISI>MIthresh & MutInfo.Rate'>MIthresh),3);
%MeanRatePlaceField = nanmean(ISIbyPOS_norm.Dist.SpikeRate(:,:,MutInfo.ISI>MIthresh & MutInfo.Rate'>MIthresh),3);
%%
%[~,sortMI_ISI] = sort(MutInfo.ISI);
%[~,sortMutInfo.Rate] = sort(MutInfo.Rate);

%% Figure Information Metrics
figure
for kk = 1:3
subplot(3,2,1+(kk-1)*2)
    imagesc(ISIbyPOS_norm.Dist.Xbins(1,:,1),[0 numcells],squeeze(log10(ISIbyPOS_norm.Dist.SpikeRate(:,:,sortMutInfo.(MIkinds{kk}))))')
    hold on
    plot(xlim(gca),numcells-sum(tunedcells.(MIkinds{kk})).*[1 1],'r--','linewidth',1)
    ylabel(['Sort by I ',(MIkinds{kk})])
    bz_piTickLabel('x')
end

subplot(3,2,2)
scatter(log10(MutInfo.Skaggs(MutInfo.goodcells)),log10(MutInfo.Rate(MutInfo.goodcells)),3,MutInfo.GSweight(MutInfo.goodcells),'filled')
axis tight
hold on
plot(log10(MutInfo.Skaggs(MutInfo.hasfield)),log10(MutInfo.Rate(MutInfo.hasfield)),'ko','markersize',4)
hold on
plot(log10(MIthresh.Skaggs).*[1 1],ylim(gca),'r--')
plot(xlim(gca),log10(MIthresh.Rate).*[1 1],'r--')
xlabel('I Skaggs');ylabel('MI Rate')
LogScale('xy',10)
colorbar

subplot(3,2,4)
scatter(log10(MutInfo.Skaggs(MutInfo.goodcells)),log10(MutInfo.ISI(MutInfo.goodcells)),3,MutInfo.GSweight(MutInfo.goodcells),'filled')
hold on
plot(log10(MutInfo.Skaggs(MutInfo.hasfield)),log10(MutInfo.ISI(MutInfo.hasfield)),'ko','markersize',4)
axis tight
hold on
plot(log10(MIthresh.Skaggs).*[1 1],ylim(gca),'r--')
plot(xlim(gca),log10(MIthresh.ISI).*[1 1],'r--')
xlabel('I Skaggs');ylabel('MI ISI')
LogScale('xy',10)
colorbar

subplot(3,2,6)
scatter(log10(MutInfo.Rate(MutInfo.goodcells)),log10(MutInfo.ISI(MutInfo.goodcells)),3,MutInfo.GSweight(MutInfo.goodcells),'filled')
hold on
plot(log10(MutInfo.Rate(MutInfo.hasfield)),log10(MutInfo.ISI(MutInfo.hasfield)),'ko','markersize',4)
axis tight
hold on
plot(xlim(gca),log10(MIthresh.ISI).*[1 1],'r--')
plot(log10(MIthresh.Rate).*[1 1],ylim(gca),'r--')
ylabel('MI ISI');xlabel('MI Rate') 
LogScale('xy',10)
colorbar

NiceSave('InfoNetrics',figfolder,[])



%% Figure: Information Metrics and GS/AS
figure
for kk = 1:3
subplot(3,2,1+(kk-1)*2)
scatter(log10(MutInfo.(MIkinds{kk})(MutInfo.goodcells)),MutInfo.GSrate(MutInfo.goodcells),5,MutInfo.GSweight(MutInfo.goodcells),'filled')
axis tight
hold on
plot(log10(MIthresh.(MIkinds{kk})).*[1 1],ylim(gca),'r--')
box off
xlabel(['I ',(MIkinds{kk})]);ylabel('GS Rate (HZ)')
LogScale('xy',10)
colorbar

subplot(3,2,2+(kk-1)*2)
scatter(log10(MutInfo.(MIkinds{kk})(MutInfo.goodcells)),MutInfo.GSweight(MutInfo.goodcells),5,MutInfo.GSrate(MutInfo.goodcells),'filled')
hold on
plot(log10(MIthresh.(MIkinds{kk})).*[1 1],ylim(gca),'r--')
box off
axis tight
if kk==3
    plot(xlim(gca).*[0 1]+log10(MIthresh.(MIkinds{kk})).*[1 0],0.5.*[1 1],'k--')
end
caxis([-1 0.5])
colorbar
LogScale('x',10)
LogScale('c',10)
xlabel(['I ',(MIkinds{kk})]);ylabel('GS Weight')
end

NiceSave('PlaceGSAS',figfolder,[])


%%

figure
for kk = 1:3
subplot(3,3,kk+3)
    imagesc(ISIbyPOS_norm.Dist.Xbins(1,:,1),ISIbyPOS_norm.Dist.Ybins(1,:,1),MeanPlaceField.(MIkinds{kk}).pISI')
    hold on
    plot(ISIbyPOS_norm.Dist.Xbins(1,:,1),-log10(MeanPlaceField.(MIkinds{kk}).Rate),'r')
    LogScale('y',10,'nohalf',true,'exp',true)
    ylabel('ISI (s)')
    bz_AddRightRateAxis
    xlabel('Position relative to PF Peak (m)')
    xlim([-1.3 1.3])

subplot(3,3,kk)
    imagesc(ISIbyPOS_norm.Dist.Xbins(1,:,1),[0 1],squeeze(log10(ISIbyPOS_norm.Dist.SpikeRate(:,:,sortAR.(MIkinds{kk}))))')
    title([(MIkinds{kk}),'-tuned cells (',num2str(length(sortAR.(MIkinds{kk}))),')'])
    ylabel('Sorted by AR')
    xlim([-1.3 1.3])

    
end
    
    NiceSave('PlaceCoding',figfolder,[])
    
    
%% Groups

hilowAR = 0.5;
groups = {'ISINotRate','RateNotISI','TunedHiAR','TunedLoAR','ISInoPF'};
tunedcells.ISINotRate = tunedcells.ISI' & ~tunedcells.Skaggs';
tunedcells.RateNotISI = ~tunedcells.ISI' & tunedcells.Skaggs';
tunedcells.TunedHiAR = tunedcells.Rate' & MutInfo.GSweight<hilowAR;
tunedcells.TunedLoAR = tunedcells.Rate' & MutInfo.GSweight>hilowAR;
tunedcells.ISINoPF = tunedcells.ISI' & ~tunedcells.Skaggs';
%tunedcells.TunedLoAR = ~tunedcells.ISI & MutInfo.GSweight<0.2;


for kk = 1:4
    %tunedcells.(groups{kk}) = MutInfo.(groups{kk})>MIthresh.(groups{kk});
    MeanPlaceField.(groups{kk}).pISI = nanmean(ISIbyPOS_norm.Dist.pYX(:,:,tunedcells.(groups{kk})),3);
    MeanPlaceField.(groups{kk}).Rate = nanmean(ISIbyPOS_norm.Dist.SpikeRate(:,:,tunedcells.(groups{kk})),3);

    %[~,sortMutInfo.(groups{kk})] = sort(MutInfo.(groups{kk}));
    %[~,sortAR.(groups{kk})] = sort(MutInfo.GSweight);
    %sortAR.(groups{kk}) = sortAR.(groups{kk})(ismember(sortAR.(groups{kk}),find(tunedcells.(groups{kk}))));
end
%Non-activating (pGS>0.95)
%Divide into high and low GS rate

%ISI-tuned
%Divide into high AR, low AR, non-activating?
%%
% figure
% ScatterWithLinFit(MutInfo.GSweight(tunedcells.ISI),MutInfo.peakwidth(tunedcells.ISI))
% box off
% xlabel('GS Weight');ylabel('Peak Width')
%%
figure
for kk = 1:4
    subplot(3,2,kk)
        imagesc(ISIbyPOS_norm.Dist.Xbins(1,:,1),ISIbyPOS_norm.Dist.Ybins(1,:,1),MeanPlaceField.(groups{kk}).pISI')
        hold on
        imagesc(ISIbyPOS_norm.Dist.Xbins(1,:,1)+2*pi,ISIbyPOS_norm.Dist.Ybins(1,:,1),MeanPlaceField.(groups{kk}).pISI')
        imagesc(ISIbyPOS_norm.Dist.Xbins(1,:,1)-2*pi,ISIbyPOS_norm.Dist.Ybins(1,:,1),MeanPlaceField.(groups{kk}).pISI')
        plot(ISIbyPOS_norm.Dist.Xbins(1,:,1),-log10(MeanPlaceField.(groups{kk}).Rate),'r')
        plot(ISIbyPOS_norm.Dist.Xbins(1,:,1)+2*pi,-log10(MeanPlaceField.(groups{kk}).Rate),'r')
        plot(ISIbyPOS_norm.Dist.Xbins(1,:,1)-2*pi,-log10(MeanPlaceField.(groups{kk}).Rate),'r')
        LogScale('y',10,'nohalf',true)
        ylabel('ISI (s)')
        bz_AddRightRateAxis
        xlabel('Position relative to HD Peak (m)')
        xlim([-1.25 1.25])
        %bz_piTickLabel('x')
        title(groups{kk})



    
end  
NiceSave('ARGroups',figfolder,[])
