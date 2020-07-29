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
MIkinds = {'Skaggs','Rate','ISI'};

MIthresh.Rate = 0.02;
MIthresh.ISI = 0.02;
MIthresh.ISI = 0.05;
MIthresh.Skaggs = 0.5;

MIthresh.numspks = 500;

MutInfo.goodcells = MutInfo.numspks>MIthresh.numspks & MutInfo.cellclass.pE';

for kk = 1:3

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
legend(cellISIStats.statenames{2:3},'location','southoutside')

subplot(3,3,8)
hold on
for ss = 4:5
    plot(meanISIhist.logbins,meanISIhist.(cellISIStats.statenames{ss}).log)
end
legend(cellISIStats.statenames{4:5},'location','southoutside')

subplot(3,3,9)
hold on
for ss = 6:7
    plot(meanISIhist.logbins,meanISIhist.(cellISIStats.statenames{ss}).log)
end
legend(cellISIStats.statenames{6:7},'location','southoutside')

for ss = 2:3
subplot(3,3,(ss-2)*3+1)
hold on
    imagesc(meanISIhist.logbins,[0 1],squeeze(cellISIStats.allISIhist.(cellISIStats.statenames{ss}).log)')
end
%legend(cellISIStats.statenames{1:3})

for ss = 4:5
subplot(3,3,(ss-4)*3+2)
hold on
    imagesc(meanISIhist.logbins,[0 1],squeeze(cellISIStats.allISIhist.(cellISIStats.statenames{ss}).log)')
end

for ss = 6:7
subplot(3,3,(ss-6)*3+3)
hold on
    imagesc(meanISIhist.logbins,[0 1],squeeze(cellISIStats.allISIhist.(cellISIStats.statenames{ss}).log)')
    colormap(gca,NREMhistcolors)
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
legend(cellISIStats.statenames{2:3},'location','southoutside')

subplot(3,3,8)
hold on
for ss = 4:5
    plot(meanISIhist.logbins,meanISIhist.(cellISIStats.statenames{ss}).log)
end
legend(cellISIStats.statenames{4:5},'location','southoutside')

subplot(3,3,9)
hold on
for ss = 6:7
    plot(meanISIhist.logbins,meanISIhist.(cellISIStats.statenames{ss}).log)
end
legend(cellISIStats.statenames{6:7},'location','southoutside')

for ss = 2:3
subplot(3,3,(ss-2)*3+1)
    imagesc(meanISIhist.logbins,meanISIhist.logbins,meanISIhist.(cellISIStats.statenames{ss}).return)
    axis xy
end
%legend(cellISIStats.statenames{1:3})

for ss = 4:5
subplot(3,3,(ss-4)*3+2)
    imagesc(meanISIhist.logbins,meanISIhist.logbins,meanISIhist.(cellISIStats.statenames{ss}).return)
    axis xy
end

for ss = 6:7
subplot(3,3,(ss-6)*3+3)
    imagesc(meanISIhist.logbins,meanISIhist.logbins,meanISIhist.(cellISIStats.statenames{ss}).return)
    axis xy
    colormap(gca,NREMhistcolors)
end
NiceSave('InOutField_Return',figfolder,[])

%%
figure
subplot(2,2,1)
plot(log10(MutInfo.numspks(MutInfo.goodcells)),log10(MutInfo.ISI(MutInfo.goodcells)),'.')
xlabel('numspks');ylabel('MI ISI')

subplot(2,2,2)
plot(log10(MutInfo.numspks(MutInfo.goodcells)),log10(MutInfo.Skaggs(MutInfo.goodcells)),'.')
xlabel('# Spks on Track');ylabel('Skaggs Info')

subplot(2,2,3)
scatter(log10(MutInfo.Skaggs(MutInfo.goodcells)),log10(MutInfo.ISI(MutInfo.goodcells)),5,log10(MutInfo.numspks(MutInfo.goodcells)))
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
plot(log10(MIthresh.Skaggs).*[1 1],ylim(gca),'r--')
plot(xlim(gca),log10(MIthresh.Rate).*[1 1],'r--')
xlabel('I Skaggs');ylabel('MI Rate')
LogScale('xy',10)
colorbar

subplot(3,2,4)
scatter(log10(MutInfo.Skaggs(MutInfo.goodcells)),log10(MutInfo.ISI(MutInfo.goodcells)),3,MutInfo.GSweight(MutInfo.goodcells),'filled')
axis tight
hold on
plot(log10(MIthresh.Skaggs).*[1 1],ylim(gca),'r--')
plot(xlim(gca),log10(MIthresh.ISI).*[1 1],'r--')
xlabel('I Skaggs');ylabel('MI ISI')
LogScale('xy',10)
colorbar

subplot(3,2,6)
scatter(log10(MutInfo.Rate(MutInfo.goodcells)),log10(MutInfo.ISI(MutInfo.goodcells)),3,MutInfo.GSweight(MutInfo.goodcells),'filled')
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
    LogScale('y',10,'nohalf',true)
    ylabel('ISI (s)')
    bz_AddRightRateAxis
    xlabel('Position relative to PF Peak (m)')
    xlim([-1.25 1.25])

subplot(3,3,kk)
    imagesc(ISIbyPOS_norm.Dist.Xbins(1,:,1),[0 1],squeeze(log10(ISIbyPOS_norm.Dist.SpikeRate(:,:,sortAR.(MIkinds{kk}))))')
    title([(MIkinds{kk}),'-tuned cells (',num2str(length(sortAR.(MIkinds{kk}))),')'])
    ylabel('Sorted by AR')
    xlim([-1.25 1.25])

    
end
    
    NiceSave('PlaceCoding',figfolder,[])
    
    
%% Groups

hilowAR = 0.55;
groups = {'ISINotRate','TunedHiAR','NonLoGS','TunedLoAR'};
tunedcells.ISINotRate = tunedcells.ISI' & ~tunedcells.Skaggs';
tunedcells.NonLoGS = MutInfo.GSweight>0.95 & MutInfo.GSrate<-0.5;
tunedcells.TunedHiAR = tunedcells.ISI' & MutInfo.GSweight<hilowAR;
tunedcells.TunedLoAR = tunedcells.ISI' & MutInfo.GSweight>hilowAR;
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
        xlim([-1 1])
        %bz_piTickLabel('x')



    
end  
NiceSave('ARGroups',figfolder,[])
