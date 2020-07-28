reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/PlaceTuningAnalysis'];


[PlaceALL,baseNames] = GetMatResults(figfolder,'PlaceTuningAnalysis');
PlaceALL = bz_CollapseStruct(PlaceALL);

%%
ISIbyPOS_norm = bz_CollapseStruct(PlaceALL.ISIbyPOS_norm,'match','justcat',true);
MutInfo = bz_CollapseStruct(PlaceALL.MutInfo,'match','justcat',true);

%%
%MutInfo.Skaggs = MutInfo.SkaggsInf;
MIkinds = {'Skaggs','Rate','ISI'};

MIthresh.Rate = 0.02;
MIthresh.ISI = 0.02;
MIthresh.Skaggs = 0.5;

MIthresh.numspks = 500;

MutInfo.goodcells = MutInfo.numspks>MIthresh.numspks & MutInfo.cellclass.pE';

for kk = 1:3
    if kk ==2
        tunedcells.(MIkinds{kk}) = MutInfo.goodcells' & MutInfo.(MIkinds{kk})>MIthresh.(MIkinds{kk});
    else
        tunedcells.(MIkinds{kk}) = MutInfo.goodcells & MutInfo.(MIkinds{kk})>MIthresh.(MIkinds{kk});
    end
    MeanPlaceField.(MIkinds{kk}).pISI = nanmean(ISIbyPOS_norm.Dist.pYX(:,:,tunedcells.(MIkinds{kk})),3);
    MeanPlaceField.(MIkinds{kk}).Rate = nanmean(ISIbyPOS_norm.Dist.SpikeRate(:,:,tunedcells.(MIkinds{kk})),3);

    [~,sortMutInfo.(MIkinds{kk})] = sort(MutInfo.(MIkinds{kk}));
    [~,sortAR.(MIkinds{kk})] = sort(MutInfo.GSweight);
    sortAR.(MIkinds{kk}) = sortAR.(MIkinds{kk})(ismember(sortAR.(MIkinds{kk}),find(tunedcells.(MIkinds{kk}))));
end

numcells = length(MutInfo.GSrate);


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
    xlim([-0.9 0.9])

subplot(3,3,kk)
    imagesc(ISIbyPOS_norm.Dist.Xbins(1,:,1),[0 1],squeeze(log10(ISIbyPOS_norm.Dist.SpikeRate(:,:,sortAR.(MIkinds{kk}))))')
    title([(MIkinds{kk}),'-tuned cells (',num2str(length(sortAR.(MIkinds{kk}))),')'])
    ylabel('Sorted by AR')

    
end
    
    NiceSave('PlaceCoding',figfolder,[])
    
    
%% Groups

hilowAR = 0.55;
groups = {'NonHiGS','TunedHiAR','NonLoGS','TunedLoAR'};
tunedcells.NonHiGS = MutInfo.GSweight>0.95 & MutInfo.GSrate>-0.5;
tunedcells.NonLoGS = MutInfo.GSweight>0.95 & MutInfo.GSrate<-0.5;
tunedcells.TunedHiAR = tunedcells.ISI & MutInfo.GSweight<hilowAR;
tunedcells.TunedLoAR = tunedcells.ISI & MutInfo.GSweight>hilowAR;
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
figure
ScatterWithLinFit(MutInfo.GSweight(tunedcells.ISI),MutInfo.peakwidth(tunedcells.ISI))
box off
xlabel('GS Weight');ylabel('Peak Width')
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
        xlim([-1.5*pi 1.5.*pi])
        bz_piTickLabel('x')


    subplot(3,2,mod(kk+1,2)+5)
        plot(ISIbyPOS_normGam.Dist.Xbins(1,:,1),MeanPlaceField.(groups{kk}).pGS,'k')
        hold on
        plot(ISIbyPOS_normGam.Dist.Xbins(1,:,1)+2*pi,MeanPlaceField.(groups{kk}).pGS,'k')
        plot(ISIbyPOS_normGam.Dist.Xbins(1,:,1)-2*pi,MeanPlaceField.(groups{kk}).pGS,'k')
        box off
        xlim([-1.5*pi 1.5.*pi])
        bz_piTickLabel('x')
    
end  
NiceSave('ARGroups',figfolder,[])
