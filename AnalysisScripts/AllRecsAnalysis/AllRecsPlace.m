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
MIthresh.Skaggs = 1;

for kk = 1:3
    tunedcells.(MIkinds{kk}) = MutInfo.(MIkinds{kk})>MIthresh.(MIkinds{kk});
    MeanPlaceField.(MIkinds{kk}).pISI = nanmean(ISIbyPOS_norm.Dist.pYX(:,:,tunedcells.(MIkinds{kk})),3);
    MeanPlaceField.(MIkinds{kk}).Rate = nanmean(ISIbyPOS_norm.Dist.SpikeRate(:,:,tunedcells.(MIkinds{kk})),3);

    [~,sortMutInfo.(MIkinds{kk})] = sort(MutInfo.(MIkinds{kk}));
    [~,sortAR.(MIkinds{kk})] = sort(MutInfo.GSweight);
    sortAR.(MIkinds{kk}) = sortAR.(MIkinds{kk})(ismember(sortAR.(MIkinds{kk}),find(tunedcells.(MIkinds{kk}))));
end

numcells = length(MutInfo.GSrate);

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
scatter(log10(MutInfo.Skaggs),log10(MutInfo.Rate),3,MutInfo.GSweight,'filled')
axis tight
hold on
plot(log10(MIthresh.Skaggs).*[1 1],ylim(gca),'r--')
plot(xlim(gca),log10(MIthresh.Rate).*[1 1],'r--')
xlabel('I Skaggs');ylabel('MI Rate')
LogScale('xy',10)
colorbar

subplot(3,2,4)
scatter(log10(MutInfo.Skaggs),log10(MutInfo.ISI),3,MutInfo.GSweight,'filled')
axis tight
hold on
plot(log10(MIthresh.Skaggs).*[1 1],ylim(gca),'r--')
plot(xlim(gca),log10(MIthresh.ISI).*[1 1],'r--')
xlabel('I Skaggs');ylabel('MI ISI')
LogScale('xy',10)
colorbar

subplot(3,2,6)
scatter(log10(MutInfo.Rate),log10(MutInfo.ISI),3,MutInfo.GSweight,'filled')
axis tight
hold on
plot(xlim(gca),log10(MIthresh.ISI).*[1 1],'r--')
plot(log10(MIthresh.Rate).*[1 1],ylim(gca),'r--')
ylabel('MI ISI');xlabel('MI Rate') 
LogScale('xy',10)
colorbar

NiceSave('InfoNetrics',figfolder,[])

%%

figure
subplot(2,2,4)
    plot(log10(MutInfo.Rate),log10(MutInfo.ISI),'.')
    hold on
    UnityLine
    xlim([-4 0]);ylim([-4 0])


subplot(2,2,3)
    imagesc(ISIbyPOS_norm.Dist.Xbins(1,:,1),ISIbyPOS_norm.Dist.Ybins(1,:,1),MeanISIPlaceField')
    hold on
    plot(ISIbyPOS_norm.Dist.Xbins(1,:,1),-log10(MeanRatePlaceField),'r')
    LogScale('y',10,'nohalf',true)
    ylabel('ISI (s)')
    bz_AddRightRateAxis
    xlabel('Position relative to PF Peak (m)')
    xlim([-0.9 0.9])

subplot(2,2,1)
    imagesc(ISIbyPOS_norm.Dist.Xbins(1,:,1),[0 1],squeeze(log10(ISIbyPOS_norm.Dist.SpikeRate(:,:,sortMI_ISI)))')
    ylabel('Sort by MIISI')
subplot(2,2,2)
    imagesc(ISIbyPOS_norm.Dist.Xbins(1,:,1),[0 1],squeeze(log10(ISIbyPOS_norm.Dist.SpikeRate(:,:,sortMutInfo.Rate)))')
    ylabel('Sort by MIRate')
    
    NiceSave('PlaceCoding',figfolder,[])
