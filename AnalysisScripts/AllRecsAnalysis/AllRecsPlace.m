reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/PlaceTuningAnalysis'];


[PlaceALL,baseNames] = GetMatResults(figfolder,'PlaceTuningAnalysis');
PlaceALL = bz_CollapseStruct(PlaceALL);

%%
ISIbyPOS_norm = bz_CollapseStruct(PlaceALL.ISIbyPOS_norm,'match','justcat',true);
MutInfo = bz_CollapseStruct(PlaceALL.MutInfo,'match','justcat',true);

%%
MIthresh = 0.03;
MeanISIPlaceField = nanmean(ISIbyPOS_norm.Dist.pYX(:,:,MutInfo.ISI>MIthresh & MutInfo.Rate'>MIthresh),3);
MeanRatePlaceField = nanmean(ISIbyPOS_norm.Dist.SpikeRate(:,:,MutInfo.ISI>MIthresh & MutInfo.Rate'>MIthresh),3);
%%
[~,sortMI_ISI] = sort(MutInfo.ISI);
[~,sortMutInfo.Rate] = sort(MutInfo.Rate);

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
