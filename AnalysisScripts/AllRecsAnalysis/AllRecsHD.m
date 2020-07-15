reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/HeadDirectionTuningAnalysis'];


[HDALL,baseNames] = GetMatResults(figfolder,'HeadDirectionTuningAnalysis');
HDALL = bz_CollapseStruct(HDALL);

%%
ISIbyHD_align = bz_CollapseStruct(HDALL.ISIbyHD_align,'match','justcat',true);
MutInfo = bz_CollapseStruct(HDALL.MutInfo,'match','justcat',true);

%%
MIthresh_rate = 0.05;
MIthresh_ISI = 0.05;
MeanISIPlaceField = nanmean(ISIbyHD_align.Dist.pYX(:,:,MutInfo.ISI>MIthresh_ISI & MutInfo.Rate>MIthresh_rate),3);
MeanRatePlaceField = nanmean(ISIbyHD_align.Dist.SpikeRate(:,:,MutInfo.ISI>MIthresh_ISI & MutInfo.Rate>MIthresh_rate),3);
%%
[~,sortMI_ISI] = sort(MutInfo.ISI);
[~,sortMutInfo.Rate] = sort(MutInfo.Rate);
%%

%%

figure
subplot(2,2,4)
    plot(log10(MutInfo.Rate),log10(MutInfo.ISI),'.')
    hold on
    UnityLine
    xlabel('MI (Rate)');ylabel('MI (ISI)')

subplot(2,2,3)
    imagesc(ISIbyHD_align.Dist.Xbins(1,:,1),ISIbyHD_align.Dist.Ybins(1,:,1),MeanISIPlaceField')
    hold on
    imagesc(ISIbyHD_align.Dist.Xbins(1,:,1)+2*pi,ISIbyHD_align.Dist.Ybins(1,:,1),MeanISIPlaceField')
    imagesc(ISIbyHD_align.Dist.Xbins(1,:,1)-2*pi,ISIbyHD_align.Dist.Ybins(1,:,1),MeanISIPlaceField')
    plot(ISIbyHD_align.Dist.Xbins(1,:,1),-log10(MeanRatePlaceField),'r')
    plot(ISIbyHD_align.Dist.Xbins(1,:,1)+2*pi,-log10(MeanRatePlaceField),'r')
    plot(ISIbyHD_align.Dist.Xbins(1,:,1)-2*pi,-log10(MeanRatePlaceField),'r')
    LogScale('y',10,'nohalf',true)
    ylabel('ISI (s)')
    bz_AddRightRateAxis
    xlabel('Position relative to HD Peak (m)')
    xlim([-1.5*pi 1.5.*pi])
    bz_piTickLabel('x')
    

subplot(2,2,1)
    imagesc(ISIbyHD_align.Dist.Xbins(1,:,1),[0 1],squeeze(log10(ISIbyHD_align.Dist.SpikeRate(:,:,sortMI_ISI)))')
    ylabel('Sort by MIISI')
subplot(2,2,2)
    imagesc(ISIbyHD_align.Dist.Xbins(1,:,1),[0 1],squeeze(log10(ISIbyHD_align.Dist.SpikeRate(:,:,sortMutInfo.Rate)))')
    ylabel('Sort by MIRate')
    
    NiceSave('HDCoding',figfolder,[])
