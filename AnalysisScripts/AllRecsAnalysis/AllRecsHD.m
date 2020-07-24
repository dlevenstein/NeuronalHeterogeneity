reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/HeadDirectionTuningAnalysis'];


[HDALL,baseNames] = GetMatResults(figfolder,'HeadDirectionTuningAnalysis');
HDALL = bz_CollapseStruct(HDALL);

%%
ISIbyHD_align = bz_CollapseStruct(HDALL.ISIbyHD_align,3,'justcat',true);
ISIbyHD_alignGam = bz_CollapseStruct(HDALL.ISIbyHD_alignGam,3,'justcat',true);
MutInfo = bz_CollapseStruct(HDALL.MutInfo,'match','justcat',true);

%%
MutInfo.Skaggs = MutInfo.SkaggsInf;
MIkinds = {'Skaggs','Rate','ISI'};

MIthresh.Rate = 0.03;
MIthresh.ISI = 0.02;
MIthresh.Skaggs = 2;

for kk = 1:3
    tunedcells.(MIkinds{kk}) = MutInfo.(MIkinds{kk})>MIthresh.(MIkinds{kk});
    MeanPlaceField.(MIkinds{kk}).pISI = nanmean(ISIbyHD_align.Dist.pYX(:,:,tunedcells.(MIkinds{kk})),3);
    MeanPlaceField.(MIkinds{kk}).Rate = nanmean(ISIbyHD_align.Dist.SpikeRate(:,:,tunedcells.(MIkinds{kk})),3);
    MeanPlaceField.(MIkinds{kk}).pGS = nanmean(ISIbyHD_alignGam.GammaModes.GSweights(:,:,tunedcells.(MIkinds{kk})),3);

    [~,sortMutInfo.(MIkinds{kk})] = sort(MutInfo.(MIkinds{kk}));
end

%% Figure Information Metrics
figure
subplot(3,2,1)
    imagesc(ISIbyHD_align.Dist.Xbins(1,:,1),[0 1],squeeze(log10(ISIbyHD_align.Dist.SpikeRate(:,:,sortMutInfo.Skaggs)))')
    ylabel('Sort by I Skaggs')
subplot(3,2,3)
    imagesc(ISIbyHD_align.Dist.Xbins(1,:,1),[0 1],squeeze(log10(ISIbyHD_align.Dist.SpikeRate(:,:,sortMutInfo.Rate)))')
    ylabel('Sort by MI Rate')
subplot(3,2,5)
    imagesc(ISIbyHD_align.Dist.Xbins(1,:,1),[0 1],squeeze(log10(ISIbyHD_align.Dist.SpikeRate(:,:,sortMutInfo.ISI)))')
    ylabel('Sort by MI ISI')

subplot(3,2,2)
scatter(log10(MutInfo.SkaggsInf),log10(MutInfo.Rate),2,log10(MutInfo.ISI))
xlabel('I Skaggs');ylabel('MI Rate')
subplot(3,2,4)
scatter(log10(MutInfo.SkaggsInf),log10(MutInfo.ISI),2,log10(MutInfo.Rate))
xlabel('I Skaggs');ylabel('MI ISI')
subplot(3,2,6)
scatter(log10(MutInfo.ISI),log10(MutInfo.Rate),2,log10(MutInfo.SkaggsInf))
xlabel('MI ISI');ylabel('MI Rate') 
%% Figure: Information Metrics and GS/AS
figure
subplot(3,2,1)
scatter(log10(MutInfo.SkaggsInf),MutInfo.GSrate,3,MutInfo.GSweight,'filled')
hold on
plot(log10(MIthresh.Skaggs).*[1 1],ylim(gca),'k--')
box off
xlabel('I Skaggs');ylabel('GS Rate (HZ)')
colorbar

subplot(3,2,2)
plot(log10(MutInfo.SkaggsInf),MutInfo.GSweight,'.')
hold on
plot(log10(MIthresh.Skaggs).*[1 1],ylim(gca),'k--')
box off
xlabel('I Skaggs');ylabel('GS Weight')


subplot(3,2,3)
scatter(log10(MutInfo.Rate),MutInfo.GSrate,3,MutInfo.GSweight,'filled')
hold on
plot(log10(MIthresh.rate).*[1 1],ylim(gca),'k--')
box off
xlabel('MI Rate') ;ylabel('GS Rate (HZ)')
colorbar

subplot(3,2,4)
plot(log10(MutInfo.Rate),MutInfo.GSweight,'.')
hold on
plot(log10(MIthresh.rate).*[1 1],ylim(gca),'k--')
box off
xlabel('MI Rate') ;ylabel('GS Weight')

subplot(3,2,5)
scatter(log10(MutInfo.ISI),MutInfo.GSrate,3,MutInfo.GSweight,'filled')
hold on
plot(log10(MIthresh.ISI).*[1 1],ylim(gca),'k--')
box off
xlabel('MI ISI');ylabel('GS Rate (HZ)')
colorbar

subplot(3,2,6)
plot(log10(MutInfo.ISI),MutInfo.GSweight,'.')
hold on
plot(log10(MIthresh.ISI).*[1 1],ylim(gca),'k--')
box off
xlabel('MI ISI');ylabel('GS Weight')
NiceSave('HDGSAS',figfolder,[])

%%

figure

subplot(2,2,1)
    imagesc(ISIbyHD_align.Dist.Xbins(1,:,1),ISIbyHD_align.Dist.Ybins(1,:,1),MeanPlaceField.ISI.pISI')
    hold on
    imagesc(ISIbyHD_align.Dist.Xbins(1,:,1)+2*pi,ISIbyHD_align.Dist.Ybins(1,:,1),MeanPlaceField.ISI.pISI')
    imagesc(ISIbyHD_align.Dist.Xbins(1,:,1)-2*pi,ISIbyHD_align.Dist.Ybins(1,:,1),MeanPlaceField.ISI.pISI')
    plot(ISIbyHD_align.Dist.Xbins(1,:,1),-log10(MeanPlaceField.ISI.Rate),'r')
    plot(ISIbyHD_align.Dist.Xbins(1,:,1)+2*pi,-log10(MeanPlaceField.ISI.Rate),'r')
    plot(ISIbyHD_align.Dist.Xbins(1,:,1)-2*pi,-log10(MeanPlaceField.ISI.Rate),'r')
    LogScale('y',10,'nohalf',true)
    ylabel('ISI (s)')
    bz_AddRightRateAxis
    xlabel('Position relative to HD Peak (m)')
    xlim([-1.5*pi 1.5.*pi])
    bz_piTickLabel('x')
    
subplot(2,2,3)
    plot(ISIbyHD_alignGam.Dist.Xbins(1,:,1),MeanPlaceField.ISI.pGS,'k')
    hold on
    plot(ISIbyHD_alignGam.Dist.Xbins(1,:,1)+2*pi,MeanPlaceField.ISI.pGS,'k')
    plot(ISIbyHD_alignGam.Dist.Xbins(1,:,1)-2*pi,MeanPlaceField.ISI.pGS,'k')
    box off
    xlim([-1.5*pi 1.5.*pi])
    bz_piTickLabel('x')
    
    
    
    NiceSave('HDCoding',figfolder,[])
