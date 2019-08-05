basePath = '/mnt/proraidDL/Database/AGData/Buddy/Buddy_06272013';
baseName = bz_BasenameFromBasepath(basePath);
load(fullfile(basePath,[baseName,'.SleepScoreLFP.LFP.mat']))
load(fullfile(basePath,[baseName,'.SleepState.states.mat']))
load(fullfile(basePath,[baseName,'.EMGFromLFP.LFP.mat']))
figfolder = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/AnalysisScripts/AnalysisFigs/Illustrations';
%%
[SleepScoreMetrics,StatePlotMaterials] = ClusterStates_GetMetrics(...
                                           basePath,SleepScoreLFP,EMGFromLFP,true,...
                                           'window',8,'smoothfact',10,...
                                           'IRASA',true,'ThIRASA',true);
                                       
%% Figure
xwin = bz_RandomWindowInIntervals(SleepScoreMetrics(1,1).t_clus([1 end]),2000);
xwin = [15800 16800];
shortwins = 16195;
fftidx = find(SleepScoreMetrics.t_clus==shortwins);
shortwins = shortwins +[-4 4];

shortwins2 = 16317;
fftidx2 = find(SleepScoreMetrics.t_clus==shortwins2);
shortwins2 = shortwins2 +[-4 4];

figure

subplot(4,1,1)
    imagesc(SleepScoreMetrics.t_clus,...
        log2(StatePlotMaterials.swFFTfreqs),...
        log10(StatePlotMaterials.swFFTspec))
    axis xy
        hold all
        StateScorePlot(SleepState.ints,{'k','b','r'})
        plot(SleepScoreMetrics.t_clus,...
            bz_NormToRange(SleepScoreMetrics.broadbandSlowWave),'r')

    xlim(xwin)
    SpecColorRange(log10(StatePlotMaterials.swFFTspec),[1 1])
    LogScale('y',2)
    bz_ScaleBar('s')
    ylabel('PSS (f, Hz)')



subplot(4,1,2)
    imagesc(SleepScoreMetrics.t_clus,...
        log2(StatePlotMaterials.thFFTfreqs),...
        (StatePlotMaterials.thFFTspec))
    axis xy
        hold all
        plot(SleepScoreMetrics.t_clus,...
            bz_NormToRange(SleepScoreMetrics.thratio),'r')

    xlim(xwin)
    %SpecColorRange((StatePlotMaterials(2,15).thFFTspec),[1.5 1.5])
    %colorbar
    LogScale('y',2)
    bz_ScaleBar('s')
    ylabel('Theta (f, Hz)')
    
    
subplot(6,2,7)
plot(SleepScoreLFP.t,SleepScoreLFP.swLFP,'k')
box off
set(gca,'yticklabel',[])

xlim(shortwins)
bz_ScaleBar('s')

subplot(3,4,10)
plot(log2(StatePlotMaterials.swFFTfreqs),log10(StatePlotMaterials.swFFTspec(:,fftidx)),'k')
hold on
plot(log2(StatePlotMaterials.swFFTfreqs),...
    log10(StatePlotMaterials.swFFTfreqs).*StatePlotMaterials.IRASAslope(fftidx)+StatePlotMaterials.IRASAintercept(fftidx))
%plot(log10(StatePlotMaterials.swFFTfreqs),(StatePlotMaterials.IRASAsmooth(:,fftidx)))
axis tight
box off
LogScale('x',2)
ylim([3 5])
xlabel('f (Hz)');ylabel('log(power)')


subplot(6,2,8)
plot(SleepScoreLFP.t,SleepScoreLFP.thLFP,'k')
box off
set(gca,'yticklabel',[])
xlim(shortwins2)
bz_ScaleBar('s')

subplot(3,4,11)
plot(log2(StatePlotMaterials.swFFTfreqs),log10(StatePlotMaterials.swFFTspec(:,fftidx2)),'k')
hold on
plot(log2(StatePlotMaterials.swFFTfreqs),...
    log10(StatePlotMaterials.swFFTfreqs).*StatePlotMaterials.IRASAslope(fftidx2)+StatePlotMaterials.IRASAintercept(fftidx2))
%plot(log10(StatePlotMaterials.swFFTfreqs),(StatePlotMaterials.IRASAsmooth(:,fftidx2)))
axis tight
box off
ylim([3 5])
xlabel('f (Hz)');ylabel('log(power)')
LogScale('x',2)

 NiceSave('IllustrateIRASA',figfolder,baseName)
