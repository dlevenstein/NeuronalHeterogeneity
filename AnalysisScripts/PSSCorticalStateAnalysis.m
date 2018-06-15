function [ output_args ] = PSSCorticalStateAnalysis( basePath,figfolder )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% DEV
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
basePath = '/mnt/proraidDL/Database/BWCRCNS/JennBuzsaki22/20140526_277um';
figfolder = '/Users/dlevenstein/Dropbox/Research/Current Projects/FRHetAndDynamics/AnalysisScripts/AnalysisFigs';
figfolder = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/AnalysisScripts/AnalysisFigs/PSSCorticalStateAnalysis';

%%
baseName = bz_BasenameFromBasepath(basePath);


spikes = bz_GetSpikes('basePath',basePath);
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
lfp = bz_GetLFP(SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID,...
    'basepath',basePath);
%%
downsamplefactor = 5;
lfp_down = bz_DownsampleLFP(lfp,downsamplefactor);
%%
dt = 0.5;

numwins = 21;
winsizes = logspace(0,1.5,numwins);
for ww = 1:numwins
    ww
    [specslope_temp] = bz_PowerSpectrumSlope(lfp_down,winsizes(ww),dt,'showfig',false);
    if ww == 1
        specslope_wins=specslope_temp;
        specslope_wins.winsizes = winsizes;
    else
        specslope_wins.data(:,ww)= interp1(specslope_temp.timestamps,specslope_temp.data,specslope_wins.timestamps);
        specslope_wins.rsq(:,ww) = interp1(specslope_temp.timestamps,specslope_temp.rsq,specslope_wins.timestamps);
        specslope_wins.intercept(:,ww) = interp1(specslope_temp.timestamps,specslope_temp.intercept,specslope_wins.timestamps);
    end
end

%% Compare different windows

states = fieldnames(SleepState.ints);
statecolors = {'k','b','r'};

maxlag = 200;
clear PSSautocorr
clear PSSwincorr
for ss = 1:length(states)
    specslope_wins.timeidx.(states{ss}) = ...
        InIntervals(specslope_wins.timestamps,SleepState.ints.(states{ss})) &...
        ~isnan(specslope_wins.data(:,end));
    
    PSSwincorr.(states{ss}) = corr(specslope_wins.data(specslope_wins.timeidx.(states{ss}),:),...
        'type','spearman','rows','pairwise');
    
    [PSSxcorr,lags] = xcov(specslope_wins.data(specslope_wins.timeidx.(states{ss}),:),maxlag./dt,'coeff');
    PSSxcorr_mat = reshape(PSSxcorr,length(lags),numwins,numwins);
    for ll = 1:numwins
        PSSautocorr.(states{ss})(:,ll) = PSSxcorr_mat(:,ll,ll);
    end
end


%% Pick the "best" window for further analysis
dt = 0.5;
winsize = 5;
[specslope,specgram] = bz_PowerSpectrumSlope(lfp,winsize,dt,'showfig',true);


numbins = 20;
PSShist.bins = linspace(-2,0,numbins);
for ss = 1:length(states)
    specslope.timeidx.(states{ss}) = InIntervals(specslope.timestamps,SleepState.ints.(states{ss}));
    
    PSShist.(states{ss}) = hist(specslope.data(specslope.timeidx.(states{ss})),PSShist.bins);
end


%% Figure: Time Window Comparison
skipnum = 5;
colororder = makeColorMap([0 0.5 0],[0.8 0.5 0],ceil(numwins/skipnum));
exwinsize = 800;
exwin = bz_RandomWindowInIntervals(specslope_wins.timestamps([1 end])',exwinsize);


figure
subplot(4,1,1)
    imagesc(specgram.timestamps,log2(specgram.freqs),specgram.amp)
    hold on
    StateScorePlot({SleepState.ints.NREMstate,SleepState.ints.REMstate,SleepState.ints.WAKEstate},...
        {'b','r','k'})
    axis xy
    xlim(exwin)
    ylabel({'Specgram','f (Hz)'})
    LogScale('y',2)

subplot(4,1,2)
    set(gca,'colororder',colororder)
    hold all
    plot(specslope_wins.timestamps,specslope_wins.data(:,1:skipnum:end),'linewidth',1)
    xlim(exwin);ylim([-2 0])
    legend(num2str(winsizes(1:skipnum:end)'))
    ylabel('PSS')
for ss = 1:length(states)
    subplot(4,4,ss+8)
        imagesc(log10(winsizes),log10(winsizes),PSSwincorr.(states{ss}))
        LogScale('xy',10)
        colorbar
        caxis([0.5 1])
        title(states{ss})
        
    subplot(4,4,ss+12)
    set(gca,'colororder',colororder)
    hold all
        plot((lags).*dt,PSSautocorr.(states{ss})(:,1:skipnum:end))
        
        plot(maxlag.*[-1 1],[0 0],'k')
        ylim([-0.2 1])
        xlim(150*[-1 1])
        title(states{ss})
        xlabel('t lag (s)');ylabel('corr')
end
NiceSave('CompareWindow',figfolder,baseName)
%% Figure: Zoom in to compare state

figure    
for ss = 1:length(states)
    
exwinsize = 30;
sexwin.(states{ss}) = bz_RandomWindowInIntervals(SleepState.ints.(states{ss}),exwinsize);
    
subplot(6,3,(ss-1).*6+[1 2])
    bz_MultiLFPPlot(lfp,'timewin',sexwin.(states{ss}))
    ylabel(states{ss})

subplot(6,3,(ss-1).*6+[4 5])
    set(gca,'colororder',colororder)
    hold all
    plot(specslope_wins.timestamps,specslope_wins.data(:,1:skipnum:end),'linewidth',1)
    xlim(sexwin.(states{ss}));ylim([-2 0])
    ylabel('PSS')
    

    subplot(3,3,ss*3)
    set(gca,'colororder',colororder)
    hold all
        plot((lags).*dt,PSSautocorr.(states{ss})(:,1:skipnum:end))
        plot(maxlag.*[-1 1],[0 0],'k')
        ylim([-0.2 1])
        xlim(20*[-1 1])
        title(states{ss})
        xlabel('t lag (s)')
legend(num2str(winsizes(1:skipnum:end)'),'location','northoutside')
end
NiceSave('CompareWindow_states',figfolder,baseName)

%%
figure
hist(specslope.rsq)



%%
rescolormap = makeColorMap([0 0 0.8],[1 1 1],[0.8 0 0]);
exwinsize = 2000;
exwin = bz_RandomWindowInIntervals(specslope.timestamps([1 end])',exwinsize);
figure
subplot(4,1,1)
    imagesc(specgram.timestamps,log2(specgram.freqs),specgram.amp)
    hold on
    StateScorePlot({SleepState.ints.NREMstate,SleepState.ints.REMstate,SleepState.ints.WAKEstate},...
        {'b','r','k'})
    axis xy
    xlim(exwin)
    ylabel({'Specgram','f (Hz)'})
    LogScale('y',2)
subplot(4,1,2)
    imagesc(specslope.timestamps,log2(specslope.freqs),specslope.resid')
    %colorbar
    axis xy
    xlim(exwin)
    LogScale('y',2)
    ylabel({'Resid.','f (Hz)'})
    colormap(gca,rescolormap)
    caxis(1.*[-1 1])
subplot(6,1,4)
    plot(specslope.timestamps,specslope.data,'k')
    xlim(exwin)
    ylabel('PSS');
subplot(6,1,5)
    plot(specslope.timestamps,specslope.rsq,'k')
    xlim(exwin)
    ylabel('Rsq')
% subplot(6,1,6)
%     plot(specslope.timestamps,specslope.intercept,'k')
%     xlim(exwin)
%     %ylabel('Intercept')
subplot(6,2,11)
    for ss = 1:length(states)
    plot(PSShist.bins,PSShist.(states{ss}),statecolors{ss},'linewidth',2)
    hold on
    end
    xlabel('PSS')
    legend(states{:},'location','eastoutside')
    axis tight
    box off
    xlim([-1.75 -0.25])
    
NiceSave('PSSbyState',figfolder,baseName)
%% Relate PSS and Spiking

ISIStats.allspikes.PSS = cellfun(@(X) interp1(specslope.timestamps,specslope.data,X,'nearest'),ISIStats.allspikes.times,'UniformOutput',false);
spkwinsize = 1;
overlap = winsize./dt;
spikemat = bz_SpktToSpkmat(spikes,'binsize',winsize,'overlap',overlap);
spikemat.PSS = interp1(specslope.timestamps,specslope.data,spikemat.timestamps);

cellclasses = {'pE','pI'};
classcolors = {'k','r'};
for tt = 1:length(cellclasses)
    spikemat.poprate.(cellclasses{tt}) = sum(spikemat.data(:,CellClass.(cellclasses{tt})),2)./winsize./sum(CellClass.(cellclasses{tt}));
end
%%
ratePSScorr.ALL = corr(spikemat.PSS,spikemat.data,'type','spearman','rows','complete');
CV2PSScorr.ALL = cellfun(@(X,Y) corr(X,Y,'type','spearman','rows','complete'),ISIStats.allspikes.PSS,ISIStats.allspikes.CV2);
%% Mean binned CV2...
clear CV2mat
CV2mat.winsize = 1;
CV2mat.timestamps = spikemat.timestamps;
CV2mat.binedges = bsxfun(@(X,Y) X+Y,spikemat.timestamps,[-0.5 0.5].*CV2mat.winsize);
for tt = 1:length(cellclasses)
    allspikes.CV2.(cellclasses{tt}) = cat(1,ISIStats.allspikes.CV2{CellClass.(cellclasses{tt})});
    allspikes.times.(cellclasses{tt}) = cat(1,ISIStats.allspikes.times{CellClass.(cellclasses{tt})});
    [CV2mat.timestamps,CV2mat.(cellclasses{tt})] = ...
        BinDataTimes(allspikes.CV2.(cellclasses{tt}),allspikes.times.(cellclasses{tt}),CV2mat.binedges);
    CV2mat.rate.(cellclasses{tt}) = interp1(spikemat.timestamps,spikemat.poprate.(cellclasses{tt}),CV2mat.timestamps);
end
CV2mat.PSS = interp1(specslope.timestamps,specslope.data,CV2mat.timestamps);

%%
figure
subplot(2,2,1)
plot(CV2mat.PSS,CV2mat.pI,'r.')
subplot(2,2,2)
plot(CV2mat.PSS,CV2mat.pE,'k.')

subplot(2,2,3)
plot(log10(CV2mat.rate.pI),CV2mat.pI,'r.')
subplot(2,2,4)
plot(log10(CV2mat.rate.pE),(CV2mat.pE),'k.')

%%
exwinsize = 100;
exwin = bz_RandomWindowInIntervals(specslope.timestamps([1 end])',exwinsize);
figure
for tt = 1:length(cellclasses)
    subplot(4,1,tt)
    plot(allspikes.times.(cellclasses{tt}),allspikes.CV2.(cellclasses{tt}),'.','color',classcolors{tt})
    hold on
    plot(CV2mat.timestamps,CV2mat.(cellclasses{tt}),'color',classcolors{tt})
    xlim(exwin)
end

%%
PSScorrhist.bins = linspace(-1,1,40);
for ss = 1:length(states)
%ss = 1;
    spikemat.timeidx.(states{ss}) = InIntervals(spikemat.timestamps,SleepState.ints.(states{ss}));
    CV2mat.timeidx.(states{ss}) = InIntervals(CV2mat.timestamps,SleepState.ints.(states{ss}));
    instatespikes = cellfun(@(X) InIntervals(X,SleepState.ints.(states{ss})),ISIStats.allspikes.times,'UniformOutput',false);

    ratePSScorr.(states{ss}) = corr(spikemat.PSS(spikemat.timeidx.(states{ss})),...
        spikemat.data(spikemat.timeidx.(states{ss}),:),'type','spearman','rows','complete');
    CV2PSScorr.(states{ss}) = cellfun(@(X,Y,Z) corr(X(Z),Y(Z),'type','spearman','rows','complete'),...
        ISIStats.allspikes.PSS,ISIStats.allspikes.CV2,instatespikes);
    
    for tt = 1:length(cellclasses)
        PSScorrhist.rate.(states{ss}).(cellclasses{tt}) = ...
            hist(ratePSScorr.(states{ss})(CellClass.(cellclasses{tt})),PSScorrhist.bins);
        PSScorrhist.CV2.(states{ss}).(cellclasses{tt}) = ...
            hist(CV2PSScorr.(states{ss})(CellClass.(cellclasses{tt})),PSScorrhist.bins);
    end
 
end    
%% Figure: PSS and Spiking
    figure
for ss = 1:length(states)
    subplot(4,4,ss)
    for tt = 1:length(cellclasses)
    	plot(ISIStats.summstats.(states{ss}).meanCV2(CellClass.(cellclasses{tt})),ratePSScorr.(states{ss})(CellClass.(cellclasses{tt})),'.','color',classcolors{tt})
        hold on
    end
    plot(get(gca,'xlim'),[0 0],'k')
    title(states{ss})
    ylim([-0.3 0.3])
    xlabel('<CV2>');ylabel('Rate-PSS Corr')
    

    subplot(4,4,ss+4)
    for tt = 1:length(cellclasses)
    	plot(ISIStats.summstats.(states{ss}).meanCV2(CellClass.(cellclasses{tt})),CV2PSScorr.(states{ss})(CellClass.(cellclasses{tt})),'.','color',classcolors{tt})
        hold on
    end
    ylim([-0.3 0.3])
    plot(get(gca,'xlim'),[0 0],'k')
    xlabel('<CV2>');ylabel('CV2-PSS Corr')
    

    subplot(4,4,ss+8)
    for tt = 1:length(cellclasses)
    	plot(log10(ISIStats.summstats.(states{ss}).meanrate(CellClass.(cellclasses{tt}))),ratePSScorr.(states{ss})(CellClass.(cellclasses{tt})),'.','color',classcolors{tt})
        hold on
    end
    ylim([-0.3 0.3])
    plot(get(gca,'xlim'),[0 0],'k')
    xlabel('Mean Rate');ylabel('Rate-PSS Corr')
    
    subplot(4,4,ss+12)
    for tt = 1:length(cellclasses)
    	plot(log10(ISIStats.summstats.(states{ss}).meanrate(CellClass.(cellclasses{tt}))),CV2PSScorr.(states{ss})(CellClass.(cellclasses{tt})),'.','color',classcolors{tt})
        hold on
    end
    ylim([-0.3 0.3])
    plot(get(gca,'xlim'),[0 0],'k')
    xlabel('Mean Rate');ylabel('CV2-PSS Corr')
end

NiceSave('PSSandCells',figfolder,baseName)
    %% Example Cell
    cc = randsample(spikes.numcells,1);
    figure
    subplot(2,2,1)
        plot(ISIStats.allspikes.PSS{cc},log10(ISIStats.allspikes.ISIs{cc}),'.')
        
    subplot(2,2,2)
        plot(ISIStats.allspikes.PSS{cc},(ISIStats.allspikes.CV2{cc}),'.')
        
	subplot(2,2,3)
        plot(spikemat.PSS,log10(spikemat.data(:,cc)),'.')
    
    

%%
figure
subplot(3,3,1)
 for ss = 1:length(states)  
        ScatterWithLinFit(spikemat.PSS(spikemat.timeidx.(states{ss})),...
            log2(spikemat.poprate.pE(spikemat.timeidx.(states{ss}))),statecolors{ss})
        hold on
 end
        xlabel('PSS');ylabel('pE Rate (spk/cell/s)')
       xlim([-1.8 -0.2])
        LogScale('y',2)
        
subplot(3,3,2)
 for ss = 1:length(states)  
        ScatterWithLinFit(spikemat.PSS(spikemat.timeidx.(states{ss})),...
            log2(spikemat.poprate.pI(spikemat.timeidx.(states{ss}))),statecolors{ss})
        hold on
 end
        xlabel('PSS');ylabel('pI Rate (spk/cell/s)')
        xlim([-1.8 -0.2])
        LogScale('y',2)
        
subplot(3,3,3)
 for ss = 1:length(states)  
        ScatterWithLinFit(spikemat.PSS(spikemat.timeidx.(states{ss})),...
            log10(spikemat.poprate.pE(spikemat.timeidx.(states{ss}))./spikemat.poprate.pI(spikemat.timeidx.(states{ss}))),statecolors{ss})
        hold on
 end
        xlabel('PSS');ylabel('E/I SpkRatio')
        xlim([-1.8 -0.2])
        LogScale('y',2)
        
for ss = 1:length(states)
    subplot(6,3,6+ss)
    for tt = 1:length(cellclasses)
        plot(PSScorrhist.bins,PSScorrhist.rate.(states{ss}).(cellclasses{tt}),...
            classcolors{tt},'linewidth',2)
        hold on
    end
    plot([0 0],get(gca,'ylim'),'k')
    xlabel('PSS-Rate Corr');
    ylabel('# cells')
    title(states{ss})
    xlim([-0.5 0.5])
end

for ss = 1:length(states)
    subplot(6,3,9+ss)
    for tt = 1:length(cellclasses)
        plot(PSScorrhist.bins,PSScorrhist.CV2.(states{ss}).(cellclasses{tt}),...
            classcolors{tt},'linewidth',2)
        hold on
    end
    plot([0 0],get(gca,'ylim'),'k')
    xlabel('PSS-CV2 Corr');
    ylabel('# cells')
    xlim([-0.5 0.5])
end

for tt = 1:length(cellclasses)
    subplot(3,3,6+tt)
    for ss = 1:length(states)
    	plot(CV2mat.PSS(CV2mat.timeidx.(states{ss})),CV2mat.(cellclasses{tt})(CV2mat.timeidx.(states{ss})),'.','color',statecolors{ss})
        hold on
    end
    xlabel('PSS');ylabel(['<CV2>, ',cellclasses{tt}])

end

NiceSave('PSSandSpiking',figfolder,baseName)
%% PSS and UP/DOWN

SlowWaves = bz_LoadEvents(basePath,'SlowWaves');


%%
updown = {'DOWN','UP'};
UDcolor = {'b','r'};
for ss = 1:2
    SlowWaves.dur.(updown{ss}) = diff(SlowWaves.ints.(updown{ss}),1,2);
    SlowWaves.midpoint.(updown{ss}) = mean(SlowWaves.ints.(updown{ss}),2);
    SlowWaves.PSS.(updown{ss}) = interp1(specslope.timestamps,specslope.data,SlowWaves.midpoint.(updown{ss}));
end

%%
figure
subplot(2,2,1)
for ss = 1:2
    plot(SlowWaves.PSS.(updown{ss}),log10(SlowWaves.dur.(updown{ss})),'.','color',UDcolor{ss})
    hold on
end
xlabel('PSS');ylabel('Dur (s)')
LogScale('y',10)
legend(updown{:},'location','northwest')
NiceSave('PSSandUPDOWN',figfolder,baseName)

