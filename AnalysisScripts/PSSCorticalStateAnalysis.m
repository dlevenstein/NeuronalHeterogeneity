function [ output_args ] = PSSCorticalStateAnalysis( basePath,figfolder )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% DEV
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/mnt/proraidDL/Database/BWCRCNS/JennBuzsaki22/20140526_277um';
figfolder = '/Users/dlevenstein/Dropbox/Research/Current Projects/FRHetAndDynamics/AnalysisScripts/AnalysisFigs';
%figfolder = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/AnalysisScripts/AnalysisFigs';
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

numwins = 10;
winsizes = logspace(-0.301,1.5,numwins);
for ww = numwins:-1:1
    ww
    [specslope_temp] = bz_PowerSpectrumSlope(lfp_down,winsizes(ww),dt,'showfig',false);
    if ww == numwins
        specslope=specslope_temp;
        specslope.winsizes = winsizes;
    else
        specslope.data(:,ww)= interp1(specslope_temp.timestamps,specslope_temp.data,specslope.timestamps,'nearest');
        specslope.rsq(:,ww) = interp1(specslope_temp.timestamps,specslope_temp.rsq,specslope.timestamps,'nearest');
        specslope.intercept(:,ww) = interp1(specslope_temp.timestamps,specslope_temp.intercept,specslope.timestamps,'nearest');
    end
end

%%
figure
hist(specslope.rsq)

%% PSS and State
dt = 0.5;
winsize = 4;
[specslope,specgram] = bz_PowerSpectrumSlope(lfp,winsize,dt,'showfig',true);

states = fieldnames(SleepState.ints);
statecolors = {'k','b','r'};
numbins = 20;
PSShist.bins = linspace(-2,0,numbins);
for ss = 1:length(states)
    specslope.timeidx.(states{ss}) = InIntervals(specslope.timestamps,SleepState.ints.(states{ss}));
    
    PSShist.(states{ss}) = hist(specslope.data(specslope.timeidx.(states{ss})),PSShist.bins);
end


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
    legend(states{:})
%% Relate PSS and Spiking

ISIStats.allspikes.PSS = cellfun(@(X) interp1(specslope.timestamps,specslope.data,X,'nearest'),ISIStats.allspikes.times,'UniformOutput',false);

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
%%
for ss = 1:length(states)
%ss = 1;
    spikemat.timeidx.(states{ss}) = InIntervals(spikemat.timestamps,SleepState.ints.(states{ss}));
    instatespikes = cellfun(@(X) InIntervals(X,SleepState.ints.(states{ss})),ISIStats.allspikes.times,'UniformOutput',false);

    ratePSScorr.(states{ss}) = corr(spikemat.PSS(spikemat.timeidx.(states{ss})),...
        spikemat.data(spikemat.timeidx.(states{ss}),:),'type','spearman','rows','complete');
    CV2PSScorr.(states{ss}) = cellfun(@(X,Y,Z) corr(X(Z),Y(Z),'type','spearman','rows','complete'),...
        ISIStats.allspikes.PSS,ISIStats.allspikes.CV2,instatespikes);
 
    

    figure
    subplot(2,2,1)
    for tt = 1:length(cellclasses)
    	plot(ISIStats.summstats.(states{ss}).meanCV2(CellClass.(cellclasses{tt})),ratePSScorr.(states{ss})(CellClass.(cellclasses{tt})),'.','color',classcolors{tt})
        hold on
    end
    xlabel('<CV2>');ylabel('Rate-PSS Corr')
    
    subplot(2,2,2)
    for tt = 1:length(cellclasses)
    	plot(log10(ISIStats.summstats.(states{ss}).meanrate(CellClass.(cellclasses{tt}))),ratePSScorr.(states{ss})(CellClass.(cellclasses{tt})),'.','color',classcolors{tt})
        hold on
    end
    xlabel('Mean Rate');ylabel('Rate-PSS Corr')
    
    subplot(2,2,3)
    for tt = 1:length(cellclasses)
    	plot(ISIStats.summstats.(states{ss}).meanCV2(CellClass.(cellclasses{tt})),CV2PSScorr.(states{ss})(CellClass.(cellclasses{tt})),'.','color',classcolors{tt})
        hold on
    end
    xlabel('<CV2>');ylabel('CV2-PSS Corr')
    subplot(2,2,4)
    for tt = 1:length(cellclasses)
    	plot(log10(ISIStats.summstats.(states{ss}).meanrate(CellClass.(cellclasses{tt}))),CV2PSScorr.(states{ss})(CellClass.(cellclasses{tt})),'.','color',classcolors{tt})
        hold on
    end
    xlabel('Mean Rate');ylabel('CV2-PSS Corr')
end
    %%
    cc = 15;
    figure
    subplot(2,2,1)
        plot(ISIStats.allspikes.PSS{cc},log10(ISIStats.allspikes.ISIs{cc}),'.')
        
    subplot(2,2,2)
        plot(ISIStats.allspikes.PSS{cc},(ISIStats.allspikes.CV2{cc}),'.')
        
	subplot(2,2,3)
        plot(spikemat.PSS,log10(spikemat.data(:,cc)),'.')
    
    
    %%
    figure
    subplot(2,2,1)
        plot(spikemat.PSS(spikemat.timeidx.(states{ss})),...
            log10(spikemat.poprate.pE(spikemat.timeidx.(states{ss}))),'k.')
    
    subplot(2,2,2)
        plot(spikemat.PSS(spikemat.timeidx.(states{ss})),...
            log10(spikemat.poprate.pI(spikemat.timeidx.(states{ss}))),'r.')
        
    subplot(2,2,3)
        ScatterWithLinFit(spikemat.PSS(spikemat.timeidx.(states{ss})),...
            log10(spikemat.poprate.pE(spikemat.timeidx.(states{ss}))./spikemat.poprate.pI(spikemat.timeidx.(states{ss}))),'r')
%end


%% PSS and UP/DOWN

SlowWaves = bz_LoadEvents(basePath,'SlowWaves');


%%
updown = {'UP','DOWN'};
UDcolor = {'r','b'};
for ss = 1:2
    SlowWaves.dur.(updown{ss}) = diff(SlowWaves.ints.(updown{ss}),1,2);
    SlowWaves.midpoint.(updown{ss}) = mean(SlowWaves.ints.(updown{ss}),2);
    SlowWaves.PSS.(updown{ss}) = interp1(specslope.timestamps,specslope.data,SlowWaves.midpoint.(updown{ss}));
end

%%
figure
for ss = 1:2
    plot(SlowWaves.PSS.(updown{ss}),log10(SlowWaves.dur.(updown{ss})),'.','color',UDcolor{ss})
    hold on
end
xlabel('PSS');ylabel('Dur (s)')
LogScale('y',10)
legend(updown{:},'location','northwest')
