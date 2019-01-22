function [ ] = SpikeStatsbyLFPAnalysis(basePath,figfolder)

%% DEV
%reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity';
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/mnt/proraidDL/Database/BWCRCNS/JennBuzsaki22/20140526_277um';
%basePath = '/mnt/proraidDL/Database/BWCRCNS/Dino_mPFC/Dino_061814';
%basePath = pwd;
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/SpikeStatsbyLFPAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
states = fieldnames(SleepState.ints);
states{4} = 'ALL';
SleepState.ints.ALL = [0 Inf];
statecolors = {'k','b','r',[0.6 0.6 0.6]};

[celltypes,~,typeidx] = unique(CellClass.label);
cellcolor = {'k','r'};

% for reformatting SleepState
%SleepState = SleepScoreMaster(basePath,'noPrompts',true);
%LFP
%%
%Pick channel with most cells close to it...... for now;
%usechannel = mode(spikes.maxWaveformCh);
downsamplefactor = 2;
lfp = bz_GetLFP(SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID,...
    'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);
lfp.data = NormToInt(single(lfp.data),'modZ', SleepState.ints.NREMstate,lfp.samplingRate);

ISIStats = bz_LoadCellinfo(basePath,'ISIStats');


%%
ISIStats.allspikes.LFPval = cellfun(@(X) interp1(lfp.timestamps,single(lfp.data),X),...
    ISIStats.allspikes.times,'UniformOutput',false);

%%
binsize = 0.005;
spikemat = bz_SpktToSpkmat(spikes,'binsize',binsize);
spikemat.lfp = interp1(lfp.timestamps,single(lfp.data),spikemat.timestamps);

for ss = 1:length(states)
    spikemat.timeidx.(states{ss}) = InIntervals(spikemat.timestamps,SleepState.ints.(states{ss}));
    rateLFPcorr.(states{ss}) = corr(spikemat.lfp(spikemat.timeidx.(states{ss})),...
        spikemat.data(spikemat.timeidx.(states{ss}),:),'type','spearman','rows','complete');
end

%%
mintime = 5;
ratebyLFP.nbins = 30;
ratebyLFP.lfpbins = linspace(-3.5,3.5,ratebyLFP.nbins);

for ss = 1:4
    instatelfp = InIntervals(lfp.timestamps,SleepState.ints.(states{ss}));
    instatespikes = cellfun(@(X) InIntervals(X,SleepState.ints.(states{ss})),...
        ISIStats.allspikes.times,'UniformOutput',false);
    
    ratebyLFP.(states{ss}).lfphist = hist(single(lfp.data(instatelfp)),ratebyLFP.lfpbins);
    ratebyLFP.(states{ss}).lfphist = ratebyLFP.(states{ss}).lfphist./lfp.samplingRate;
    ratebyLFP.spikecounts = cellfun(@(X,Y) hist(X(Y),ratebyLFP.lfpbins),...
        ISIStats.allspikes.LFPval,instatespikes,'UniformOutput',false);
    ratebyLFP.(states{ss}).rate = cellfun(@(X) X./ratebyLFP.(states{ss}).lfphist,...
        ratebyLFP.spikecounts,'UniformOutput',false);
    ratebyLFP.(states{ss}).rate = cat(1,ratebyLFP.(states{ss}).rate{:});
    ratebyLFP.(states{ss}).rate(:,ratebyLFP.(states{ss}).lfphist<mintime) = nan;
    
    ratebyLFP.(states{ss}).lfphist = ratebyLFP.(states{ss}).lfphist./lfp.timestamps(end);
    
    for cc = 1:length(celltypes)
        ratebyLFP.(states{ss}).meanrate.(celltypes{cc}) = ...
            nanmean(ratebyLFP.(states{ss}).rate(CellClass.(celltypes{cc}),:),1);
        ratebyLFP.(states{ss}).std.(celltypes{cc}) = ...
            nanstd(ratebyLFP.(states{ss}).rate(CellClass.(celltypes{cc}),:),[],1);
    end
end
clear instatelfp
%%
%Conditional probability of instantaneous rate given LFP value
% [ CONDXY ] = ConditionalHist(ISIStats.allspikes.LFPval{cellex},...
%     log10(1./ISIStats.allspikes.meanISI{cellex}),'numXbins',ratebyLFP.nbins);

%%
figure

subplot(6,4,12)
    hold on
    for ss = 1:3
    plot(ratebyLFP.lfpbins,ratebyLFP.(states{ss}).lfphist,...
        'linewidth',2,'color',statecolors{ss})
    end
    xlabel('LFP (Z_N_R_E_M)');       axis tight
    xlim(ratebyLFP.lfpbins([1 end]))
    
for cc = 1:length(celltypes)
subplot(6,4,cc*4)
    hold on
%     for ss = 1:3
%     errorshade(ratebyLFP.lfpbins,(ratebyLFP.(states{ss}).meanrate),...
%         (ratebyLFP.(states{ss}).std),(ratebyLFP.(states{ss}).std),...
%         statecolors{ss},'scalar')
%     end
    for ss = 4:-1:1
    plot(ratebyLFP.lfpbins,(ratebyLFP.(states{ss}).meanrate.(celltypes{cc})),...
        'linewidth',2,'color',statecolors{ss})
    end
    axis tight
    xlim(ratebyLFP.lfpbins([1 end]));
    set(gca,'xtick',[]);  ylabel({celltypes{cc},'Rate (Hz)'})    
   % LogScale('y',10)
end

subplot(3,3,1)
    hold on
    for ss = 1:3
        [histcounts,histbins] = hist(rateLFPcorr.(states{ss}));
        plot(histbins,histcounts,'color',statecolors{ss},'linewidth',2)
    end
    xlabel('LFP-Rate Corr');ylabel('# Cells')
    
    
subplot(6,4,3)
    for tt = 1:length(celltypes)
%         plot(PSScorrhist.bins,PSScorrhist.CV2.(states{ss}).(cellclasses{tt}),...
%             classcolors{tt},'linewidth',2)
            HistWithMean(rateLFPcorr.ALL(CellClass.(celltypes{tt})),...
                'numbins',8,'color',cellcolor{tt},'showtext',false)
        hold on
    end
    plot([0 0],get(gca,'ylim'),'k')
    xlabel('LFP-Rate Corr');
    ylabel('# cells')
    xlim([-0.1 0.1])
    %legend('pE cells','pI ce
    
subplot(3,3,4)
    plot(log10(ISIStats.summstats.(states{ss}).meanrate),rateLFPcorr.(states{ss}),'.')
    LogScale('x',10)
    xlabel('sFR (Hz)');ylabel('LFP-Rate Corr')
    
% subplot(2,2,3)
%     plot(rateLFPcorr.WAKEstate,rateLFPcorr.NREMstate,'.')
%     hold on
%     axis tight
%     UnityLine
%     xlabel('WAKE LFP-Rate Corr');ylabel('NREM LFP-Rate Corr')  
    
    NiceSave('CellLFPStats',figfolder,baseName)

%%
% cellex = randsample(spikes.numcells,1);
% [ CONDXY ] = ConditionalHist(ISIStats.allspikes.LFPval{cellex},...
%     log10(1./ISIStats.allspikes.meanISI{cellex}),'numXbins',ratebyLFP.nbins);
% 
% 
% figure
% subplot(2,2,1)
% plot(ISIStats.allspikes.LFPval{cellex},log10(1./ISIStats.allspikes.meanISI{cellex}),'.')
% subplot(2,2,2)
% imagesc(CONDXY.Xbins,CONDXY.Ybins,CONDXY.pYX')
% hold on
% plot(ratebyLFP.lfpbins,log10(ratebyLFP.rate{cellex}),'w')
% plot(ratebyLFP.lfpbins([1 end]),log10(ISIStats.summstats.ALL.meanrate(cellex).*[1 1]),'r')
% axis tight
% axis xy
% LogScale('y',10)
% subplot(2,2,3)
% plot(ratebyLFP.lfpbins,log10(ratebyLFP.rate{cellex}))
% 
% subplot(2,2,4)
% imagesc(log10(cat(1,ratebyLFP.rate{:})))




%% GLM for coupling
%excell = 9;
thisstate = 'ALL';
stateints = SleepState.ints.(thisstate);



clear GLMmodelfit
for cc = 1:spikes.numcells
    cc
    [ GLMmodelfit(cc) ] = GLMLFP_raw(spikes.times(cc),lfp,...
        'intervals',stateints );
end





%% Simulate spikes from GLM
clear Poissmodelfit
clear simspikes
clear simspikes_poiss
for cc = 1:spikes.numcells
    cc
    GLMspkmat = rand(size(GLMmodelfit(cc).predRate))<=GLMmodelfit(cc).predRate;
    Poissspkmat  = rand(size(GLMmodelfit(cc).predRate))<=ISIStats.summstats.(thisstate).meanrate(cc).*GLMmodelfit(cc).dt;
%     for tt = 1:length(GLMmodelfit(cc).timestamps)
%         GLMspkmat(tt) = rand(1)<=GLMmodelfit(cc).predRate(tt);
%         Poissspkmat(tt) = rand(1)<=ISIStats.summstats.NREMstate.meanrate(cc).*GLMmodelfit(cc).dt;
%     end
    simspikes.times{cc} = GLMmodelfit(cc).timestamps(GLMspkmat);
    simspikes_poiss.times{cc} = GLMmodelfit(cc).timestamps(Poissspkmat);
    clear GLMspkmat
end
simspikes.UID = spikes.UID;
simspikes_poiss.UID = spikes.UID;

%%
viewwin = bz_RandomWindowInIntervals(SleepState.ints.NREMstate,5);
figure
subplot(5,1,4)
    hold on
    for cc = 1:spikes.numcells
        whichcell = ISIStats.sorts.(thisstate).rate(cc);
        plot(simspikes.times{whichcell},cc.*ones(size(simspikes.times{whichcell})),'k.')
    end
    xlim(viewwin);ylim([0 spikes.numcells])
    ylabel({'Simulated', 'Spikes'});
    box off
    set(gca,'xticklabels',[])
    
subplot(5,1,5)
    hold on
    for cc = 1:spikes.numcells
        whichcell = ISIStats.sorts.(thisstate).rate(cc);
        plot(simspikes_poiss.times{whichcell},cc.*ones(size(simspikes_poiss.times{whichcell})),'k.')
    end
    xlim(viewwin);ylim([0 spikes.numcells])
    ylabel({'Poisson','Spikes'});
    box off
    set(gca,'xticklabels',[])
%  

% subplot(5,1,3)
% plot([GLMmodelfit(:).timestamps],log([GLMmodelfit(:).predRate]))
% xlim(viewwin)

subplot(5,1,2)
    bz_MultiLFPPlot(lfp,'timewin',viewwin)
    
subplot(5,1,1)
    hold on
    for cc = 1:spikes.numcells
        whichcell = ISIStats.sorts.(thisstate).rate(cc);
        plot(spikes.times{whichcell},cc.*ones(size(spikes.times{whichcell})),'k.')
    end
%     for tt = 1:length(celltypes)
%         plot(spikemat.timestamps,spikemat.poprate.(celltypes{tt}),cellcolor{tt})
%     end
    xlim(viewwin);ylim([0 spikes.numcells])
    ylabel({'Observed', 'Spikes'});
    box off
    set(gca,'xticklabels',[])
    
 NiceSave('SimPopLFP',figfolder,baseName)

%% Calculate ISI stats for simulated spikes
[ ISIstats_sim ] = bz_ISIStats( simspikes,'ints',SleepState.ints,...
    'cellclass',CellClass.label);%,'figfolder',figfolder);

[ ISIstats_poiss ] = bz_ISIStats( simspikes_poiss,'ints',SleepState.ints,...
    'cellclass',CellClass.label);%,'figfolder',figfolder);

%Mean ISI histograms
for tt= 1:length(celltypes)
    ISIstats_sim.ISIhist.NREMstate.popmean.(celltypes{tt}) = ...
        mean(ISIstats_sim.ISIhist.NREMstate.log(CellClass.(celltypes{tt}),:),1);
    ISIstats_poiss.ISIhist.NREMstate.popmean.(celltypes{tt}) = ...
        mean(ISIstats_poiss.ISIhist.NREMstate.log(CellClass.(celltypes{tt}),:),1);
    ISIStats.ISIhist.NREMstate.popmean.(celltypes{tt}) = ...
        mean(ISIStats.ISIhist.NREMstate.log(CellClass.(celltypes{tt}),:),1);
end



%%
figure
subplot(3,3,1)
    plot(log10(ISIStats.summstats.(thisstate).meanrate),log10(ISIstats_sim.summstats.(thisstate).meanrate),'.')
subplot(3,3,2)
    plot(ISIStats.summstats.(thisstate).meanCV2,ISIstats_sim.summstats.(thisstate).meanCV2,'.')
    hold on;box off
        plot([0.8 1.4],[0.8 1.4],'k')
        plot([0.8 1.4],[1 1],'k:')
        xlabel('CV2 - Observed');ylabel('CV2 - Simulated')

subplot(3,3,3)
plot(log10(ISIStats.summstats.(thisstate).meanrate),ISIstats_sim.summstats.(thisstate).meanCV2,'.')
hold on;box off
LogScale('x',10)
    plot(get(gca,'xlim'),[1 1],'k:')
    xlabel('FR');ylabel('CV2')
% subplot(2,2,4)
% plot(log10(ISIStats.summstats.(thisstate).meanrate),ISIstats_poiss.summstats.(thisstate).meanCV2,'.')
% hold on
% LogScale('x',10)
%     plot(get(gca,'xlim'),[1 1],'k:')


subplot(3,3,4)
plot(log10(ISIStats.summstats.(thisstate).meanrate),...
    (ISIstats_sim.summstats.(thisstate).meanCV2-ISIStats.summstats.(thisstate).meanCV2)./...
    (ISIStats.summstats.(thisstate).meanCV2-1),'.')
hold on
    plot(get(gca,'xlim'),[0 0],'k-')
    plot(get(gca,'xlim'),[1 1],'k:')
    plot(get(gca,'xlim'),[-1 -1],'k:')

LogScale('x',10)
    xlabel('FR');ylabel('FracCV2')


for tt = 1:length(celltypes)    
subplot(4,3,9+tt)
    plot(ISIstats_sim.ISIhist.logbins,...
        ISIstats_sim.ISIhist.NREMstate.popmean.(celltypes{tt}),'--','color',cellcolor{tt})
    hold on
    plot(ISIstats_poiss.ISIhist.logbins,...
        ISIstats_poiss.ISIhist.NREMstate.popmean.(celltypes{tt}),':','color',cellcolor{tt})
    plot(ISIStats.ISIhist.logbins,...
        ISIStats.ISIhist.NREMstate.popmean.(celltypes{tt}),'-','linewidth',2,'color',cellcolor{tt})
    xlabel('ISI (s)');
    axis tight
    box off
    LogScale('x',10)
    
end

NiceSave('GLMSimCV2',figfolder,baseName)


%%
ISIstats_sim.allspikes.LFPval = cellfun(@(X) interp1(lfp.timestamps,single(lfp.data),X),...
    ISIstats_sim.allspikes.times,'UniformOutput',false);

%%
spikemat = bz_SpktToSpkmat(simspikes,'binsize',binsize);
spikemat.lfp = interp1(lfp.timestamps,single(lfp.data),spikemat.timestamps);

for ss = 1:length(states)
    spikemat.timeidx.(states{ss}) = InIntervals(spikemat.timestamps,SleepState.ints.(states{ss}));
    rateLFPcorr_sim.(states{ss}) = corr(spikemat.lfp(spikemat.timeidx.(states{ss})),...
        spikemat.data(spikemat.timeidx.(states{ss}),:),'type','spearman','rows','complete');
end

%%

for ss = 1:4
    instatelfp = InIntervals(lfp.timestamps,SleepState.ints.(states{ss}));
    instatespikes = cellfun(@(X) InIntervals(X,SleepState.ints.(states{ss})),...
        ISIstats_sim.allspikes.times,'UniformOutput',false);
    
    ratebyLFP_sim.(states{ss}).lfphist = hist(single(lfp.data(instatelfp)),ratebyLFP.lfpbins);
    ratebyLFP_sim.(states{ss}).lfphist = ratebyLFP_sim.(states{ss}).lfphist./lfp.samplingRate;
    ratebyLFP_sim.spikecounts = cellfun(@(X,Y) hist(X(Y),ratebyLFP.lfpbins),...
        ISIstats_sim.allspikes.LFPval,instatespikes,'UniformOutput',false);
    ratebyLFP_sim.(states{ss}).rate = cellfun(@(X) X./ratebyLFP_sim.(states{ss}).lfphist,...
        ratebyLFP_sim.spikecounts,'UniformOutput',false);
    ratebyLFP_sim.(states{ss}).rate = cat(1,ratebyLFP_sim.(states{ss}).rate{:});
    ratebyLFP_sim.(states{ss}).rate(:,ratebyLFP_sim.(states{ss}).lfphist<mintime) = nan;
    
    ratebyLFP_sim.(states{ss}).lfphist = ratebyLFP_sim.(states{ss}).lfphist./lfp.timestamps(end);
    
    for cc = 1:length(celltypes)
        ratebyLFP_sim.(states{ss}).meanrate.(celltypes{cc}) = ...
            nanmean(ratebyLFP_sim.(states{ss}).rate(CellClass.(celltypes{cc}),:),1);
        ratebyLFP_sim.(states{ss}).std.(celltypes{cc}) = ...
            nanstd(ratebyLFP_sim.(states{ss}).rate(CellClass.(celltypes{cc}),:),[],1);
    end
end
clear instatelfp


%%
figure

subplot(6,4,12)
    hold on
    for ss = 1:3
    plot(ratebyLFP.lfpbins,ratebyLFP_sim.(states{ss}).lfphist,...
        'linewidth',2,'color',statecolors{ss})
    end
    xlabel('LFP (Z_N_R_E_M)');       axis tight
    xlim(ratebyLFP.lfpbins([1 end]))
    
for cc = 1:length(celltypes)
subplot(6,4,cc*4)
    hold on
%     for ss = 1:3
%     errorshade(ratebyLFP.lfpbins,(ratebyLFP.(states{ss}).meanrate),...
%         (ratebyLFP.(states{ss}).std),(ratebyLFP.(states{ss}).std),...
%         statecolors{ss},'scalar')
%     end
    for ss = 4:-1:1
    plot(ratebyLFP.lfpbins,(ratebyLFP_sim.(states{ss}).meanrate.(celltypes{cc})),...
        'linewidth',2,'color',statecolors{ss})
    end
    axis tight
    xlim(ratebyLFP.lfpbins([1 end]));
    set(gca,'xtick',[]);  ylabel({celltypes{cc},'Rate (Hz)'})    
   % LogScale('y',10)
end

subplot(3,3,1)
    hold on
    for ss = 1:3
        [histcounts,histbins] = hist(rateLFPcorr_sim.(states{ss}));
        plot(histbins,histcounts,'color',statecolors{ss},'linewidth',2)
    end
    xlabel('LFP-Rate Corr');ylabel('# Cells')
    
    
subplot(6,4,3)
    for tt = 1:length(celltypes)
%         plot(PSScorrhist.bins,PSScorrhist.CV2.(states{ss}).(cellclasses{tt}),...
%             classcolors{tt},'linewidth',2)
            HistWithMean(rateLFPcorr_sim.ALL(CellClass.(celltypes{tt})),...
                'numbins',8,'color',cellcolor{tt},'showtext',false)
        hold on
    end
    plot([0 0],get(gca,'ylim'),'k')
    xlabel('LFP-Rate Corr');
    ylabel('# cells')
    xlim([-0.1 0.1])
    %legend('pE cells','pI ce
    
subplot(3,3,4)
    plot(log10(ISIStats.summstats.(states{ss}).meanrate),rateLFPcorr_sim.(states{ss}),'.')
    LogScale('x',10)
    xlabel('sFR (Hz)');ylabel('LFP-Rate Corr')
    
% subplot(2,2,3)
%     plot(rateLFPcorr.WAKEstate,rateLFPcorr.NREMstate,'.')
%     hold on
%     axis tight
%     UnityLine
%     xlabel('WAKE LFP-Rate Corr');ylabel('NREM LFP-Rate Corr')  
    
    NiceSave('CellLFPStats_sim',figfolder,baseName)

%%

[specslope,specgram] = bz_PowerSpectrumSlope([],[],[],'saveMat',basePath);

numbins = 20;
PSShist.bins = linspace(-2,0,numbins);
for ss = 1:length(states)
    specslope.timeidx.(states{ss}) = InIntervals(specslope.timestamps,SleepState.ints.(states{ss}));
    
    PSShist.(states{ss}) = hist(specslope.data(specslope.timeidx.(states{ss})),PSShist.bins);
end

ISIStats.allspikes.PSS = cellfun(@(X) interp1(specslope.timestamps,specslope.data,X,'nearest'),ISIStats.allspikes.times,'UniformOutput',false);
ISIstats_sim.allspikes.PSS = cellfun(@(X) interp1(specslope.timestamps,specslope.data,X,'nearest'),ISIstats_sim.allspikes.times,'UniformOutput',false);
ISIstats_poiss.allspikes.PSS = cellfun(@(X) interp1(specslope.timestamps,specslope.data,X,'nearest'),ISIstats_poiss.allspikes.times,'UniformOutput',false);

CV2PSScorr.ALL = cellfun(@(X,Y) corr(X,Y,'type','spearman','rows','complete'),ISIStats.allspikes.PSS,ISIStats.allspikes.CV2);
CV2PSScorr.ALL_sim = cellfun(@(X,Y) corr(X,Y,'type','spearman','rows','complete'),ISIstats_sim.allspikes.PSS,ISIstats_sim.allspikes.CV2);
CV2PSScorr.ALL_poiss = cellfun(@(X,Y) corr(X,Y,'type','spearman','rows','complete'),ISIstats_poiss.allspikes.PSS,ISIstats_poiss.allspikes.CV2);

%% Mean binned CV2...
clear CV2mat
CV2mat.winsize = specslope.detectionparms.winsize;
CV2mat.timestamps = specslope.timestamps;
CV2mat.binedges = bsxfun(@(X,Y) X+Y,specslope.timestamps,[-0.5 0.5].*CV2mat.winsize);
for tt = 1:length(celltypes)
    allspikes.CV2.(celltypes{tt}) = cat(1,ISIstats_sim.allspikes.CV2{CellClass.(celltypes{tt})});
    allspikes.times.(celltypes{tt}) = cat(1,ISIstats_sim.allspikes.times{CellClass.(celltypes{tt})});
    [CV2mat.timestamps,CV2mat.sim.(celltypes{tt})] = ...
        BinDataTimes(allspikes.CV2.(celltypes{tt}),allspikes.times.(celltypes{tt}),CV2mat.binedges);
    
    allspikes.CV2.(celltypes{tt}) = cat(1,ISIStats.allspikes.CV2{CellClass.(celltypes{tt})});
    allspikes.times.(celltypes{tt}) = cat(1,ISIStats.allspikes.times{CellClass.(celltypes{tt})});
    [CV2mat.timestamps,CV2mat.(celltypes{tt})] = ...
        BinDataTimes(allspikes.CV2.(celltypes{tt}),allspikes.times.(celltypes{tt}),CV2mat.binedges);
end
CV2mat.PSS = interp1(specslope.timestamps,specslope.data,CV2mat.timestamps);
%%
numXbins = 60;
numYbins = 150;

Xbounds = [-2 0];
Ybounds = [0 2];

minX = 50;

[ PSSpEhist ] = ConditionalHist(CV2mat.PSS,CV2mat.sim.pE,...
    'numXbins',60,'numYbins',150,'Xbounds',[-2 0],'Ybounds',[0 2]);
[ PSSpIhist ] = ConditionalHist(CV2mat.PSS,CV2mat.sim.pI,...
    'numXbins',60,'numYbins',150,'Xbounds',[-2 0],'Ybounds',[0 2]);
%%
figure
subplot(5,4,1)
imagesc(PSSpEhist.Xbins,PSSpEhist.Ybins,PSSpEhist.pYX')
axis xy
hold on
xlim([-1.6 -0.3])
ylabel({'CV_2', 'pE Pop.'})
ylim([0.9 1.4])
plot(get(gca,'xlim'),[1 1],'w--')

subplot(5,4,5)
imagesc(PSSpIhist.Xbins,PSSpIhist.Ybins,PSSpIhist.pYX')
axis xy
hold on
xlim([-1.6 -0.3])
ylim([0.75 1.15])
ylabel({'CV_2',' pI Pop.'})
plot(get(gca,'xlim'),[1 1],'w--')

subplot(10,4,17)
    for ss = 1:3
    plot(PSShist.bins,PSShist.(states{ss}),'color',statecolors{ss},'linewidth',2)
    hold on
    end
    xlabel('PSS')
    ylabel({'Time', 'Occupancy'})
    set(gca,'ytick',[])
    %legend(states{:},'location','eastoutside')
    axis tight
    box off
    xlim([-1.6 -0.3])
    
subplot(6,4,2)
    for tt = 1:length(celltypes)
%         plot(PSScorrhist.bins,PSScorrhist.CV2.(states{ss}).(cellclasses{tt}),...
%             classcolors{tt},'linewidth',2)
            HistWithMean(CV2PSScorr.ALL_sim(CellClass.(celltypes{tt})),...
                'numbins',8,'color',cellcolor{tt},'showtext',false)
        hold on
    end
    plot([0 0],get(gca,'ylim'),'k')
    xlabel('PSS-CV2 Corr');
    ylabel('# cells')
    xlim([-0.25 0.25])
    title({'PSS-CV2 Corr','LFP GLM'})
    
%subplot(
    
    NiceSave('simCV2byPSSstats',figfolder,baseName,'tiff')
    
    
    %%

figure

for ss = 1:3
    exwin = randsample(Restrict(CV2mat.timestamps,SleepState.ints.(states{ss})),1);


    subplot(3,3,0+ss)
        bz_MultiLFPPlot(lfp,'spikes',spikes,'timewin',...
            exwin+specslope.detectionparms.winsize.*[-0.2 0.2],...
            'cellgroups',{CellClass.pE,CellClass.pI},...
            'sortmetric',ISIStats.summstats.ALL.meanrate,...
            'scaleLFP',0.4,'scalespikes',0.1,'spikeside','bottom')
        xlabel('')
        box off
        bz_ScaleBar('s')
        title({['PSS: ',num2str(round(CV2mat.PSS(CV2mat.timestamps==exwin),3))],...
            ['CV2pE: ',num2str(round(CV2mat.pE(CV2mat.timestamps==exwin),3))],...
            ['CV2pI: ',num2str(round(CV2mat.pI(CV2mat.timestamps==exwin),3))]})
        set(gca,'ytick',[]);ylabel('Observed');

    subplot(6,3,9+ss)
        hold on
        for cc = 1:spikes.numcells
            whichcell = ISIStats.sorts.ALL.rate(cc);
            plot(simspikes.times{whichcell},cc.*ones(size(simspikes.times{whichcell})),...
                'k.','markersize',0.1)
        end
        xlim(exwin+specslope.detectionparms.winsize.*[-0.2 0.2]);ylim([0 spikes.numcells])
        ylabel('LFP GLM');
        title({['CV2pE: ',num2str(round(CV2mat.sim.pE(CV2mat.timestamps==exwin),3))],...
            ['CV2pI: ',num2str(round(CV2mat.sim.pI(CV2mat.timestamps==exwin),3))]})    
        box off
        set(gca,'xticklabels',[]);set(gca,'ytick',[])

    subplot(6,3,12+ss)
        hold on
        for cc = 1:spikes.numcells
            whichcell = ISIStats.sorts.ALL.rate(cc);
            plot(simspikes_poiss.times{whichcell},cc.*ones(size(simspikes_poiss.times{whichcell})),...
                'k.','markersize',0.1)
        end
        xlim(exwin+specslope.detectionparms.winsize.*[-0.2 0.2]);ylim([0 spikes.numcells])
        ylabel('Poisson');
        box off
        set(gca,'xticklabels',[]);set(gca,'ytick',[])
end
    NiceSave('PSSandSpikingExample2',figfolder,baseName,'tiff')
    
%%
figure
subplot(2,2,1)
hold on
for tt = 1:length(celltypes)
plot(log10(ISIStats.summstats.(thisstate).meanrate(CellClass.(celltypes{tt}))),...
    [GLMmodelfit(CellClass.(celltypes{tt})).Rlfp],'.','color',cellcolor{tt})
LogScale('x',10);
xlabel('sFR (Hz)');ylabel('LFP-coupling')

end

subplot(2,2,2)
hold on
for tt = 1:length(celltypes)
plot((ISIStats.summstats.(thisstate).meanCV2(CellClass.(celltypes{tt}))),...
    [GLMmodelfit(CellClass.(celltypes{tt})).Rlfp],'.','color',cellcolor{tt})
xlabel('<CV2>');ylabel('LFP-coupling')

end


%%

plotstates = {'WAKEstate','REMstate','WAKEstate'};
plotstates2 = {'NREMstate','NREMstate','REMstate'};

figure
figure

%Rate
for ss=1:3
    subplot(4,4,ss)
        plot(log10(ISIstats_sim.summstats.(plotstates{ss}).meanrate(CellClass.pE)),...
            log10(ISIstats_sim.summstats.(plotstates2{ss}).meanrate(CellClass.pE)),...
            'k.','markersize',2)
        hold on
        plot(log10(ISIstats_sim.summstats.(plotstates{ss}).meanrate(CellClass.pI)),...
            log10(ISIstats_sim.summstats.(plotstates2{ss}).meanrate(CellClass.pI)),...
            'r.','markersize',2)
        %plot(log10([0.03 30]),log10([0.03 30]),'k')
        xlabel([plotstates{ss},' Rate']);ylabel([plotstates2{ss},' Rate'])
        axis tight
        LogScale('xy',10)
        UnityLine('linetype','-')
end

%CV
for ss=1:3
    subplot(4,4,ss+4)
        plot(log2(ISIstats_sim.summstats.(plotstates{ss}).ISICV(CellClass.pE)),...
            log2(ISIstats_sim.summstats.(plotstates2{ss}).ISICV(CellClass.pE)),...
            'k.','markersize',2)
        hold on
        plot(log2(ISIstats_sim.summstats.(plotstates{ss}).ISICV(CellClass.pI)),...
            log2(ISIstats_sim.summstats.(plotstates2{ss}).ISICV(CellClass.pI)),...
            'r.','markersize',2)
        %plot(log2([1 6]),log2([1 6]),'k')
        xlabel([plotstates{ss},' CV']);ylabel([plotstates2{ss},' CV'])
        axis tight
        LogScale('xy',2)
        UnityLine('linetype','-')
end

for ss=1:3
    subplot(4,4,ss+8)
        plot((ISIstats_sim.summstats.(plotstates{ss}).meanCV2(CellClass.pE)),...
            (ISIstats_sim.summstats.(plotstates2{ss}).meanCV2(CellClass.pE)),...
            'k.','markersize',2)
        hold on
        plot((ISIstats_sim.summstats.(plotstates{ss}).meanCV2(CellClass.pI)),...
            (ISIstats_sim.summstats.(plotstates2{ss}).meanCV2(CellClass.pI)),...
            'r.','markersize',2)
        %plot(([0 2]),([0 2]),'k')
        
        xlabel([plotstates{ss},' CV2']);ylabel([plotstates2{ss},' CV2'])
        xlim([0.4 1.6]);ylim([0.4 1.6])
       % LogScale('xy',2)
       axis tight
        UnityLine('linetype','-')
end
NiceSave('ISIstatsbystate_sim',figfolder,baseName)

%%

figure
for ss = 1:3
subplot(3,3,ss)
    plot(log10(ISIstats_sim.summstats.(states{ss}).meanrate(CellClass.pE)),...
        ISIstats_sim.summstats.(states{ss}).meanCV2(CellClass.pE),'k.','markersize',4)
    hold on
    plot(log10(ISIstats_sim.summstats.(states{ss}).meanrate(CellClass.pI)),...
        ISIstats_sim.summstats.(states{ss}).meanCV2(CellClass.pI),'r.','markersize',4)
    xlim([-2.2 1.7]); ylim([0.9 1.1])
	LogScale('x',10)
    plot(get(gca,'xlim'),[1 1],'k')
    title(states{ss})
    xlabel('FR (Hz)');ylabel('<CV2>')
    
    
subplot(3,3,ss+3)
    plot(log10(ISIstats_sim.summstats.(states{ss}).meanrate(CellClass.pE)),...
        log2(ISIstats_sim.summstats.(states{ss}).ISICV(CellClass.pE)),'k.','markersize',4)
    hold on
    plot(log10(ISIstats_sim.summstats.(states{ss}).meanrate(CellClass.pI)),...
        log2(ISIstats_sim.summstats.(states{ss}).ISICV(CellClass.pI)),'r.','markersize',4)
    
    xlim([-2.2 1.7]); ylim([-0.25 0.25])
    LogScale('x',10);LogScale('y',2);
    plot(get(gca,'xlim'),[0 0],'k')
    title(states{ss})
    xlabel('FR (Hz)');ylabel('CV')
end

NiceSave('RateandCV2',figfolder,baseName)
