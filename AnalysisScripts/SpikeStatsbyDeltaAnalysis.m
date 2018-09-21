function [  ] = SpikeStatsbyDeltaAnalysis(basePath,figfolder)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% DEV
reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
basePath = '/mnt/proraidDL/Database/BWCRCNS/JennBuzsaki22/20140526_277um';
%basePath = pwd;
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/SpikeStatsbyDeltaAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
SleepState.ints.ALL = [0 Inf];
statenames = fieldnames(SleepState.ints);
numstates = length(statenames);
statecolors = {'k','b','r'};
classnames = unique(CellClass.label);
numclasses = length(classnames);
classcolors = {'k','r'};
SlowWaves = bz_LoadEvents(basePath,'SlowWaves');

downsamplefactor = 1;
lfp = bz_GetLFP(SlowWaves.detectorinfo.detectionchannel,...
     'basepath',basePath,'downsample',downsamplefactor);
sessionInfo = bz_getSessionInfo(basePath);

%% Filter the LFP for delta activity

deLFP = bz_Filter(lfp,'passband',[0.5 6],'order',1,'filter','fir1');

%gaLFP = bz_Filter(lfp,'passband',[50 120],'order',4,'filter','fir1');

%% Calculate the power-phase ratemap, cv2map

[PowerPhaseRatemap,spikebinIDs] = bz_PowerPhaseRatemap(ISIStats.allspikes,deLFP,...
    'ints',SleepState.ints.NREMstate,'metric',ISIStats.allspikes.CV2);

%%
% PowerPhaseRatemap_ga = bz_PowerPhaseRatemap(ISIStats.allspikes,gaLFP,...
%     'ints',SleepState.ints.NREMstate,'metric',ISIStats.allspikes.CV2);

%CV2map
%%
figure
plot(PowerPhaseRatemap.phasebins, PowerPhaseRatemap.meanrate(1,:))
hold on
plot(PowerPhaseRatemap.phasebins+2*pi, PowerPhaseRatemap.meanrate(1,:))
plot(PowerPhaseRatemap.phasebins, PowerPhaseRatemap.meanrate(end,:))
hold on
plot(PowerPhaseRatemap.phasebins+2*pi, PowerPhaseRatemap.meanrate(end,:))
plot(PowerPhaseRatemap.phasebins, PowerPhaseRatemap.meanrate(end./2,:))
hold on
plot(PowerPhaseRatemap.phasebins+2*pi, PowerPhaseRatemap.meanrate(end/2,:))
%% GLM for coupling
excell = 9;

predLFP = deLFP;
predLFP.data = predLFP.hilb;
predLFP.freqs = 4;
clear GLMmodelfit
for cc = 1:spikes.numcells
    cc
% [ GLMmodelfit(cc) ] = GLMLFP_param(spikes.times(cc),predLFP,...
%     'intervals',SleepState.ints.NREMstate );
    [ GLMmodelfit(cc) ] = GLMLFP(spikes.times(cc),predLFP,...
        'intervals',SleepState.ints.NREMstate );
end
%% 

%% Nonparametric coupling for comparison
% GLMmodel_nonparam = GLMLFP(spikes.times(cc),predLFP,...
%     'intervals',SleepState.ints.NREMstate );
%% Take a look at some parameters
figure
plot(log10(ISIStats.summstats.NREMstate.meanrate),log10([GLMmodelfit.R0]),'.')

%% Simulate spikes from GLM
for cc = 1:spikes.numcells
    cc
    for tt = 1:length(GLMmodelfit(cc).timestamps)
        GLMmodelfit(cc).spkmat(tt) = rand(1)<=GLMmodelfit(cc).predRate(tt);
        Poissmodelfit(cc).spkmat(tt) = rand(1)<=GLMmodelfit(excell).R0;
    end
    simspikes.times{cc} = GLMmodelfit(cc).timestamps(GLMmodelfit(cc).spkmat);
    simspikes_poiss.times{cc} = GLMmodelfit(cc).timestamps(Poissmodelfit(cc).spkmat);
end
simspikes.UID = spikes.UID;
simspikes_poiss.UID = spikes.UID;
%% Calculate ISI stats for simulated spikes
[ ISIstats_sim ] = bz_ISIStats( simspikes,'ints',SleepState.ints,...
    'cellclass',CellClass.label);%,'figfolder',figfolder);

[ ISIstats_poiss ] = bz_ISIStats( simspikes_poiss,'ints',SleepState.ints,...
    'cellclass',CellClass.label);%,'figfolder',figfolder);
%% Calculate the power-phase ratemap, cv2map of simulated spikes

[PowerPhaseRatemap_sim] = bz_PowerPhaseRatemap(ISIstats_sim.allspikes,deLFP,...
    'ints',SleepState.ints.NREMstate,'metric',ISIstats_sim.allspikes.CV2);
%%
for tt= 1:numclasses
    PowerPhaseRatemap_sim.classCV2.(classnames{tt}) = nanmean(cat(3,PowerPhaseRatemap_sim.metricmap{CellClass.(classnames{tt})}),3);
    PowerPhaseRatemap_sim.classrate.(classnames{tt}) = nanmean(cat(3,PowerPhaseRatemap_sim.ratemap{CellClass.(classnames{tt})}),3);

    PowerPhaseRatemap.classCV2.(classnames{tt}) = nanmean(cat(3,PowerPhaseRatemap.metricmap{CellClass.(classnames{tt})}),3);
    PowerPhaseRatemap.classrate.(classnames{tt}) = nanmean(cat(3,PowerPhaseRatemap.ratemap{CellClass.(classnames{tt})}),3);

  %  PowerPhaseRatemap_ga.classCV2.(classnames{tt}) = nanmean(cat(3,PowerPhaseRatemap_ga.metricmap{CellClass.(classnames{tt})}),3);
   % PowerPhaseRatemap_ga.classrate.(classnames{tt}) = nanmean(cat(3,PowerPhaseRatemap_ga.ratemap{CellClass.(classnames{tt})}),3);

    
    ISIstats_sim.ISIhist.NREMstate.popmean.(classnames{tt}) = ...
        mean(ISIstats_sim.ISIhist.NREMstate.log(CellClass.(classnames{tt}),:),1);
    ISIstats_poiss.ISIhist.NREMstate.popmean.(classnames{tt}) = ...
        mean(ISIstats_poiss.ISIhist.NREMstate.log(CellClass.(classnames{tt}),:),1);
    ISIStats.ISIhist.NREMstate.popmean.(classnames{tt}) = ...
        mean(ISIStats.ISIhist.NREMstate.log(CellClass.(classnames{tt}),:),1);
end

%%
figure
subplot(3,3,1)
    for tt= 1:numclasses
    plot(ISIStats.summstats.NREMstate.meanCV2(CellClass.(classnames{tt})),...
        ISIstats_sim.summstats.NREMstate.meanCV2(CellClass.(classnames{tt})),'.','color',classcolors{tt})
    hold on
    end
    plot([0.8 1.4],[0.8 1.4],'k')
    plot([0.8 1.4],[1 1],'k:')
    xlabel('<CV2> observed');ylabel('<CV2> delta model');
    
subplot(3,3,2)
    for tt= 1:numclasses
    plot(log10(ISIStats.summstats.NREMstate.meanrate(CellClass.(classnames{tt}))),...
        ISIstats_sim.summstats.NREMstate.meanCV2(CellClass.(classnames{tt})),'.','color',classcolors{tt})
    hold on
    end
    plot(get(gca,'xlim'),[1 1],'k:')
    xlabel('Rate');ylabel('<CV2> delta model');
    LogScale('x',10)
 
for tt = 1:numclasses    
subplot(4,3,6+tt)
    plot(ISIstats_sim.ISIhist.logbins,...
        ISIstats_sim.ISIhist.NREMstate.popmean.(classnames{tt}),'--','color',classcolors{tt})
    hold on
%    plot(ISIstats_poiss.ISIhist.logbins,...
%        ISIstats_poiss.ISIhist.NREMstate.popmean.(classnames{tt}),':','color',classcolors{tt})
    plot(ISIStats.ISIhist.logbins,...
        ISIStats.ISIhist.NREMstate.popmean.(classnames{tt}),'-','linewidth',2,'color',classcolors{tt})
    xlabel('ISI (s)');
    axis tight
    box off
    LogScale('x',10)
    
end
NiceSave('PopStats',figfolder,baseName)

%%
cv2color = [makeColorMap([0.5 0.5 1],[0 0 0.8],[0 0 0]);makeColorMap([0 0 0],[0.8 0 0],[1 0.5 0.5])];

maxrate.pE = 2.5;
maxrate.pI = 20;

figure
for tt= 1:numclasses
    subplot(4,4,tt)
        imagesc(PowerPhaseRatemap.phasebins,PowerPhaseRatemap.powerbins,...
            PowerPhaseRatemap_sim.classrate.(classnames{tt}))
        hold on
        imagesc(PowerPhaseRatemap.phasebins+2*pi,PowerPhaseRatemap.powerbins,...
            PowerPhaseRatemap_sim.classrate.(classnames{tt}))
        xlim([-pi 3*pi])
        caxis([0 maxrate.(classnames{tt})])
        axis xy
        colorbar
        xlabel('Phase');ylabel('Power')
        title('Simulated Delta Ratemap')
end

for tt= 1:numclasses
    subplot(4,4,2+tt)
        h=imagesc(PowerPhaseRatemap.phasebins,PowerPhaseRatemap.powerbins,...
            PowerPhaseRatemap_sim.classCV2.(classnames{tt}))
        hold on
        h2=imagesc(PowerPhaseRatemap.phasebins+2*pi,PowerPhaseRatemap.powerbins,...
            PowerPhaseRatemap_sim.classCV2.(classnames{tt}))
        set(h,'AlphaData',~isnan(PowerPhaseRatemap_sim.classCV2.(classnames{tt})));
        set(h2,'AlphaData',~isnan(PowerPhaseRatemap_sim.classCV2.(classnames{tt})));
        xlim([-pi 3*pi])
        caxis([0.7 1.3])
        axis xy
        colorbar
        colormap(gca,cv2color)
        xlabel('Phase');ylabel('Power')
        title('Simulated Delta CV2map')
end

for tt= 1:numclasses
    subplot(4,4,4+tt)
        imagesc(PowerPhaseRatemap.phasebins,PowerPhaseRatemap.powerbins,...
            PowerPhaseRatemap.classrate.(classnames{tt}))
        hold on
        imagesc(PowerPhaseRatemap.phasebins+2*pi,PowerPhaseRatemap.powerbins,...
            PowerPhaseRatemap.classrate.(classnames{tt}))
        xlim([-pi 3*pi])
        caxis([0 maxrate.(classnames{tt})])
        axis xy
        colorbar
        xlabel('Phase');ylabel('Power')
        title('Real Delta Ratemap')
end

for tt= 1:numclasses
    subplot(4,4,6+tt)
        h=imagesc(PowerPhaseRatemap.phasebins,PowerPhaseRatemap.powerbins,...
            PowerPhaseRatemap.classCV2.(classnames{tt}))
        hold on
        h2=imagesc(PowerPhaseRatemap.phasebins+2*pi,PowerPhaseRatemap.powerbins,...
            PowerPhaseRatemap.classCV2.(classnames{tt}))
        set(h,'AlphaData',~isnan(PowerPhaseRatemap.classCV2.(classnames{tt})));
        set(h2,'AlphaData',~isnan(PowerPhaseRatemap.classCV2.(classnames{tt})));
        xlim([-pi 3*pi])
        caxis([0.7 1.3])
        axis xy
        colorbar
        colormap(gca,cv2color)
        xlabel('Phase');ylabel('Power')
        title('Real Delta CV2map')
end


% for tt= 1:numclasses
%     subplot(4,4,12+tt)
%         imagesc(PowerPhaseRatemap_ga.phasebins,PowerPhaseRatemap_ga.powerbins,...
%             PowerPhaseRatemap_ga.classrate.(classnames{tt}))
%         hold on
%         imagesc(PowerPhaseRatemap_ga.phasebins+2*pi,PowerPhaseRatemap_ga.powerbins,...
%             PowerPhaseRatemap_ga.classrate.(classnames{tt}))
%         xlim([-pi 3*pi])
%         caxis([0 maxrate.(classnames{tt})])
%         axis xy
%         colorbar
%         xlabel('Phase');ylabel('Power')
%         title('Real Gamma Ratemap')
% end


% for tt= 1:numclasses
%     subplot(4,4,14+tt)
%         imagesc(PowerPhaseRatemap_ga.phasebins,PowerPhaseRatemap_ga.powerbins,...
%             PowerPhaseRatemap_ga.classCV2.(classnames{tt}))
%         hold on
%         imagesc(PowerPhaseRatemap_ga.phasebins+2*pi,PowerPhaseRatemap_ga.powerbins,...
%             PowerPhaseRatemap_ga.classCV2.(classnames{tt}))
%         xlim([-pi 3*pi])
%         caxis([0.7 1.3])
%         axis xy
%         colorbar
%         colormap(gca,cv2color)
%         xlabel('Phase');ylabel('Power')
%         title('Real Gamma CV2map')
% end

NiceSave('DeltaCoupling',figfolder,baseName)

%%
% figure
% plot(PowerPhaseRatemap.phasebins, PowerPhaseRatemap.ratemap{excell}(end,:),'k')
% hold on
% plot(PowerPhaseRatemap.phasebins+2*pi, PowerPhaseRatemap.ratemap{excell}(end,:),'k')
% 
% plot(PowerPhaseRatemap_sim.phasebins, PowerPhaseRatemap_sim.meanrate(end,:),'r')
% hold on
% plot(PowerPhaseRatemap_sim.phasebins+2*pi, PowerPhaseRatemap_sim.meanrate(end,:),'r')
% % plot(PowerPhaseRatemap.phasebins, PowerPhaseRatemap.meanrate(end,:))
% % hold on
% % plot(PowerPhaseRatemap.phasebins+2*pi, PowerPhaseRatemap.meanrate(end,:))
% % plot(PowerPhaseRatemap.phasebins, PowerPhaseRatemap.meanrate(end./2,:))
% % hold on
% % plot(PowerPhaseRatemap.phasebins+2*pi, PowerPhaseRatemap.meanrate(end/2,:))
%% Figure
viewwin = bz_RandomWindowInIntervals(SleepState.ints.NREMstate,5);
figure
subplot(5,1,3)
    hold on
    for cc = 1:spikes.numcells
        whichcell = ISIStats.sorts.NREMstate.rate(cc);
        plot(simspikes.times{whichcell},cc.*ones(size(simspikes.times{whichcell})),'k.')
    end
    xlim(viewwin);ylim([0 spikes.numcells])
    ylabel('Simulated Spikes');
    box off
    set(gca,'xticklabels',[])
    
% subplot(5,1,4)
%     hold on
%     for cc = 1:spikes.numcells
%         whichcell = ISIStats.sorts.NREMstate.rate(cc);
%         plot(simspikes_poiss.times{whichcell},cc.*ones(size(simspikes_poiss.times{whichcell})),'k.')
%     end
%     xlim(viewwin);ylim([0 spikes.numcells])
%     ylabel('Poisson Spikes');
%     box off
%     set(gca,'xticklabels',[])
%   
subplot(5,1,2)
    plot(lfp.timestamps,lfp.data,'k')
    hold on
    plot(deLFP.timestamps,deLFP.data,'b')
    xlim(viewwin)
    ylabel({'LFP','Real Spikes'});
    box off
    set(gca,'xticklabels',[]);set(gca,'yticklabels',[])
    
subplot(5,1,1)
    hold on
    for cc = 1:spikes.numcells
        whichcell = ISIStats.sorts.NREMstate.rate(cc);
        plot(spikes.times{whichcell},cc.*ones(size(spikes.times{whichcell})),'k.')
    end
    xlim(viewwin);ylim([0 spikes.numcells])
    ylabel('Observed Spikes');
    box off
    set(gca,'xticklabels',[])
NiceSave('SimPopDelta',figfolder,baseName)
%% Figure
viewwin = bz_RandomWindowInIntervals(SleepState.ints.NREMstate,10);

figure
subplot(4,1,1)
    plot(lfp.timestamps,lfp.data,'k')
    hold on
    plot(deLFP.timestamps,deLFP.data,'b')
    plot(spikes.times{excell},max(double(lfp.data)).*ones(size(spikes.times{excell})),'k.')
    xlim(viewwin)
    ylabel('LFP');
    box off
    set(gca,'xticklabels',[]);set(gca,'yticklabels',[])
    
subplot(4,1,2)
    plot(GLMmodelfit(excell).timestamps,GLMmodelfit(excell).predRate,'k')
    hold on
    plot(simspikes.times{excell},max(GLMmodelfit(excell).predRate).*ones(size(simspikes.times{excell})),'k.')
    xlim(viewwin)
    ylabel({'Predicted Rate','Simulated Spikes'});
    box off
    set(gca,'xticklabels',[])
 

subplot(4,3,7)
    imagesc(PowerPhaseRatemap.phasebins,PowerPhaseRatemap.powerbins,...
        PowerPhaseRatemap.ratemap{excell})
    hold on
    imagesc(PowerPhaseRatemap.phasebins+2*pi,PowerPhaseRatemap.powerbins,...
        PowerPhaseRatemap.ratemap{excell})
    xlim([-pi 3*pi])
    caxis([0 15])
    axis xy
    colorbar
    xlabel('Phase');ylabel('Power')
    title('Real Delta Ratemap')
    
subplot(4,3,8)
    imagesc(PowerPhaseRatemap_sim.phasebins,PowerPhaseRatemap_sim.powerbins,...
        PowerPhaseRatemap_sim.ratemap{excell})
    hold on
    imagesc(PowerPhaseRatemap_sim.phasebins+2*pi,PowerPhaseRatemap_sim.powerbins,...
        PowerPhaseRatemap_sim.ratemap{excell})
    xlim([-pi 3*pi])
    caxis([0 15])
    axis xy
    colorbar
    xlabel('Phase');ylabel('Power')
    title('Simulated Delta Ratemap')

subplot(4,3,11)
    plot(GLMmodelfit(excell).powerbins,GLMmodelfit(excell).Rpower,'k','linewidth',2)
    hold on
    text(-1.5,1,['R0 = ',num2str(round(GLMmodelfit(excell).R0*lfp.samplingRate,1)),' Hz'])
    axis tight
    box off
    xlabel('Power');ylabel('GLM: Phase Coupling')
    %title('GLM Kernel')
    
subplot(4,3,12)
    plot(GLMmodelfit(excell).powerbins,GLMmodelfit(excell).Rratepower  ,'k','linewidth',2)
    hold on
    %text(-1.5,1,['R0 = ',num2str(round(GLMmodelfit(excell).R0*lfp.samplingRate,1)),' Hz'])
    axis tight
    box off
    xlabel('Power');ylabel('GLM: Phase Coupling')
    %title('GLM Kernel')

subplot(4,3,9)
    plot(ISIStats.ISIhist.logbins,ISIStats.ISIhist.NREMstate.log(excell,:),'k','linewidth',2)
    hold on
    plot(ISIstats_sim.ISIhist.logbins,ISIstats_sim.ISIhist.NREMstate.log(excell,:),'k--','linewidth',2)
   % plot(ISIstats_poiss.ISIhist.logbins,ISIstats_poiss.ISIhist.NREMstate.log(excell,:),'k:','linewidth',2)

    axis tight
    box off
    xlabel('ISI (s)')
    LogScale('x',10)
    legend('location','northeast','Real','Sim.')
    
NiceSave('ExCellDelta',figfolder,baseName)

%% Fit the 
expfit = fittype('a.*k.^(p+p0)','dependent',{'A'},'independent',{'p'},...
    'coefficients',{'a','k','p0'});

PLfit = fittype('a.*max(p-p0,0).^k','dependent',{'A'},'independent',{'p'},...
    'coefficients',{'a','k','p0'});

f_phasecoupling = fit(GLMmodelfit(excell).powerbins',GLMmodelfit(excell).Rpower,PLfit)
f_powercoupling = fit(GLMmodelfit(excell).powerbins',GLMmodelfit(excell).Rratepower-max(GLMmodelfit(excell).Rratepower),PLfit)

figure
subplot(2,2,1)
plot(f_phasecoupling,GLMmodelfit(excell).powerbins,GLMmodelfit(excell).Rpower)

subplot(2,2,2)
plot(f_powercoupling,GLMmodelfit(excell).powerbins,GLMmodelfit(excell).Rratepower-max(GLMmodelfit(excell).Rratepower))

NiceSave('PLFitPowerKernels',figfolder,baseName)


%% Check fit between parametric and nonparametric

figure
plot(GLMmodel_nonparam(excell).powerbins,GLMmodel_nonparam(excell).Rpower)
hold on

end

