function [PeriSWISIDist_next,PeriSWISIDist,SW_ISIstats,SWRConditionalISI_gamma ] = SharpWaveISIAnalysis(basePath,figfolder)
% Date XX/XX/20XX
%
%Question: 
%
%Plots
%-
%-
%
%% Load Header
%Initiate Paths
reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/';
% basePath = '/Users/dl2820/Dropbox/Research/Datasets/20140526_277um';
basePath = '/Users/dl2820/Dropbox/Research/Datasets/Cicero_09102014';
% %basePath = pwd;
% %basePath = fullfile(reporoot,'Datasets/onProbox/AG_HPC/Achilles_11012013');
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%Load Stuff
sessionInfo = bz_getSessionInfo(basePath);
spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
%ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
%states = fieldnames(SleepState.ints);
%states{4} = 'ALL';
%SleepState.ints.ALL = [0 Inf];
%statecolors = {'k','b','r',[0.6 0.6 0.6]};

try
    celltypes = CellClass.celltypes;
catch
    celltypes = unique(CellClass.label);
end
cellcolor = {'k','r'};

%%
SharpWaves = bz_LoadEvents(basePath,'SWR');
%FindbestRippleChannel
%bz_GetBestRippleChan

%%

%ripples = bz_FindRipples(basePath,rpchan)

%%
eventimes = SharpWaves.peaktimes;

%% ISI dist/return map: in and out of SWR
SWRints.SWR = SharpWaves.times + [-0.025 0.025];
SWRints.iSWR = [SharpWaves.times(1:end-1,2) SharpWaves.times(2:end,1)];
SWRints.iSWR = RestrictInts(SWRints.iSWR,SleepState.ints.NREMstate);
SWRints.SWR = RestrictInts(SWRints.SWR,SleepState.ints.NREMstate);
SW_ISIstats = bz_ISIStats(spikes,'ints',SWRints,'showfig',true,'cellclass',CellClass.label);
SW_ISIstats = rmfield(SW_ISIstats,'allspikes');

SW_ISIstats.cellinfo.celltype = CellClass;
%%
swrlabels = {'SWR','iSWR'};
%%

figure
subplot(2,2,1)
    hold on
    for tt = 1:length(celltypes)
        plot(log10(SW_ISIstats.summstats.SWR.meanrate(CellClass.(celltypes{tt}))),...
            log10(SW_ISIstats.summstats.iSWR.meanrate(CellClass.(celltypes{tt}))),...
            '.','color',cellcolor{tt})
    end
    hold on
    UnityLine
    xlabel('SWR Rate');ylabel('iSWR Rate')
    

for tt = 1:length(celltypes) 
    
    subplot(4,4,tt+6)
        plot(SW_ISIstats.ISIhist.logbins,SW_ISIstats.meandists.SWR.(celltypes{tt}).ISIdist,'k')
        hold on
        plot(SW_ISIstats.ISIhist.logbins,SW_ISIstats.meandists.iSWR.(celltypes{tt}).ISIdist,'r')
        LogScale('x',10,'exp',true)
        title(celltypes{tt})
        box off 
        
    for ss = 1:length(swrlabels)
        subplot(4,4,10+tt+(ss-1)*4)
            imagesc(SW_ISIstats.ISIhist.logbins,SW_ISIstats.ISIhist.logbins,...
            SW_ISIstats.meandists.(swrlabels{ss}).(celltypes{tt}).Return)
            LogScale('xy',10,'exp',true)
            axis xy
    end
end

NiceSave('SW_ISIStats',figfolder,baseName)


%%
[PeriSWISIDist_next] = bz_PeriEventISIDist(spikes.times,eventimes,...
    'numXbins',160,'minX',40,'whichISIs','next','winsize',[-0.5 0.5],...
    'cellclass','load','basePath',basePath);

[PeriSWISIDist] = bz_PeriEventISIDist(spikes.times,eventimes,...
    'numXbins',160,'minX',40,'whichISIs','both','winsize',[-0.5 0.5],...
    'cellclass','load','basePath',basePath);

%%
figure
for tt = 1:length(celltypes)
subplot(2,2,tt)
    imagesc(PeriSWISIDist.pop.(celltypes{tt}).Xbins,...
        PeriSWISIDist.pop.(celltypes{tt}).Ybins,PeriSWISIDist.pop.(celltypes{tt}).pYX')
    hold on
    plot(PeriSWISIDist.pop.(celltypes{tt}).Xbins,...
        log10(1./PeriSWISIDist.pop.(celltypes{tt}).rate),cellcolor{tt},'linewidth',1)
    axis tight
    plot([0 0],ylim(gca),'w--')
    LogScale('y',10,'nohalf',true)
    xlabel('t (s) - relative to SW');ylabel('ISI (s)')
    title(celltypes{tt})
    
    
subplot(2,2,tt+2)
    imagesc(PeriSWISIDist_next.pop.(celltypes{tt}).Xbins,...
        PeriSWISIDist_next.pop.(celltypes{tt}).Ybins,PeriSWISIDist_next.pop.(celltypes{tt}).pYX')
    hold on
    plot(PeriSWISIDist_next.pop.(celltypes{tt}).Xbins,...
        log10(1./PeriSWISIDist_next.pop.(celltypes{tt}).rate),cellcolor{tt},'linewidth',1)
    axis tight
    plot([0 0],ylim(gca),'w--')
    LogScale('y',10,'nohalf',true)
    xlabel('t (s) - relative to SW');ylabel('ISI (s)')
    title(celltypes{tt})
end
NiceSave('PeriSWISI',figfolder,baseName)


%% SWR/iSWR Gamma Fits

GammaFit = bz_LoadCellinfo(basePath,'GammaFit');
%%
SWRints.SWRstate = SWRints.SWR;
SWRints.iSWRstate = SWRints.iSWR;
clear SWRConditionalISI
for cc = 1:spikes.numcells
    bz_Counter(cc,spikes.numcells,'Cell')
    
	cellUID(cc) = spikes.UID(cc);
    GFIDX = find(GammaFit.NREMstate.cellstats.UID==cellUID(cc));
    if isempty(GFIDX)
        continue
    end
    
    cellGamma = GammaFit.NREMstate.singlecell(GFIDX);
    %stateIDX(cc).data = stateIDX(cc).states;
    %try
    [SWRConditionalISI(cc)] = bz_ConditionalISI(spikes.times(cc),SWRints,...
        'ints','input','normtype','none',...
        'GammaFitParms',cellGamma,'GammaFit',true,...
        'showfig',false);
    %catch

     %   continue
    %end
    
    %Number of spikes....
    %outfieldspikes(cc) = StateConditionalISI(cc).Dist.Xhist(2);
    
    SWRConditionalISI(cc).GammaModes.GSrate_all = GammaFit.NREMstate.sharedfit.GSlogrates(GFIDX);
    SWRConditionalISI(cc).GammaModes.GSweight_all  = GammaFit.NREMstate.sharedfit.GSweights(GFIDX);
    %SWRConditionalISI(cc).GammaModes.ASweight_all  = GammaFit.NREMstate.sharedfit.GSweights(GFIDX);


    %keyboard
end

%%
SWRConditionalISI_gamma = bz_CollapseStruct(SWRConditionalISI,3,'justcat');
SWRConditionalISI_gamma.GammaModes = bz_CollapseStruct(SWRConditionalISI_gamma.GammaModes,3,'justcat',true);

%%
figure
subplot(3,3,1)
plot(squeeze(1-SWRConditionalISI_gamma.GammaModes.GSweights(1,2,:)),squeeze(1-SWRConditionalISI_gamma.GammaModes.GSweights(1,1,:)),'k.')
hold on
xlim([0 1]);ylim([0 1])
UnityLine
xlabel('Activation Ratio: iSWR');ylabel('Activation Ratio: SWR')

subplot(3,3,2)
plot(squeeze(SWRConditionalISI_gamma.GammaModes.GSlogrates(1,2,:)),squeeze(SWRConditionalISI_gamma.GammaModes.GSlogrates(1,1,:)),'k.')
hold on
axis tight
UnityLine
ylabel('GS Rate: SWR');xlabel('GS Rate: iSWR')
LogScale('xy',10)
box off

diffASweight = (SWRConditionalISI_gamma.GammaModes.ASweights(1,:,:)-SWRConditionalISI_gamma.GammaModes.ASweights(2,:,:));
diffGS = (SWRConditionalISI_gamma.GammaModes.GSweights(1,1,:))-(SWRConditionalISI_gamma.GammaModes.GSweights(1,2,:));
subplot(3,3,3)
hold on
for aa = 1:5
scatter(-SWRConditionalISI_gamma.GammaModes.ASlogrates(1,aa,:),...
    log10(SWRConditionalISI_gamma.GammaModes.ASCVs(1,aa,:)),...
    60*mean(SWRConditionalISI_gamma.GammaModes.ASweights(:,aa,:),1)+eps,...
    squeeze(diffASweight(1,aa,:)),'filled')
end
scatter(-SWRConditionalISI_gamma.GammaModes.GSlogrates(1,1,:),...
    log10(SWRConditionalISI_gamma.GammaModes.GSCVs(1,1,:)),...
    10,...
    squeeze(diffGS))
axis tight
hold on
plot([-2.5 1],[0 0],'k--')
colorbar
LogScale('x',10,'exp',true,'nohalf',true)
LogScale('y',10,'exp',false,'nohalf',true)
caxis([-0.1 0.1])
crameri('vik','pivot',0)
xlabel('Mean'); ylabel('CV')


%diffAR = (1-cellISIStats.GammaModes.GSweights(1,3,:))-(1-cellISIStats.GammaModes.GSweights(1,2,:));

% subplot(2,2,4)
% plot(squeeze(cellISIStats.MIskaggs),squeeze(diffAR),'.')
NiceSave('GSAS_SWR',figfolder,baseName)

