function [TH_ISIstats,ISIbythetaphase,ISIbytheta,ThetaISImodes] = ThetaISIAnalysis(basePath,figfolder)
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
%reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/';
% basePath = '/Users/dl2820/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/Cicero_09102014';
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/YMV12_171211';
% %basePath = pwd;
% %basePath = fullfile(reporoot,'Datasets/onProbox/AG_HPC/Achilles_11012013');
%figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%Load Stuff
sessionInfo = bz_getSessionInfo(basePath);
spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
states = fieldnames(SleepState.ints);
states{4} = 'ALL';
SleepState.ints.ALL = [0 Inf];
statecolors = {[0 0 0],[0 0 1],[1 0 0],[0.6 0.6 0.6]};

try
celltypes = CellClass.celltypes;
catch
    celltypes = unique(CellClass.label);
end
cellcolor = {'k','r'};


%% Load the LFP if needed
LFPMapFolder = [reporoot,'AnalysisScripts/AnalysisFigs/ISILFPMap'];

%Check for an LFP Map
try
    LFPMapFolder
    [ISILFPMap] = GetMatResults(LFPMapFolder,'ISILFPMap','baseNames',baseName);
    try
    region = 'vCTX';
    lfpchannel = ISILFPMap.MIMap.selectedchans.(region).channel;
    catch
            region = 'CA1';
            lfpchannel = ISILFPMap.MIMap.selectedchans.(region).channel;
    end
catch
    display('No ISILFP CHannel')
    lfpchannel = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.THchanID;
end

%%
downsamplefactor = 1;
lfp = bz_GetLFP(lfpchannel,...
    'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);
%Noralize the LFP
%lfp.data = NormToInt(single(lfp.data),'modZ', SleepState.ints.NREMstate,lfp.samplingRate);


%%
dt = 0.01;
winsize = 10;
frange = [1 312];
nfreqs = 150;
[specslope] = bz_PowerSpectrumSlope(lfp,winsize,dt,'spectype','wavelet',...
    'nfreqs',nfreqs,'showfig',true,'showprogress',true,'frange',frange);%,...
    %'saveMat',basePath,'saveName',['Chan',num2str(lfpchannel)],...
    %'saveFolder','WavPSS');


%%

for ss = 1:3
    specslope.instate.(states{ss}) = InIntervals(specslope.timestamps,SleepState.ints.(states{ss}));
    specstats.(states{ss}).meanresid = mean(specslope.resid(specslope.instate.(states{ss}),:),1);
    %specstats.(states{ss}).meanIRASA = mean(spec.IRASAsmooth(specslope.instate.(states{ss}),:),1);
    %specstats.(states{ss}).meanosci = mean(spec.osci(specslope.instate.(states{ss}),:),1);

    specstats.(states{ss}).meanspec = mean(specslope.specgram(specslope.instate.(states{ss}),:),1);
    
    
end
%%
figure
subplot(2,2,1)
hold on
for ss = 1:2
    plot(log2(specslope.freqs),log10(specstats.(states{ss}).meanspec),'color',statecolors{ss})
    %plot(log2(specslope.freqs),log10(specstats.(states{ss}).meanIRASA),'--','color',statecolors{ss})
end
LogScale('x',2)

subplot(2,2,2)
hold on
for ss = 1:2
    %plot(log2(specslope.freqs),specstats.(states{ss}).meanosci,'--','color',statecolors{ss})
    plot(log2(specslope.freqs),specstats.(states{ss}).meanresid,'color',statecolors{ss})
    
end
LogScale('x',2)


 NiceSave('Spec',figfolder,baseName)

%%
f_theta = [5.5 12];
thfreqs = (specslope.freqs>=f_theta(1) & specslope.freqs<=f_theta(2));
thratio = max((specslope.resid(:,thfreqs)),[],2);


%% Normalize the brain state metrics
ThetaPower.timestamps = specslope.timestamps;
ThetaPower.data = thratio;


%% Normalize the brain state metrics
%ThetaPower.timestamps = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.t_clus;
%ThetaPower.data = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.thratio;

for ss = 1:length(states)
   ThetaPower.instatetime.(states{ss}) = InIntervals(ThetaPower.timestamps,SleepState.ints.(states{ss}));
end
%% Divide into high and low theta time
ThetaPower.hilo_percentile = prctile(ThetaPower.data(ThetaPower.instatetime.WAKEstate),[35 65]);
ThetaIDX.threshs.hitheta = ThetaPower.hilo_percentile(2);
ThetaIDX.threshs.lotheta = ThetaPower.hilo_percentile(1);
ThetaIDX.timestamps = ThetaPower.timestamps;
ThetaIDX.states = zeros(size(ThetaIDX.timestamps));
ThetaIDX.states(ThetaPower.data>=ThetaIDX.threshs.hitheta & ThetaPower.instatetime.WAKEstate) = 1;
ThetaIDX.states(ThetaPower.data<=ThetaIDX.threshs.lotheta & ThetaPower.instatetime.WAKEstate) = 2;

ThetaIDX.statenames = {'hiTheta','loTheta'};
ThetaIDX.timestamps = ThetaPower.timestamps;

ThetaINT = bz_IDXtoINT(ThetaIDX);
%%
TH_ISIstats = bz_ISIStats(spikes,'ints',ThetaINT,'showfig',true,'cellclass',CellClass.label);
TH_ISIstats = rmfield(TH_ISIstats,'allspikes');
TH_ISIstats.cellinfo.celltype = CellClass;
%% Get the Theta LFP for phase

% lfpchan = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.THchanID;
% downsamplefactor = 5;
% lfp = bz_GetLFP(lfpchan,...
%     'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);
thetalfp = bz_Filter(lfp,'passband',f_theta);

ThetaPhase.data = thetalfp.phase;
ThetaPhase.timestamps = thetalfp.timestamps;

%% Calculate the conditional distributions

[ISIbytheta] = bz_ConditionalISI(spikes.times,ThetaPower,...
    'ints',SleepState.ints.WAKEstate,...
    'showfig',false,'GammaFit',false,'numXbins',20,'numISIbins',100);

[ISIbythetaphase] = bz_ConditionalISI(spikes.times,ThetaPhase,...
    'ints',ThetaINT.hiThetastate,...
    'showfig',false,'GammaFit',false,...
    'normtype','none','Xwin',[-pi pi],'numXbins',20,'numISIbins',100);

%% Population averages
for tt = 1:length(celltypes)
    ISIbytheta.pop.(celltypes{tt}) = nanmean(ISIbytheta.Dist.pYX(:,:,CellClass.(celltypes{tt})),3);
    ISIbythetaphase.pop.(celltypes{tt}) = nanmean(ISIbythetaphase.Dist.pYX(:,:,CellClass.(celltypes{tt})),3);
    ISIbythetaphase.popRate.(celltypes{tt}) = nanmean(ISIbythetaphase.Dist.SpikeRate(:,:,CellClass.(celltypes{tt})),3);
    [~,ISIbythetaphase.peakPhaseIdx.(celltypes{tt})] = max(ISIbythetaphase.popRate.(celltypes{tt}));
    ISIbythetaphase.peakPhase.(celltypes{tt}) = ISIbythetaphase.Dist.Xbins(1,ISIbythetaphase.peakPhaseIdx.(celltypes{tt}),1);
    
    ISIbythetaphase.celltypeidx.(celltypes{tt}) = CellClass.(celltypes{tt});
    ISIbytheta.celltypeidx.(celltypes{tt}) = CellClass.(celltypes{tt});
end
if length(celltypes)==1
    ISIbytheta.celltypeidx.pI = false(size(CellClass.pE));
end

%Shift to Peak first
ISIbythetaphase.shiftDist.pYX = circshift(ISIbythetaphase.Dist.pYX,...
    ISIbythetaphase.peakPhaseIdx.(celltypes{2}),1);
ISIbythetaphase.shiftDist.SpikeRate = circshift(ISIbythetaphase.Dist.SpikeRate,...
    ISIbythetaphase.peakPhaseIdx.(celltypes{2}),2);


%%
phasex = linspace(-pi,3*pi,100);

THlabels = {'hiThetastate','loThetastate'};
figure

for tt = 1:length(celltypes)
subplot(4,3,tt*3-2+6)
    imagesc(ISIbythetaphase.Dist.Xbins(1,:,1),ISIbythetaphase.Dist.Ybins(1,:,1),...
        ISIbythetaphase.pop.(celltypes{tt})')
    hold on
    imagesc(ISIbythetaphase.Dist.Xbins(1,:,1)+2*pi,ISIbythetaphase.Dist.Ybins(1,:,1), ...
        ISIbythetaphase.pop.(celltypes{tt})')
    
    %axis xy
    plot(phasex,-cos(phasex),'k')
    
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    %axis xy
    %xlim([-pi 3*pi])
    LogScale('y',10,'exp',true)
    ylabel('ISI (s)');xlabel('Theta Phase')
    %ylabel({(celltypes{tt}),'ISI (s)'});xlabel('Theta Ratio')
    %title((celltypes{tt}))
   % colorbar
%     if tt ==1 
%         caxis([0 0.02])
%     elseif tt==2
%          caxis([0 0.03])
%     end
    xlim([-pi 3*pi])
    plot(xlim(gca),log10(1/8).*[1 1],'w--')
    
    bounds = ylim(gca);
    yyaxis right
    plot(ISIbythetaphase.Dist.Xbins(1,:,1),log10(ISIbythetaphase.popRate.(celltypes{tt})),'r-','linewidth',1)
    plot(ISIbythetaphase.Dist.Xbins(1,:,1)+2*pi,log10(ISIbythetaphase.popRate.(celltypes{tt})),'r-','linewidth',1)
    ylim(-fliplr(bounds))
    LogScale('y',10,'nohalf',true)
    ylabel('Rate (Hz)')
end 

% subplot(3,3,7)
%     hold on
%     for tt = 1:length(celltypes)
%         plot(log10(TH_ISIstats.summstats.hiThetastate.meanrate(CellClass.(celltypes{tt}))),...
%             log10(TH_ISIstats.summstats.loThetastate.meanrate(CellClass.(celltypes{tt}))),...
%             '.','color',cellcolor{tt})
%     end
%     hold on
%     UnityLine
%     xlabel('hiTheta Rate');ylabel('loTheta Rate')
    

for tt = 1:length(celltypes) 
    
    subplot(4,4,tt+6)
        plot(TH_ISIstats.ISIhist.logbins,TH_ISIstats.meandists.hiThetastate.(celltypes{tt}).ISIdist,'k')
        hold on
        plot(TH_ISIstats.ISIhist.logbins,TH_ISIstats.meandists.loThetastate.(celltypes{tt}).ISIdist,'r')
        LogScale('x',10,'exp',true)
        title(celltypes{tt})
        box off 
        
    for ss = 1:length(THlabels)
        subplot(4,4,10+tt+(ss-1)*4)
            imagesc(TH_ISIstats.ISIhist.logbins,TH_ISIstats.ISIhist.logbins,...
            TH_ISIstats.meandists.(THlabels{ss}).(celltypes{tt}).Return)
            LogScale('xy',10,'exp',true)
            axis xy
    end
end


for tt = 1:length(celltypes)
subplot(4,3,tt*3-2)
    imagesc(ISIbytheta.Dist.Xbins(1,:,1),ISIbytheta.Dist.Ybins(1,:,1), ...
        ISIbytheta.pop.(celltypes{tt})')
    %hold on
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    %axis xy
    LogScale('y',10,'exp',true)
    ylabel('ISI (s)');xlabel('Theta Power')
    %ylabel({(celltypes{tt}),'ISI (s)'});xlabel('Theta Ratio')
    %title((celltypes{tt}))
   % colorbar
%     if tt ==1 
%         caxis([0 0.02])
%     elseif tt==2
%          caxis([0 0.03])
%     end
    xlim(ISIbytheta.Dist.Xbins(1,[1 end],1))
end 


    
    
% subplot(6,3,10)
%     hold on
%     for ss = [1]
%     plot(BShist.bins,BShist.(states{ss}).thratio,'color',statecolors{ss})
%     end
%     xlabel('Theta Ratio')
%     %xlim(ISIbytheta.Dist.Xbins(1,[1 end],1))

NiceSave('TH_ISIstats',figfolder,baseName)



%%
 GammaFit = bz_LoadCellinfo(basePath,'GammaFit');

 %% Example cell
 excell = 3;
 excellUID = spikes.UID(excell);
 GFIDX = find(GammaFit.WAKEstate.cellstats.UID==excellUID);
cellGamma = GammaFit.WAKEstate.singlecell(GFIDX);
 [ExThetaISIModes] = bz_ConditionalISI(spikes.times{excell},ThetaPower,...
    'ints',SleepState.ints.WAKEstate,'GammaFitParms',cellGamma,...
    'showfig',true,'GammaFit',true,'MutInf',false,'minX',0,...
    'figfolder',figfolder,'figname',['Theta',num2str(excellUID)],...
    'basePath',basePath,'numISIbins',100);
%% Theta ISI modulation - all cells

numAS = GammaFit.WAKEstate.parms.numAS;
GSModulation = nan(spikes.numcells,1);
ASModulation = nan(spikes.numcells,numAS);
ASlogRates = nan(spikes.numcells,numAS);
GSlogRates = nan(spikes.numcells,10);
ASweight = nan(spikes.numcells,numAS);
GSrate = nan(spikes.numcells,1);

for cc = 1:spikes.numcells
    bz_Counter(cc,spikes.numcells,'Cell')
    excellUID = spikes.UID(cc);
    %Find the UID of the cell in the Gamma fit so match...
    %Put the gamma fit parms to conditional dist in as initial parms
    GFIDX = find(GammaFit.WAKEstate.cellstats.UID==excellUID);
    cellGamma = GammaFit.WAKEstate.singlecell(GFIDX);
    if length(GFIDX)~=1

    else
        try
        [ThetaConditionalISI] = bz_ConditionalISI(spikes.times{cc},ThetaPower,...
            'ints',SleepState.ints.WAKEstate,'GammaFitParms',cellGamma,...
            'showfig',false,'GammaFit',true,'MutInf',false,'minX',0);
        
            GSModulation(cc) = ThetaConditionalISI.GammaModes.GS_R;
            ASModulation(cc,:) = ThetaConditionalISI.GammaModes.AS_R;
            ASlogRates(cc,:) = ThetaConditionalISI.GammaModes.ASlogrates;
            GSlogRates(cc,:) = ThetaConditionalISI.GammaModes.GSlogrates;
            ASweight(cc,:) = cellGamma.ASweights;
            GSrate(cc) = cellGamma.GSlogrates;
        catch
            cellGamma
            error(['Cell',num2str(cc)])
        end
    end
end

%%
ThetaISImodes.GSModulation = GSModulation;
ThetaISImodes.ASModulation = ASModulation;
ThetaISImodes.ASlogRates = ASlogRates;
ThetaISImodes.GSlogRates = GSlogRates;
ThetaISImodes.ASweight = ASweight;
ThetaISImodes.GSrate = GSrate;
%%
%AllThetaISIModes = bz_CollapseStruct(ThetaConditionalISI);
%AllThetaISIModes = bz_CollapseStruct(AllThetaISIModes.GammaModes,2);
%% Figure: Theta ISI modulation all cells
GScolor = [0.6 0.4 0];

figure
subplot(2,2,1)
    hist(ThetaISImodes.GSModulation)
subplot(2,2,2)
hold on
    for rr = 1:5
        scatter(ThetaISImodes.ASlogRates(:,rr),...
            ThetaISImodes.ASModulation(:,rr),20*ThetaISImodes.ASweight(:,rr)+0.00001,'k','filled')
    end
   % hold on

    
   plot(mean(ThetaISImodes.GSlogRates,2),...
       ThetaISImodes.GSModulation,'.','color',GScolor)
    axis tight
    box off
        plot(xlim(gca),[0 0],'k--')
        LogScale('x',10)
        xlabel('Mode Rate (Hz)')
        ylabel('Weight-Power Corr')
    
    
subplot(2,2,3)
    plot([1:10],log10(ThetaISImodes.GSlogRates),'.')
    
    
 NiceSave('ThetaMod',figfolder,baseName)
 

end
