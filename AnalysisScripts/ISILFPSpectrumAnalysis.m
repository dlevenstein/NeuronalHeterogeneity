function [PSSConditionalISIDist,MutInf,HiLowISIStats,...
    PSSConditionalGamma,GammaParms,PSScorr] = ISILFPSpectrumAnalysis(basePath,figfolder)
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
reporoot = '/gpfs/data/buzsakilab/DL/NeuronalHeterogeneity/';
%reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/';
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/Cicero_09102014';
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


%%

% LFPMapFolder = [reporoot,'AnalysisScripts/AnalysisFigs/ISILFPMap'];
% 
% %Check for an LFP Map
% try
%     [ISILFPMap] = GetMatResults(LFPMapFolder,'ISILFPMap','baseNames',baseName);
% catch
%     error('No Channel selected')
% end
% 
% region = 'fCTX';
% %If One exists: bz_tagChannel
% lfpchannel = ISILFPMap.MIMap.selectedchans.(region).channel;
% usecells = ISILFPMap.MIMap.(ISILFPMap.MIMap.selectedchans.(region).regname).UIDs;
% %If it doesnt, skip all the mode stuff and just run bz_ISILFPMap
% %Question: what to do about pir/bla
% %Note: need to only use UIDs from the right region

    %%
    LFPMapFolder = [reporoot,'AnalysisScripts/AnalysisFigs/ISILFPMap'];
    %Check for an LFP Map
    try
        [ISILFPMap] = GetMatResults(LFPMapFolder,'ISILFPMap','baseNames',baseName);
        try
            region = 'fCTX';
            lfpchannel = ISILFPMap.MIMap.selectedchans.(region).channel;
        catch
            try
                region = 'vCTX';
                lfpchannel = ISILFPMap.MIMap.selectedchans.(region).channel;
            catch
                region = 'THAL';
                lfpchannel = ISILFPMap.MIMap.selectedchans.(region).channel;
            end
        end
        %usecells = ISILFPMap.MIMap.(ISILFPMap.MIMap.selectedchans.(region).regname).UIDs;
    catch
        display('No Channel selected, using theta channel (CA1?)')
        lfpchannel = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.THchanID; 
    end


%If One exists: bz_tagChannel

    
    %Load from the analysisresults. and tag channel.
    %Then run [ISILFPMap] =
    %bz_ISILFPMap(basePath,varargin) to save again and save with channel
    %marking? (also 
    %In the future, here, tagChannel for each recording!
%% Load the LFP
% downsamplefactor = 1;
% lfp = bz_GetLFP(lfpchannel,...
%     'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);

%%
  
%% Calculate Residuals
PSSmethod = 'FFT';
switch PSSmethod
    case 'Wavelet'
        %Wavelet, load
        dt = 0.01;
        winsize = 10;
        frange = [1 312];
        nfreqs = 150;
        [specslope] = bz_PowerSpectrumSlope([],winsize,dt,'spectype','wavelet',...
            'nfreqs',nfreqs,'showfig',true,'showprogress',true,'frange',frange,...
            'saveMat',basePath,'saveName',['Chan',num2str(lfpchannel)],...
            'saveFolder','WavPSS');
    case 'FFT'
        %% FFT
        downsamplefactor = 1;
        lfp = bz_GetLFP(lfpchannel,...
            'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);
        dt = 0.25;
        winsize = 1;
        frange = [2 128];
        nfreqs = 150;
        [specslope] = bz_PowerSpectrumSlope(lfp,winsize,dt,'spectype','fft',...
            'nfreqs',nfreqs,'showfig',true,'showprogress',true,'frange',frange);
    
end
% PSS Options: FFT (1 or 2s window). Wavelet. Wavelet -> smooth (1-2s)
%smooth osci/power too?
%%
specslope = rmfield(specslope,'phase');
specslope = rmfield(specslope,'intercept');
specslope = rmfield(specslope,'rsq');

%%
for ss = 1:3

    %%


    [PSSConditionalISIDist.(states{ss})] = bz_ConditionalISI(spikes.times,specslope,...
        'ints',SleepState.ints.(states{ss}),...
        'showfig',false,'GammaFit',false,'minX',0);
    %% Here... gamma fit for PSS...
    %%
    MutInf.(states{ss}).PSS = nan(1,spikes.numcells);
    MutInf.(states{ss}).Power = nan(length(specslope.freqs),spikes.numcells);
    MutInf.(states{ss}).Osci = nan(length(specslope.freqs),spikes.numcells);

    MutInf.(states{ss}).PSS = squeeze(PSSConditionalISIDist.(states{ss}).MutInf);

    clear fPower
    fPower.timestamps = specslope.timestamps;
    for ff = 1:length(specslope.freqs)
        bz_Counter(ff,length(specslope.freqs),'Freq')
        %display(num2str(ff))
        %ff
        fPower(ff).timestamps = downsample(specslope.timestamps,3);
        fPower(ff).data = downsample(specslope.resid(:,ff),3);

            [cdist] = bz_ConditionalISI(spikes.times,fPower(ff),...
                'ints',SleepState.ints.(states{ss}),...
                'showfig',false,'GammaFit',false,'minX',0,'ISIDist',false);

            tempOsci(ff,:) = squeeze(cdist.MutInf);

        fPower(ff).data = downsample(specslope.specgram(:,ff),3);
            [cdist] = bz_ConditionalISI(spikes.times,fPower(ff),...
                'ints',SleepState.ints.(states{ss}),...
                'showfig',false,'GammaFit',false,'minX',0,'ISIDist',false);

            tempPower(ff,:) = squeeze(cdist.MutInf);

          fPower(ff).timestamps = 0;
          fPower(ff).data = 0;

    end
    MutInf.(states{ss}).Osci = tempOsci;
    MutInf.(states{ss}).Power = tempPower;

end

    %%
for ss = 1:3
    for tt = 1:length(celltypes)
        MeanMI.(states{ss}).(celltypes{tt}).Osci = nanmean((MutInf.(states{ss}).Osci(:,CellClass.(celltypes{tt}))),2);
        MeanMI.(states{ss}).(celltypes{tt}).Power = nanmean((MutInf.(states{ss}).Power(:,CellClass.(celltypes{tt}))),2);

        %Conditioned on spike AND power
        PSSConditionalISIDist.(states{ss}).(celltypes{tt}) = ...
            nanmean(PSSConditionalISIDist.(states{ss}).Dist.pYX(:,:,CellClass.(celltypes{tt})),3);
        
        %Just conditioned on power
%         PSSConditionalISIDist.(states{ss}).(celltypes{tt}) = ...
%             nanmean(PSSConditionalISIDist.(states{ss}).Dist.XYprob(:,:,CellClass.(celltypes{tt})),3);
        
        PSSConditionalISIDist.(states{ss}).meanlograte.(celltypes{tt}) = ...
            nanmean(log10(PSSConditionalISIDist.(states{ss}).Dist.SpikeRate(:,:,CellClass.(celltypes{tt}))),3);
    end
end
MutInf.CellClass = CellClass;
%%
tt=1;

figure
for ss = 1:3

subplot(3,5,1+5.*(ss-1))
BoxAndScatterPlot({(MutInf.(states{ss}).PSS(CellClass.(celltypes{tt})))})
if ss == 1
PSSlim = ylim(gca); PSSlim(1) = 0;
end
ylim(PSSlim)
box off
ylabel('MI - PSS')

subplot(3,5,[2:3]+5.*(ss-1))
plot(log2(specslope.freqs),MeanMI.(states{ss}).(celltypes{tt}).Osci,'k')
hold on 
plot(log2(specslope.freqs),MeanMI.(states{ss}).(celltypes{tt}).Power,'r:')
LogScale('x',2)
axis tight
ylim(PSSlim)
box off
ylabel('MI - Power')
xlabel('freq (Hz)')


subplot(3,5,[4:5]+5.*(ss-1))
imagesc(PSSConditionalISIDist.(states{ss}).Dist.Xbins(1,:,1),...
    PSSConditionalISIDist.(states{ss}).Dist.Ybins(1,:,1),...
    PSSConditionalISIDist.(states{ss}).(celltypes{tt})')
hold on
plot(PSSConditionalISIDist.(states{ss}).Dist.Xbins(1,:,1),...
    -PSSConditionalISIDist.(states{ss}).meanlograte.(celltypes{tt}),'r','linewidth',2)
LogScale('y',10,'exp',true,'nohalf',true)
xlabel('PSS (%ile)');ylabel('ISI (s)')
bz_AddRightRateAxis

end
NiceSave('ISIModSpectrum',figfolder,baseName)

%% Gamma Modes for PSS
    %Load Gamma Modes for PSS
GammaFit = bz_LoadCellinfo(basePath,'GammaFit');



for ss = 1:3
    
    GammaParms.(states{ss}).GSrate = nan(size(MutInf.WAKEstate.PSS));
    GammaParms.(states{ss}).GSCV = nan(size(MutInf.WAKEstate.PSS));
    GammaParms.(states{ss}).GSweight = nan(size(MutInf.WAKEstate.PSS));
    GammaParms.(states{ss}).GSmod = nan(size(MutInf.WAKEstate.PSS));
    GammaParms.(states{ss}).GSmod_p = nan(size(MutInf.WAKEstate.PSS));
    
for cc = 1:spikes.numcells
    bz_Counter(cc,spikes.numcells,'Cell')
    
	cellUID(cc) = spikes.UID(cc);
    GFIDX = find(GammaFit.(states{ss}).cellstats.UID==cellUID(cc));
    if isempty(GFIDX)
        continue
    end
    
    cellGamma = GammaFit.(states{ss}).singlecell(GFIDX);
    %stateIDX(cc).data = stateIDX(cc).states;
    try
    [PSSConditionalGamma.(states{ss})(cc)] = bz_ConditionalISI(spikes.times(cc),specslope,...
        'ints',SleepState.ints.(states{ss}),...
        'GammaFitParms',cellGamma,'GammaFit',true,...
        'showfig',false);
    catch

        continue
    end
    
    GammaParms.(states{ss}).GSrate(cc) = cellGamma.GSlogrates;
    GammaParms.(states{ss}).GSCV(cc) = cellGamma.GSCVs;
    GammaParms.(states{ss}).GSweight(cc) = cellGamma.GSweights;
    GammaParms.(states{ss}).GSmod(cc) = PSSConditionalGamma.(states{ss})(cc).GammaModes.GS_R;
    GammaParms.(states{ss}).GSmod_p(cc) = PSSConditionalGamma.(states{ss})(cc).GammaModes.GScorr_p;
    
end

end
%%
for ss = 1:3
    PSSConditionalGamma.modes.(states{ss}) = bz_CollapseStruct([PSSConditionalGamma.(states{ss}).GammaModes],3,'justcat');
    PSSConditionalGamma.dist.(states{ss}) = bz_CollapseStruct([PSSConditionalGamma.(states{ss}).Dist],3,'justcat');
end

%% PSS-conditional ISI dists groups
modgroups = {'unmodcells','negmodcells','posmodcells'};
for ss = 1:3
PSSConditionalGamma.dist.(states{ss}).posmodcells = ...
    PSSConditionalGamma.modes.(states{ss}).GS_R<0 & PSSConditionalGamma.modes.(states{ss}).GScorr_p<=0.05;
PSSConditionalGamma.dist.(states{ss}).negmodcells = ...
    PSSConditionalGamma.modes.(states{ss}).GS_R>0 & PSSConditionalGamma.modes.(states{ss}).GScorr_p<=0.05;
PSSConditionalGamma.dist.(states{ss}).unmodcells = ...
    PSSConditionalGamma.modes.(states{ss}).GScorr_p>0.05;

PSSConditionalISIDist.(states{ss}).posmodcells = ...
            nanmean(PSSConditionalGamma.dist.(states{ss}).pYX(:,:,PSSConditionalGamma.dist.(states{ss}).posmodcells),3);
PSSConditionalISIDist.(states{ss}).negmodcells = ...
            nanmean(PSSConditionalGamma.dist.(states{ss}).pYX(:,:,PSSConditionalGamma.dist.(states{ss}).negmodcells),3);
        PSSConditionalISIDist.(states{ss}).unmodcells = ...
            nanmean(PSSConditionalGamma.dist.(states{ss}).pYX(:,:,PSSConditionalGamma.dist.(states{ss}).unmodcells),3);
end
%%
figure
for ss = 1:3
subplot(3,5,1+(ss-1)*5)
BoxAndScatterPlot({squeeze(-PSSConditionalGamma.modes.(states{ss}).GS_R)})
hold on
plot(xlim(gca),[0 0],'k--')
ylabel('AS Modulation')
box off

subplot(3,5,[2:3]+5.*(ss-1))
hold on
keepcells = PSSConditionalGamma.modes.(states{ss}).GScorr_p<=0.05;
%keepcells = true(size(keepcells))
for aa = 1:5
    %keepmodes = keepcells&mean(PSSConditionalGamma.modes.(states{ss}).ASweights(:,aa,:),1)>0.02;
    keepmodes = (PSSConditionalGamma.modes.(states{ss}).AScorr_p(:,aa,:))<=0.05;
    %keepmodes = true(size(keepmodes))
scatter(-PSSConditionalGamma.modes.(states{ss}).ASlogrates(1,aa,keepmodes),...
    log10(PSSConditionalGamma.modes.(states{ss}).ASCVs(1,aa,keepmodes)),...
    60*mean(PSSConditionalGamma.modes.(states{ss}).ASweights(:,aa,keepmodes),1)+eps,...
    squeeze(PSSConditionalGamma.modes.(states{ss}).AS_R(1,aa,keepmodes)),'filled')
end
scatter(-PSSConditionalGamma.modes.(states{ss}).GSlogrates(1,1,keepcells),...
    log10(PSSConditionalGamma.modes.(states{ss}).GSCVs(1,1,keepcells)),...
    10,...
    squeeze(PSSConditionalGamma.modes.(states{ss}).GS_R(1,1,keepcells)))
colorbar
axis tight
caxis([-0.1 0.1])
crameri('vik','pivot',0)
xlabel('Mean');ylabel('CV')


%tt = 1
subplot(3,5,[4:5]+5.*(ss-1))
imagesc(PSSConditionalISIDist.(states{ss}).Dist.Xbins(1,:,1),...
    PSSConditionalISIDist.(states{ss}).Dist.Ybins(1,:,1),...
    PSSConditionalISIDist.(states{ss}).pE')
end
NiceSave('GSASModPSS',figfolder,baseName)

%%
figure
for ss = 1:3
    for mm = 1:3
 subplot(3,3,mm+3.*(ss-1))
imagesc(PSSConditionalISIDist.(states{ss}).Dist.Xbins(1,:,1),...
    PSSConditionalISIDist.(states{ss}).Dist.Ybins(1,:,1),...
    PSSConditionalISIDist.(states{ss}).(modgroups{mm})')  
    xlabel('PSS')
    if ss == 1
        title(modgroups{mm})
    end
    if mm == 1
        ylabel(states{ss})
    end
    
    LogScale('y',10,'exp',true,'nohalf',true)
    xlabel('PSS (%ile)');ylabel('ISI (s)')
    bz_AddRightRateAxis
    end
    

NiceSave('GSModGroups',figfolder,baseName)
end
%% Hi/Low PSS ints
PSSthresh = [0.75 0.25];
smoothwin = 1; %s
smoothdata = smooth(specslope.data,round(specslope.samplingRate.*smoothwin));
for ss = 1:3
    instate = InIntervals(specslope.timestamps,SleepState.ints.(states{ss}));
    percentilenorm = NormToInt(smoothdata(instate),'percentile');

    IDX.timestamps = specslope.timestamps(instate);
    IDX.states = zeros(size(IDX.timestamps));
    IDX.states(percentilenorm<=PSSthresh(1)) = 1;
    IDX.states(percentilenorm>=PSSthresh(2)) = 2;
    IDX.statenames = {'LowPSS','HighPSS'};
    INT = bz_IDXtoINT(IDX);

    tempstruct = bz_ISIStats(spikes,'ints',INT,'showfig',false);
     HiLowISIStats.(states{ss}) = tempstruct.ISIhist;
%Later - put this as temp and save the things you want (i.e. no
%allspikes...)
end

%%
hilow = {'LowPSSstate','HighPSSstate'};

for ss = 1:3
    %ISIstatshist.(states{ss}) = HiLowISIStats.(states{ss}).ISIhist;
    for hl = 1:2
        meanISIhist.(states{ss}).(hilow{hl}).return =...
            mean(HiLowISIStats.(states{ss}).(hilow{hl}).return(:,:,CellClass.pE),3);
        meanISIhist.(states{ss}).(hilow{hl}).logdist =...
            mean(HiLowISIStats.(states{ss}).(hilow{hl}).log(CellClass.pE,:),1);
    end
end
meanISIhist.logbins = HiLowISIStats.(states{ss}).logbins;
%% In/Out Field: REturn Maps
histcolors = flipud(gray);
NREMhistcolors = makeColorMap([1 1 1],[0 0 0.8]);
REMhistcolors = makeColorMap([1 1 1],[0.8 0 0]);
statecolormap = {histcolors,NREMhistcolors,REMhistcolors};

figure
for ss = 1:3

subplot(3,3,6+ss)
hold on
for hl = 1:2
    plot(meanISIhist.logbins,meanISIhist.(states{ss}).(hilow{hl}).logdist)
end
    axis tight
    box off
    LogScale('x',10,'exp',true,'nohalf',true)

%legend(cellISIStats.statenames{2:3},'location','southoutside')



for hl = 1:2
subplot(3,3,(hl-1)*3+ss)
    imagesc(meanISIhist.logbins,meanISIhist.logbins,meanISIhist.(states{ss}).(hilow{hl}).return)
    axis xy
    axis tight
    LogScale('xy',10,'exp',true,'nohalf',true)
    colormap(gca,statecolormap{ss})

end
%legend(cellISIStats.statenames{1:3})
end

NiceSave('HiLowPSSReturn',figfolder,baseName)



%% PSS-Rate correlation, all cells
dt = 0.25;
binsize = 1;
spikemat = bz_SpktToSpkmat(spikes.times,'dt',dt,'binsize',binsize);
%SpkToSpkmat (same windows as PSS)
%Correlation (each E cell by GS rate), I cells
%Compare correlation to AS modulation
spikemat.PSS = interp1(specslope.timestamps,specslope.data,spikemat.timestamps);



for ss = 1:3
    spikemat.(states{ss}) = InIntervals(spikemat.timestamps,SleepState.ints.(states{ss}));
    PSScorr.(states{ss}) = corr(spikemat.data(spikemat.(states{ss}),:),spikemat.PSS(spikemat.(states{ss})),'type','spearman');
end


%%
figure
for ss = 1:3
    GammaParms.(states{ss}).meanrate = ISIStats.summstats.(states{ss}).meanrate;
    GammaParms.(states{ss}).sigmod = GammaParms.(states{ss}).GSmod_p <=0.05;
    subplot(4,3,ss)
    hold on
    for tt = 1:2
        plot(log10(GammaParms.(states{ss}).meanrate(CellClass.(celltypes{tt}))),...
            PSScorr.(states{ss})(CellClass.(celltypes{tt})),'.','color',cellcolor{tt})
    end
    axis tight; box off
    plot(xlim(gca),[0 0],'k--')
    xlabel('Mean Rate (Hz)');ylabel('PSS-Rate Corr')
    
    subplot(4,3,ss+3)
        hold on
        ScatterWithLinFit(GammaParms.(states{ss}).GSrate,PSScorr.(states{ss}))
        axis tight; box off
        plot(xlim(gca),[0 0],'k--')
        xlabel('GS Rate (Hz)');ylabel('PSS-Rate Corr')
        
    subplot(4,3,ss+6)
        hold on
        plot(-GammaParms.(states{ss}).GSmod(GammaParms.(states{ss}).sigmod),...
            PSScorr.(states{ss})(GammaParms.(states{ss}).sigmod),'k.')
        plot(-GammaParms.(states{ss}).GSmod(~GammaParms.(states{ss}).sigmod),...
            PSScorr.(states{ss})(~GammaParms.(states{ss}).sigmod),'.','color',[0.5 0.5 0.5])
        axis tight; box off
        plot(xlim(gca),[0 0],'k--')
        xlabel('AS Mod(Hz)');ylabel('PSS-Rate Corr')
        
    subplot(4,3,ss+9)
        hold on
        plot(GammaParms.(states{ss}).GSrate(GammaParms.(states{ss}).sigmod),...
            -GammaParms.(states{ss}).GSmod(GammaParms.(states{ss}).sigmod),'k.')
        plot(GammaParms.(states{ss}).GSrate(~GammaParms.(states{ss}).sigmod),...
            -GammaParms.(states{ss}).GSmod(~GammaParms.(states{ss}).sigmod),'.','color',[0.5 0.5 0.5])
        axis tight; box off
        plot(xlim(gca),[0 0],'k--')
        ylabel('AS Mod(Hz)');ylabel('GS Rate')
        
end
NiceSave('PSSRateCorr',figfolder,baseName)

end
