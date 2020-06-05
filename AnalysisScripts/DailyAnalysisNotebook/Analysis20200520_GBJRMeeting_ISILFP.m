function [PopConditional,PopConditional_phase] = ISIModulationAnalysis(basePath,figfolder)
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
    lfpchannel = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.THchanID; 

%% Load the LFP
% downsamplefactor = 2;
% lfp = bz_GetLFP(lfpchannel,...
%     'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);

%%

% ss = 1;
% state = states{ss};
% %Take only subset of time (random intervals) so wavelets doesn't break
% %computer (total 625s)
% usetime = 7200;%2500
% winsize = 25;
% if sum(diff(SleepState.ints.(state),1,2))>usetime
%     nwin = round(usetime./winsize);
%     %winsize = 30; %s
%     try
%         windows = bz_RandomWindowInIntervals( SleepState.ints.(state),winsize,nwin );
%     catch
%         windows = SleepState.ints.(state);
%     end
% else
%     windows = SleepState.ints.(state);
% end
        
%% Calculate Residuals
dt = 0.01;
winsize = 10;
frange = [1 312];
nfreqs = 150;
[specslope] = bz_PowerSpectrumSlope([],winsize,dt,'spectype','wavelet',...
    'nfreqs',nfreqs,'showfig',true,'showprogress',true,'frange',frange,...
    'saveMat',basePath,'saveName',['Chan',num2str(lfpchannel)],...
    'saveFolder','WavPSS');
    

%%
for tt = 1:length(celltypes)
    popspikes.(celltypes{tt}) = cat(1,spikes.times{CellClass.(celltypes{tt})});
end

%%
clear fPower


fPower.timestamps = specslope.timestamps;
%fPower(length(specslope.freqs)+1).data = [];
clear PopConditionalISIDist PopConditionalISIDist_phase PopConditionalISIDist_power
for ff = 1:length(specslope.freqs)
	bz_Counter(ff,length(specslope.freqs),'Freq')
    fPower.data = specslope.resid(:,ff);
    for tt = 1:length(celltypes)
        for ss = 1:2
        [PopConditionalISIDist.(states{ss}).(celltypes{tt})(ff)] = ConditionalISI(popspikes.(celltypes{tt}),fPower,...
            'ints',SleepState.ints.(states{ss}),...
            'showfig',false,'ISIDist',true);
        end
    end
    
    fPower.data = specslope.specgram(:,ff);
    for tt = 1
        for ss = 1:2
        [PopConditionalISIDist_power.(states{ss}).(celltypes{tt})(ff)] = ConditionalISI(popspikes.(celltypes{tt}),fPower,...
            'ints',SleepState.ints.(states{ss}),...
            'showfig',false,'ISIDist',false);
        end
    end
    
    fPower.data = specslope.phase(:,ff);
    for tt = 1:length(celltypes)
        for ss = 1:2
        [PopConditionalISIDist_phase.(states{ss}).(celltypes{tt})(ff)] = ConditionalISI(popspikes.(celltypes{tt}),fPower,...
            'ints',SleepState.ints.(states{ss}),...
            'showfig',false,'normtype','none','Xwin',[-pi pi],'ISIDist',true);
        end
    end
end



%%
PopConditional = bz_CollapseStruct(PopConditionalISIDist,3,'justcat',true);
PopConditional_phase = bz_CollapseStruct(PopConditionalISIDist_phase,3,'justcat',true);
PopConditionalISIDist_power = bz_CollapseStruct(PopConditionalISIDist_power,3,'justcat',true);

%%




figure
for tt = 1:length(celltypes)
subplot(3,3,1+(tt-1)*3)
hold on
for ss = 1:2
plot(log10(specslope.freqs),squeeze(PopConditional.(states{ss}).(celltypes{tt}).MutInf),'color',statecolors{ss},'linewidth',1)
end
LogScale('x',10)
xlabel('f (Hz)')
ylabel({(celltypes{tt}),'MI[ISI;Power]'})
title('ISI-Power Modulation')


subplot(3,3,2+(tt-1)*3)
hold on
for ss = 1:2
plot(log10(specslope.freqs),squeeze(PopConditional_phase.(states{ss}).(celltypes{tt}).MutInf),'color',statecolors{ss},'linewidth',1)
end
LogScale('x',10)
xlabel('f (Hz)')
ylabel('MI[ISI;Phase]')
title('ISI-Phase Modulation')


if tt == 1
    subplot(3,3,3)
    hold on
    for ss = 1
    plot(log10(specslope.freqs),squeeze(PopConditionalISIDist_power.(states{ss}).(celltypes{tt}).MutInf),'color',statecolors{ss},'linewidth',1)
    end
    LogScale('x',10)
    xlabel('f (Hz)')
    ylabel('MI[ISI;Phase]')
    title('ISI-Power (inc 1/f)')
    end
end

tt = 1;
ss = 1;
exf = 7;
exf_IDX = find(abs(specslope.freqs-exf) == min(abs(specslope.freqs-exf)));

phasex = linspace(-pi,3*pi,100);

subplot(3,2,5)
imagesc(PopConditional.(states{ss}).(celltypes{1}).Dist.Xbins(1,:,1),...
    PopConditional.(states{ss}).(celltypes{1}).Dist.Ybins(1,:,1),...
    PopConditional.(states{ss}).(celltypes{1}).Dist.pYX(:,:,exf_IDX)')
LogScale('y',10,'exp',true)
xlabel(['Power (',num2str(exf),'Hz)'])
ylabel('ISI (s)')
colorbar

subplot(3,2,6)
imagesc(PopConditional_phase.(states{ss}).(celltypes{1}).Dist.Xbins(1,:,1),...
    PopConditional_phase.(states{ss}).(celltypes{1}).Dist.Ybins(1,:,1),...
    PopConditional_phase.(states{ss}).(celltypes{1}).Dist.pYX(:,:,exf_IDX)')
hold on
imagesc(2*pi+PopConditional_phase.(states{ss}).(celltypes{1}).Dist.Xbins(1,:,1),...
    PopConditional_phase.(states{ss}).(celltypes{1}).Dist.Ybins(1,:,1),...
    PopConditional_phase.(states{ss}).(celltypes{1}).Dist.pYX(:,:,exf_IDX)')

plot(phasex,-cos(phasex),'k')
xlim([-pi 3*pi])
plot(xlim(gca),log10(1/exf).*[1 1],'k--')
ylabel('ISI (s)')
LogScale('y',10,'exp',true)
xlabel(['Phase (',num2str(exf),'Hz)'])
colorbar
bz_piTickLabel('x')

NiceSave('PopISIMod',figfolder,baseName)

%% Zoom on Theta (load?)
    %Load Gamma Modes
    GammaFit = bz_LoadCellinfo(basePath,'GammaFit');
    %Normalize the brain state metrics
    ThetaPower.timestamps = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.t_clus;
    ThetaPower.data = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.thratio;


%%
excell = 3;
excellUID = spikes.UID(excell);
%Find the UID of the cell in the Gamma fit so match...
%Put the gamma fit parms to conditional dist in as initial parms
GFIDX = find(GammaFit.WAKEstate.cellstats.UID==excellUID);
cellGamma = GammaFit.WAKEstate.singlecell(GFIDX);
[ThetaConditionalISI(excell)] = ConditionalISI(spikes.times{excell},ThetaPower,...
    'ints',SleepState.ints.WAKEstate,'GammaFitParms',cellGamma,...
    'basePath',basePath,'figname',['ThetaUID',num2str(excellUID)],...
    'figfolder',figfolder,'GammaFit',true);

%% Theta ISI modulation - all cells

parfor cc = 1:spikes.numcells
    cc
excellUID = spikes.UID(cc);
%Find the UID of the cell in the Gamma fit so match...
%Put the gamma fit parms to conditional dist in as initial parms
GFIDX = find(GammaFit.WAKEstate.cellstats.UID==excellUID);
if isempty(GFIDX)
    continue
end
cellGamma = GammaFit.WAKEstate.singlecell(GFIDX);
try
[ThetaConditionalISI(cc)] = ConditionalISI(spikes.times{cc},ThetaPower,...
    'ints',SleepState.ints.WAKEstate,'GammaFitParms',cellGamma,...
    'showfig',false,'GammaFit',true);
catch
    continue
end
end

%%
AllThetaISIModes = bz_CollapseStruct(ThetaConditionalISI);
AllThetaISIModes = bz_CollapseStruct(AllThetaISIModes.GammaModes,2);
%% Figure: Theta ISI modulation all cells
figure
subplot(2,2,1)
    hist(AllThetaISIModes.GSCorr)
subplot(2,2,2)
    plot(AllThetaISIModes.ASlogrates(AllThetaISIModes.AScorr_p<0.05),...
        AllThetaISIModes.AS_R(AllThetaISIModes.AScorr_p<0.05),'k.')
    hold on
    plot(AllThetaISIModes.ASlogrates(AllThetaISIModes.AScorr_p>=0.05),...
        AllThetaISIModes.AS_R(AllThetaISIModes.AScorr_p>=0.05),'.','color',[0.5 0.5 0.5])

    
    plot(AllThetaISIModes.GSlogrates(AllThetaISIModes.GScorr_p<0.05),...
        AllThetaISIModes.GS_R(AllThetaISIModes.GScorr_p<0.05),'.')
    axis tight
    box off
        plot(xlim(gca),[0 0],'k--')
        LogScale('x',10)
        xlabel('Mode Rate (Hz)')
        ylabel('Weight-Power Corr')
    
    
subplot(2,2,3)
    plot(AllThetaISIModes.GS_R,...
        log10(AllThetaISIModes.GScorr_p),'o')
    
    
 NiceSave('ThetaMod',figfolder,baseName)

    
%% Ex cell: all frequencies
%for cc = 1:spikes.numcells
    cc = excell;
excellUID = spikes.UID(cc);
%Find the UID of the cell in the Gamma fit so match...
%Put the gamma fit parms to conditional dist in as initial parms
GFIDX = find(GammaFit.WAKEstate.cellstats.UID==excellUID);
% if isempty(GFIDX)
%     continue
% end
cellGamma = GammaFit.WAKEstate.singlecell(GFIDX);
%try
clear fPower
fPower.timestamps = specslope.timestamps;
%fPower(length(specslope.freqs)+1).data = [];
for ff = 1:length(specslope.freqs)
    bz_Counter(ff,length(specslope.freqs),'Freq')
    fPower.data = specslope.resid(:,ff);
    [FConditionalISIDist(ff)] = ConditionalISI(spikes.times{cc},fPower,...
        'ints',SleepState.ints.WAKEstate,'GammaFitParms',cellGamma,...
        'showfig',false,'GammaFit',true);
    %fPower(ff) = [];
end
%catch
%    continue
%end
%end
%%
AllFConditionalISIDist = bz_CollapseStruct(FConditionalISIDist);
AllFConditionalISIModes = bz_CollapseStruct(AllFConditionalISIDist.GammaModes,2);
%%
%AllFConditionalISIDist = bz_CollapseStruct(AllFConditionalISIDist,2);
%%
figure
subplot(2,2,1)
plot(log10(specslope.freqs),-AllFConditionalISIModes.GS_R,'k')
hold on
%plot(log10(specslope.freqs(AllFConditionalISIModes.GScorr_p<0.05)),AllFConditionalISIModes.GS_R(AllFConditionalISIModes.GScorr_p<0.05),'o')
hold on
plot(xlim(gca),[0 0],'k--')
LogScale('x',10)
ylabel('AR-Power Corr')
xlabel('f (Hz)')
box off 

subplot(2,2,2)
plot(log10(specslope.freqs),AllFConditionalISIDist.MutInf,'k')
hold on
LogScale('x',10)
ylabel('MutInf')
xlabel('f (Hz)')
box off


NiceSave(['PowerModUID',num2str(excellUID)],figfolder,baseName)
end
