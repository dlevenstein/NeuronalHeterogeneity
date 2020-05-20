function [TH_ISIstats,ISIbythetaphase,ISIbytheta] = ThetaISIAnalysis(basePath,figfolder)
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


%% Load the LFP if needed

% lfpchan = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID;
% downsamplefactor = 1;
% lfp = bz_GetLFP(lfpchan,...
%     'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);
% %Noralize the LFP
% lfp.data = NormToInt(single(lfp.data),'modZ', SleepState.ints.NREMstate,lfp.samplingRate);


%% Load Gamma Modes
GammaFit = bz_LoadCellinfo(basePath,'GammaFit');
%% Normalize the brain state metrics
ThetaPower.timestamps = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.t_clus;
ThetaPower.data = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.thratio;


%%
excell = 3;
excellUID = spikes.UID(excell);
%Find the UID of the cell in the Gamma fit so match...
%Put the gamma fit parms to conditional dist in as initial parms
GFIDX = find(GammaFit.WAKEstate.cellstats.UID==excellUID);
cellGamma = GammaFit.WAKEstate.singlecell(GFIDX);
[ConditionalISIDist(cc)] = ConditionalISI(spikes.times{excell},ThetaPower,...
    'ints',SleepState.ints.WAKEstate,'GammaFitParms',cellGamma,...
    'basePath',basePath,'figname',['ThetaUID',num2str(excellUID)],...
    'figfolder',figfolder);
%%
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
[ConditionalISIDist(cc),ConditionalISIModes(cc)] = ConditionalISI(spikes.times{cc},ThetaPower,...
    'ints',SleepState.ints.WAKEstate,'GammaFitParms',cellGamma,...
    'showfig',false);
catch
    continue
end
end

%%
AllConditionalISIModes = bz_CollapseStruct(ConditionalISIModes);
%%
figure
subplot(2,2,1)
    hist(AllConditionalISIModes.GSCorr)
subplot(2,2,2)
    plot(AllConditionalISIModes.ASlogrates(AllConditionalISIModes.AScorr_p<0.05),...
        AllConditionalISIModes.ASCorr(AllConditionalISIModes.AScorr_p<0.05),'o')
    hold on
    plot(AllConditionalISIModes.ASlogrates(AllConditionalISIModes.AScorr_p>=0.05),...
        AllConditionalISIModes.ASCorr(AllConditionalISIModes.AScorr_p>=0.05),'.')
    
subplot(2,2,3)
    plot(AllConditionalISIModes.GSCorr,...
        log10(AllConditionalISIModes.GScorr_p),'o')
    
%%
downsamplefactor = 2;
lfp = bz_GetLFP(SleepState.detectorinfo.detectionparms.SleepScoreMetrics.THchanID,...
    'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);

%%
dt = 0.01;
winsize = 12;
frange = [1 312];
nfreqs = 200;
[specslope] = bz_PowerSpectrumSlope(lfp,winsize,dt,'spectype','wavelet',...
    'nfreqs',nfreqs,'showfig',true,'ints',SleepState.ints.WAKEstate,...
    'showprogress',true,'frange',frange);


%%
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
    [FConditionalISIDist(ff),FConditionalISIModes(ff)] = ConditionalISI(spikes.times{cc},fPower,...
        'ints',SleepState.ints.WAKEstate,'GammaFitParms',cellGamma,...
        'showfig',false);
    %fPower(ff) = [];
end
%catch
%    continue
%end
%end
%%
AllFConditionalISIModes = bz_CollapseStruct(FConditionalISIModes);
AllFConditionalISIDist = bz_CollapseStruct(FConditionalISIDist);

%%
figure
subplot(2,2,1)
plot(log10(specslope.freqs),AllFConditionalISIModes.GSCorr)
hold on
plot(xlim(gca),[0 0],'k--')
LogScale('x',10)

subplot(2,2,2)
plot(log10(specslope.freqs),AllFConditionalISIDist.MutInf)
hold on
plot(xlim(gca),[0 0],'k--')
LogScale('x',10)

NiceSave(['PowerModUID',num2str(excellUID)],figfolder,baseName)
end
