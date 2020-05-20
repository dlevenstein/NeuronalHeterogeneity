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



    
%% Load the LFP

downsamplefactor = 2;
lfp = bz_GetLFP(SleepState.detectorinfo.detectionparms.SleepScoreMetrics.THchanID,...
    'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);

%%
dt = 0.01;
winsize = 12;
frange = [1 312];
nfreqs = 150;
[specslope,spec] = bz_PowerSpectrumSlope(lfp,winsize,dt,'spectype','wavelet',...
    'nfreqs',nfreqs,'showfig',true,'showprogress',true,'frange',frange,...
    'ints',SleepState.ints.WAKEstate);
    

%% Calculate MI Map (cell x frequency)

clear fPower
fPower.timestamps = specslope.timestamps;
%fPower(length(specslope.freqs)+1).data = [];
for ff = 1:length(specslope.freqs)
    bz_Counter(ff,length(specslope.freqs),'Freq')
    fPower.data = specslope.resid(:,ff);
        parfor cc = 1:spikes.numcells
        [FConditionalISIDist(ff,cc)] = ConditionalISI(spikes.times{cc},fPower,...
            'ints',SleepState.ints.WAKEstate,...
            'showfig',false);
        end
    %fPower(ff) = [];
end
%catch
%    continue
%end
%end
%% Recombine things

for cc = 1:spikes.numcells
    CellCondISIDist(cc) = bz_CollapseStruct(FConditionalISIDist(:,cc),3);
end
MInf =bz_CollapseStruct(CellCondISIDist,2);
MInf = squeeze(MInf.MutInf);

%% Mean PXY
for ff = 1:length(specslope.freqs)
    FreqCondISIDist(ff) = bz_CollapseStruct(FConditionalISIDist(ff,CellClass.pE),3,'mean');
end

%%

%%
figure
subplot(2,2,1)
    imagesc(log10(specslope.freqs),[1 spikes.numcells],MInf(ISIStats.sorts.WAKEstate.ratebyclass,:))
    LogScale('x',10)
    xlabel('f (Hz)')
    
subplot(2,2,2)
    plot(log10(specslope.freqs),mean(MInf,1))
    LogScale('x',10)
    xlabel('f (Hz)')
NiceSave('MutInf',figfolder,baseName)

    %%
    
    exf = 100;
    exf_IDX = find(abs(specslope.freqs-exf) == min(abs(specslope.freqs-exf)));
    
    figure
    imagesc(FreqCondISIDist(exf_IDX).pYX')
    %%
    excell = 3;
    figure
    imagesc(FConditionalISIDist(exf_IDX,excell).pYX')
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
