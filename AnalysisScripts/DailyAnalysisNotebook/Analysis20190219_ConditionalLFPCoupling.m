function [ ] = Analysis20190219(basePath,figfolder)
% Date 0219/2019
%
%Goal: Develop function for conditional LFP coupling (i.e. power and phase
%coupling as function of some condition (ISI, for example)
%
%Plots
%-
%-
%
%% Load Header
%Initiate Paths
reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/Cicero_09102014';
basePath = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onDesktop/AG_HPC/Achilles_10252013'
%basePath = pwd;
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
statecolors = {'k','b','r',[0.6 0.6 0.6]};

[celltypes,~,typeidx] = unique(CellClass.label);
cellcolor = {'k','r'};


%% Load the LFP if needed

lfpchan = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.THchanID;
downsamplefactor = 5;
lfp = bz_GetLFP(lfpchan,...
    'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);
%Noralize the LFP
%lfp.data = NormToInt(single(lfp.data),'modZ', SleepState.ints.NREMstate,lfp.samplingRate);

%% Restrict to state

state = states{1};
ints = SleepState.ints.(state);

ISIStats.allspikes.instate = cellfun(@(X) InIntervals(X,ints),...
    ISIStats.allspikes.times,'UniformOutput',false);
%% Get complex valued wavelet transform at each timestamp
wavespec = bz_WaveSpec(lfp,'intervals',ints,'showprogress',true,'ncyc',12,...
    'nfreqs',100,'frange',[1 120]); 
%Mean-Normalize power within the interval
wavespec.data = bsxfun(@(X,Y) X./Y,wavespec.data,mean(abs(wavespec.data),1));

%% Get Power/Phase at each spike

%Get phase and normalized power at each spike time
ISIStats.allspikes.power = cellfun(@(X) interp1(wavespec.timestamps,abs(wavespec.data),X,...
    'nearest'),ISIStats.allspikes.times,'UniformOutput',false);
ISIStats.allspikes.phase = cellfun(@(X) interp1(wavespec.timestamps,angle(wavespec.data),X,...
    'nearest'),ISIStats.allspikes.times,'UniformOutput',false);

%%
ISIStats.spikes.X = cellfun(@(X) log10(X),ISIStats.allspikes.ISIs,'UniformOutput',false);
%Divide condition into bin (as in ConditionalHist)
Xbounds = [-2.5 1];
numXbins = 50;

Xedges = linspace(Xbounds(1),Xbounds(2),numXbins+1);
Xbins = Xedges(1:end-1)+ 0.5.*diff(Xedges([1 2]));
Xedges(1) = -inf;Xedges(end) = inf;

%First calculate the marginal probability of the conditional (X)
[Xhist,~,ISIStats.allspikes.XbinID] = cellfun(@(X) histcounts(X,Xedges),ISIStats.spikes.X,'UniformOutput',false);
%Xhist4norm = Xhist;Xhist4norm(Xhist4norm<=minX) = nan;

%%
minX = 100;
%%
%For each bin calculate - mean power of spikes, power-weighted MRL (in
%state)
%Mean Y given X
clear meanpower
clear mrl
clear mrlangle
for xx = 1:length(Xbins)
    meanpowertemp = cellfun(@(pow,binID,instate) nanmean(pow(binID==xx & instate,:)),...
        ISIStats.allspikes.power,ISIStats.allspikes.XbinID,ISIStats.allspikes.instate,'UniformOutput',false);
    
    %Mean resultant vector
    pMRVtemp = cellfun(@(pow,phs,binID,instate) nanmean(pow(binID==xx & instate,:).*exp(i.*phs(binID==xx & instate,:))),...
        ISIStats.allspikes.power,ISIStats.allspikes.phase,...
        ISIStats.allspikes.XbinID,ISIStats.allspikes.instate,'UniformOutput',false);
   
   for cc = 1:length(ISIStats.allspikes.times)
       meanpower(xx,:,cc)=meanpowertemp{cc};
       mrl(xx,:,cc)=abs(pMRVtemp{cc});
       mrlangle(xx,:,cc)=angle(pMRVtemp{cc});
       if Xhist{cc}(xx)<minX
           meanpower(xx,:,cc) = nan;
           mrl(xx,:,cc) = nan;
           mrlangle(xx,:,cc) = nan;
       end
           
   end
end

%%
for tt = 1:length(celltypes)
    conditionalpower.(celltypes{tt}) = nanmean(meanpower(:,:,CellClass.(celltypes{tt})),3);
    conditionalpMRL.(celltypes{tt}) = nanmean(mrl(:,:,CellClass.(celltypes{tt})),3);

end

%%
powermap = makeColorMap([0 0 0.8],[1 1 1],[0.8 0 0]);
figure
for tt = 1:length(celltypes)
subplot(3,3,tt*3)
colormap(gca,powermap)
    imagesc(Xbins,log2(wavespec.freqs), conditionalpower.(celltypes{tt})')
    colorbar
    caxis([0.4 1.6])
    LogScale('x',10);
    LogScale('y',2)
    axis xy
    xlabel('ISI (s)');ylabel('freq (Hz)')
    title((celltypes{tt}))
end   
    

for tt = 1:length(celltypes)
subplot(3,3,tt*3-1)
    imagesc(Xbins,log2(wavespec.freqs), conditionalpMRL.(celltypes{tt})')
    colorbar
    %caxis([0.5 1.5])
    LogScale('x',10);
    LogScale('y',2)
    axis xy
    xlabel('ISI (s)');ylabel('freq (Hz)')
    title((celltypes{tt}))
end   
  NiceSave(['ISIConditionalLFP',state],figfolder,baseName,'includeDate',true)
  
%%
figure
imagesc(Xbins,log2(wavespec.freqs),meanpower(:,:,20)')
colorbar
caxis([0.5 1.5])
LogScale('x',10);
LogScale('y',2)
axis xy
