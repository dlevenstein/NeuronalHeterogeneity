function [ISIbyHD_align,MutInfo ] = HeadDirectionTuningAnalysis(basePath,figfolder)
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
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = [reporoot,'/Datasets/onProbox/AG_HPC/Achilles_10252013'];
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/Mouse12-120807';
%basePath = [reporoot,'/Datasets/onProbox/AG_HPC/Achilles_10252013'];
%basePath = pwd;
%figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%Load Stuff
sessionInfo = bz_getSessionInfo(basePath);
spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
%ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
states = fieldnames(SleepState.ints);
states{4} = 'ALL';
SleepState.ints.ALL = [0 Inf];
statecolors = {'k','b','r',[0.6 0.6 0.6]};

%[celltypes,~,typeidx] = unique(CellClass.label);
%cellcolor = {'k','r'};


%%
headdir.samplingRate = 39.06; %Hz
hdfilename = fullfile(basePath,[baseName,'.ang']);
headdir.data = importdata(hdfilename);
headdir.data(headdir.data==-1)=nan;
headdir.timestamps = [1:length(headdir.data)]'./headdir.samplingRate; 

%%
% nantimes = isnan(headdir.position.lin);
% headdir.data = headdir.position.lin(~(nantimes));
%Possible here: remove drops for more than a... second?
%minjump = 2;
%[ headdir.data ] = NanPadJumps( headdir.timestamps(~(nantimes)),headdir.data,minjump );
%headdir.data = interp1(headdir.timestamps(~(nantimes)),headdir.data,headdir.timestamps);
[ISIbyHD] = bz_ConditionalISI(spikes.times,headdir,...
    'ints',SleepState.ints.WAKEstate,...
    'showfig',false,'GammaFit',false,'numXbins',25,'numISIbins',100,...
    'normtype','none','Xwin',[0 2.*pi]);

%%
[~,ISIbyHD.fieldpeak] = max(ISIbyHD.Dist.SpikeRate,[],2);
ISIbyHD.fieldpeak = ISIbyHD.Dist.Xbins(ISIbyHD.fieldpeak);

%% Compare ISI MI vs Rate MI
MutInfo.ISI = squeeze(ISIbyHD.MutInf);

binsizes = logspace(-2.5,1.5,25);
for bb = 1:length(binsizes)
    bz_Counter(bb,length(binsizes),'Bin')
spkmat = bz_SpktToSpkmat(spikes.times,'dt',binsizes(bb),'binsize',binsizes(bb),'units','rate');
spkmat.pos = interp1(headdir.timestamps,headdir.data,spkmat.timestamps);
spkmat.InWake = InIntervals(spkmat.timestamps,SleepState.ints.WAKEstate);

for cc = 1:spikes.numcells
    MutInfo.Rate_BinCompare(cc,bb) = mutualinfo(spkmat.data(spkmat.InWake,cc),spkmat.pos(spkmat.InWake));
end
end

%%
usebin = 0.3;
spkmat = bz_SpktToSpkmat(spikes.times,'dt',usebin,'binsize',usebin,'units','rate');
spkmat.pos = interp1(headdir.timestamps,headdir.data,spkmat.timestamps);
spkmat.InWake = InIntervals(spkmat.timestamps,SleepState.ints.WAKEstate);

for cc = 1:spikes.numcells
    MutInfo.Rate(cc) = mutualinfo(spkmat.data(spkmat.InWake,cc),spkmat.pos(spkmat.InWake));
end

%%
MIthresh = 0.05;
MIforbest = ISIbyHD.MutInf;
MIforbest(MutInfo.Rate<MIthresh) = nan;
[~,bestcell] =max(MIforbest);
[~,sortMI_ISI] = sort(MutInfo.ISI);
[~,sortMutInfo.Rate] = sort(MutInfo.Rate);

%%

for cc = 1:spikes.numcells
    position_norm = headdir;
    position_norm.data = headdir.data-ISIbyHD.fieldpeak(:,:,cc);
    position_norm.data = mod(position_norm.data+pi,2.*pi)-pi;
[ISIbyHD_align(cc)] = bz_ConditionalISI(spikes.times{cc},position_norm,...
    'ints',SleepState.ints.WAKEstate,...
    'showfig',false,'GammaFit',false,'numXbins',30,'numISIbins',100,...
    'normtype','none','Xwin',[-pi pi],'minX',10);
end




ISIbyHD_align_mean = bz_CollapseStruct( ISIbyHD_align(squeeze(ISIbyHD.MutInf)>MIthresh & MutInfo.Rate'>MIthresh),3,'mean',true);
ISIbyHD_align_mean = bz_CollapseStruct( ISIbyHD_align(squeeze(ISIbyHD.MutInf)>MIthresh),3,'mean',true);


ISIbyHD_align = bz_CollapseStruct( ISIbyHD_align,3,'justcat',true);



%%
figure
subplot(3,3,1)
    imagesc(ISIbyHD_align_mean.Dist.Xbins,[1 spikes.numcells],squeeze(log10(ISIbyHD_align.Dist.SpikeRate(:,:,sortMI_ISI)))')
    ylabel('Sort by MIISI')
subplot(3,3,4)
    imagesc(ISIbyHD_align_mean.Dist.Xbins,[1 spikes.numcells],squeeze(log10(ISIbyHD_align.Dist.SpikeRate(:,:,sortMutInfo.Rate)))')
    ylabel('Sort by MIRate')
    
subplot(3,3,7)
    imagesc(ISIbyHD_align_mean.Dist.Xbins,ISIbyHD_align_mean.Dist.Ybins,ISIbyHD_align_mean.Dist.pYX')
    hold on
    imagesc(ISIbyHD_align_mean.Dist.Xbins+2.*pi,ISIbyHD_align_mean.Dist.Ybins,ISIbyHD_align_mean.Dist.pYX')
    plot(ISIbyHD_align_mean.Dist.Xbins,-log10(ISIbyHD_align_mean.Dist.SpikeRate),'r')
    plot(ISIbyHD_align_mean.Dist.Xbins+2.*pi,-log10(ISIbyHD_align_mean.Dist.SpikeRate),'r')
    LogScale('y',10,'nohalf',true)
    ylabel('ISI (s)')
    bz_AddRightRateAxis
    xlabel('Position relative to PF Peak (m)')
    xlim([-pi 3.*pi])
    
subplot(3,3,9)
    plot(log10(MutInfo.Rate),squeeze(log10(ISIbyHD.MutInf)),'.')
    xlabel('MI - rate');ylabel('MI - ISI')
    
subplot(3,3,6)
imagesc(log10(binsizes),[1 spikes.numcells],log10(MutInfo.Rate_BinCompare(sortMutInfo.Rate,:)))
LogScale('x',10)

subplot(3,3,3)
imagesc(log10(binsizes),[1 spikes.numcells],log10(MutInfo.Rate_BinCompare(sortMI_ISI,:)))
LogScale('x',10)
    
NiceSave('HDCoding',figfolder,baseName)


%%
% [excell] = bz_ConditionalISI(spikes.times{bestcell},headdir,...
%     'ints',SleepState.ints.WAKEstate,...
%     'showfig',true,'GammaFit',false,'numXbins',25,'numISIbins',100,'minX',0,...
%     'normtype','none','Xwin',[0 2.*pi]);
