function [ISIbyHD_align,MutInfo,ISIbyHD_alignGam ] = HeadDirectionTuningAnalysis(basePath,figfolder)
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
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/Mouse24-131213';
%basePath = [reporoot,'/Datasets/onProbox/AG_HPC/Achilles_10252013'];
%basePath = pwd;
%figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%Load Stuff
sessionInfo = bz_getSessionInfo(basePath);
spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true,'region','THAL');
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
numXbins = 100;
[ISIbyHD] = bz_ConditionalISI(spikes.times,headdir,...
    'ints',SleepState.ints.WAKEstate,...
    'showfig',false,'GammaFit',false,'numXbins',numXbins,'numISIbins',100,...
    'normtype','none','Xwin',[0 2.*pi],'Xbinoverlap',3);

%% Get Tuning Specificity (Skaggs)

meanrate = nansum(ISIbyHD.Dist.SpikeRate.*ISIbyHD.Dist.pX,2);
meanrate_full = repmat(meanrate,1,numXbins,1);
MutInfo.SkaggsInf = squeeze(nansum(ISIbyHD.Dist.SpikeRate.*log2(ISIbyHD.Dist.SpikeRate./meanrate_full).*ISIbyHD.Dist.pX,2));
MutInfo.SkaggsInf_sec = MutInfo.SkaggsInf .*squeeze(meanrate);


%%
[~,ISIbyHD.fieldpeak] = max(ISIbyHD.Dist.SpikeRate,[],2);
ISIbyHD.fieldpeak = ISIbyHD.Dist.Xbins(ISIbyHD.fieldpeak);

%%

%% Compare ISI MI vs Rate MI


binsizes = logspace(-2.5,1.5,25);
for bb = 1:length(binsizes)
    bz_Counter(bb,length(binsizes),'Bin')
spkmat = bz_SpktToSpkmat(spikes.times,'dt',binsizes(bb),'binsize',binsizes(bb),'units','rate');
spkmat.pos = interp1(headdir.timestamps,headdir.data,spkmat.timestamps,'nearest');
spkmat.InWake = InIntervals(spkmat.timestamps,SleepState.ints.WAKEstate);

for cc = 1:spikes.numcells
    MutInfo.Rate_BinCompare(cc,bb) = mutualinfo(spkmat.data(spkmat.InWake,cc),spkmat.pos(spkmat.InWake));
    
end
MutInfo.Rate_SkaggsCompare(bb) = corr(MutInfo.SkaggsInf,MutInfo.Rate_BinCompare(:,bb));
end



%%
MutInfo.ISI = squeeze(ISIbyHD.MutInf)';
usebin = 0.1;
spkmat = bz_SpktToSpkmat(spikes.times,'dt',usebin,'binsize',usebin,'units','rate');
spkmat.pos = interp1(headdir.timestamps,headdir.data,spkmat.timestamps,'nearest');
spkmat.InWake = InIntervals(spkmat.timestamps,SleepState.ints.WAKEstate);

for cc = 1:spikes.numcells
    MutInfo.Rate(cc) = mutualinfo(spkmat.data(spkmat.InWake,cc),spkmat.pos(spkmat.InWake));
end


%% Get gamma mode
GammaFit = bz_LoadCellinfo(basePath,'GammaFit');



%%
MutInfo.GSrate = nan(1,spikes.numcells);
clear ISIbyHD_align
for cc = 1:spikes.numcells
    bz_Counter(cc,spikes.numcells,'Cell')    
   
    position_norm = headdir;
    position_norm.data = headdir.data-ISIbyHD.fieldpeak(:,:,cc);
    position_norm.data = mod(position_norm.data+pi,2.*pi)-pi;
[ISIbyHD_align(cc)] = bz_ConditionalISI(spikes.times{cc},position_norm,...
    'ints',SleepState.ints.WAKEstate,...
    'showfig',false,'GammaFit',false,...
    'numXbins',50,'numISIbins',100,...
    'normtype','none','Xwin',[-pi pi],'minX',10,'Xbinoverlap',1);

    cellUID(cc) = spikes.UID(cc);
    GFIDX = find(GammaFit.WAKEstate.cellstats.UID==cellUID(cc));
    if isempty(GFIDX) %Cells with no gamma....
        continue
    end
    cellGamma = GammaFit.WAKEstate.singlecell(GFIDX);
    MutInfo.GSrate(cc) = GammaFit.WAKEstate.sharedfit.GSlogrates(GFIDX);
    MutInfo.GSweight(cc) = GammaFit.WAKEstate.sharedfit.GSweights(GFIDX);
    [ISIbyHD_alignGam(cc)] = bz_ConditionalISI(spikes.times{cc},position_norm,...
        'ints',SleepState.ints.WAKEstate,...
        'showfig',false,'GammaFit',true,'GammaFitParms',cellGamma,...
        'numXbins',25,'numISIbins',100,...
        'normtype','none','Xwin',[-pi pi],'minX',10,'Xbinoverlap',2);
end

%% Average all HD cells
MIthresh = 0.3;
% try
% ISIbyHD_align_mean = bz_CollapseStruct( ISIbyHD_align(MutInfo.ISI>MIthresh & MutInfo.Rate>MIthresh),3,'mean',true);
% ISIbyHD_alignGam_mean = bz_CollapseStruct( ISIbyHD_alignGam(MutInfo.ISI>MIthresh & MutInfo.Rate>MIthresh),3,'mean',true);
% numcells = sum(MutInfo.ISI>MIthresh & MutInfo.Rate>MIthresh);
% catch
% ISIbyHD_align_mean = bz_CollapseStruct( ISIbyHD_align(MutInfo.Rate>MIthresh),3,'mean',true);
% ISIbyHD_alignGam_mean = bz_CollapseStruct( ISIbyHD_alignGam(MutInfo.Rate>MIthresh),3,'mean',true);
% numcells = sum(MutInfo.Rate>MIthresh);
% end
%%
%MutInfo.cellclass = CellClass;
%%
ISIbyHD_align_mean = bz_CollapseStruct( ISIbyHD_align(MutInfo.SkaggsInf>MIthresh),3,'mean',true);
numcells = sum(MutInfo.SkaggsInf>MIthresh);

try
ISIbyHD_alignGam_mean = bz_CollapseStruct( ISIbyHD_alignGam(MutInfo.SkaggsInf>MIthresh),3,'mean',true);
catch
    display('Error Mean')
end
%%
ISIbyHD_align = bz_CollapseStruct( ISIbyHD_align,3,'justcat',true);
ISIbyHD_alignGam = bz_CollapseStruct( ISIbyHD_alignGam,3,'justcat',true);




%% Get Peak Width
for cc = 1:spikes.numcells
[~,~,MutInfo.peakwidth(cc)] = findpeaks(ISIbyHD_align.Dist.SpikeRate(1,:,cc),...
    ISIbyHD_align.Dist.Xbins(1,:,cc),...
    'NPeaks',1,'SortStr','descend');
end

%%
%  figure
% plot(MutInfo.peakwidth,MutInfo.SkaggsInf,'.')
%%

% MIforbest = ISIbyHD.MutInf;
% MIforbest(MutInfo.Rate<MIthresh) = nan;
% [~,bestcell] =max(MIforbest);
[~,sortMutInfo.ISI] = sort(MutInfo.ISI);
[~,sortMutInfo.Rate] = sort(MutInfo.Rate);
[~,sortMutInfo.Skaggs] = sort(MutInfo.SkaggsInf);
[~,sortMutInfo.peakwidth] = sort(MutInfo.peakwidth);

%% Rate bin
figure
subplot(3,3,9)
    plot(log10(MutInfo.Rate),squeeze(log10(ISIbyHD.MutInf)),'.')
    xlabel('MI - rate');ylabel('MI - ISI')
    
subplot(3,3,6)
imagesc(log10(binsizes),[1 spikes.numcells],log10(MutInfo.Rate_BinCompare(sortMutInfo.Skaggs,:)))
LogScale('x',10)

subplot(3,3,3)
imagesc(log10(binsizes),[1 spikes.numcells],log10(MutInfo.Rate_BinCompare(sortMutInfo.ISI,:)))
LogScale('x',10)

subplot(3,3,5)    
    plot(log10(binsizes),MutInfo.Rate_SkaggsCompare)
    LogScale('x',10)
    
subplot(3,3,1)

    imagesc(ISIbyHD_align_mean.Dist.Xbins,[1 spikes.numcells],squeeze(log10(ISIbyHD_align.Dist.SpikeRate(:,:,sortMutInfo.ISI)))')
    ylabel('Sort by MIISI')
    hold on
    plot(ISIbyHD_align_mean.Dist.Xbins([1 end]),spikes.numcells-numcells-0.5.*[1 1],'r')
subplot(3,3,4)
    imagesc(ISIbyHD_align_mean.Dist.Xbins,[1 spikes.numcells],squeeze(log10(ISIbyHD_align.Dist.SpikeRate(:,:,sortMutInfo.Rate)))')
    ylabel('Sort by MIRate')
    hold on
    plot(ISIbyHD_align_mean.Dist.Xbins([1 end]),spikes.numcells-numcells-0.5.*[1 1],'r')
 NiceSave('RateBin',figfolder,baseName)   
%%
figure

subplot(3,3,3)
    imagesc(ISIbyHD_align_mean.Dist.Xbins,[1 spikes.numcells],squeeze(log10(ISIbyHD_align.Dist.SpikeRate(:,:,sortMutInfo.Skaggs)))')
    ylabel('Sort by MIRate')
    hold on
    plot(ISIbyHD_align_mean.Dist.Xbins([1 end]),spikes.numcells-numcells-0.5.*[1 1],'r')
subplot(3,3,6)
    imagesc(ISIbyHD_alignGam.Dist.Xbins(1,:,1),[1 spikes.numcells],...
        squeeze((ISIbyHD_alignGam.GammaModes.GSweights(:,:,sortMutInfo.Skaggs)))')
    ylabel('Sort by MIRate')
    hold on
    plot(ISIbyHD_align_mean.Dist.Xbins([1 end]),spikes.numcells-numcells-0.5.*[1 1],'r')
    
    
subplot(3,3,1)
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
 
try
subplot(3,3,4)
plot(ISIbyHD_alignGam_mean.Dist.Xbins,ISIbyHD_alignGam_mean.GammaModes.GSweights,'k')
hold on
plot(ISIbyHD_alignGam_mean.Dist.Xbins+2*pi,ISIbyHD_alignGam_mean.GammaModes.GSweights,'k')
box off 
axis tight
catch
    display('Error Mean')
end

    
NiceSave('HDCoding',figfolder,baseName)

%%
figure
subplot(2,2,2)
plot((MutInfo.GSrate),(MutInfo.SkaggsInf),'.')
xlabel('GS Rate');ylabel('MI')
LogScale('x',10)

subplot(2,2,4)
plot(MutInfo.GSweight,(MutInfo.SkaggsInf),'.')
xlabel('GS Weight');ylabel('MI')

subplot(2,2,1)
hold on
for cc = 1:spikes.numcells
    if MutInfo.ISI(cc)<MIthresh & MutInfo.Rate(cc)<MIthresh
        continue
    end
    for mm = 1:5
    scatter(ISIbyHD_alignGam.Dist.Xbins(1,:,1),...
        ones(size(ISIbyHD_alignGam.Dist.Xbins(1,:,1))).*ISIbyHD_alignGam.GammaModes.ASlogrates(1,mm,cc),...
        5,ISIbyHD_alignGam.GammaModes.ASweights(:,mm,cc))
    end
end
LogScale('y',10)
%imagesc(ISIbyHD_alignGam(excell).Dist.pYX')      

%plot(ISIbyHD_alignGam(excell).GammaModes.GSweights)

% subplot(2,2,2)
% plot(ISIbyHD_alignGam_mean.GammaModes.ASweights)
try
subplot(2,2,3)
plot(ISIbyHD_alignGam_mean.GammaModes.GSweights)
catch
    display('Error Mean')
end
NiceSave('HDGSAS',figfolder,baseName)
%%
% [excell] = bz_ConditionalISI(spikes.times{bestcell},headdir,...
%     'ints',SleepState.ints.WAKEstate,...
%     'showfig',true,'GammaFit',false,'numXbins',25,'numISIbins',100,'minX',0,...
%     'normtype','none','Xwin',[0 2.*pi]);
