function [ISIbyHD_align,MutInfo,cellISIStats ] = HeadDirectionEncodingAnalysis(basePath,figfolder)
% Date XX/XX/20XX %ISIbyHD_alignGam
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
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = [reporoot,'/Datasets/onProbox/AG_HPC/Achilles_10252013'];
basePath = '/Users/dl2820/Dropbox/Research/Datasets/Mouse12-120807';
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/Mouse24-131213';
%basePath = [reporoot,'/Datasets/onProbox/AG_HPC/Achilles_10252013'];
%basePath = pwd;
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
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
MutInfo.Skaggs = squeeze(nansum(ISIbyHD.Dist.SpikeRate.*log2(ISIbyHD.Dist.SpikeRate./meanrate_full).*ISIbyHD.Dist.pX,2));
MutInfo.Skaggs_sec = MutInfo.Skaggs .*squeeze(meanrate);


%%
[~,ISIbyHD.fieldpeak] = max(ISIbyHD.Dist.SpikeRate,[],2);
ISIbyHD.fieldpeak = ISIbyHD.Dist.Xbins(ISIbyHD.fieldpeak);

%%




%%
MutInfo.ISI = squeeze(ISIbyHD.MutInf)';
binsize = 0.1;
dt = 0.05;
spkmat = bz_SpktToSpkmat(spikes.times,'dt',dt,'binsize',binsize,'units','counts');
spkmat.pos = interp1(headdir.timestamps,headdir.data,spkmat.timestamps,'nearest');
spkmat.InWake = InIntervals(spkmat.timestamps,SleepState.ints.WAKEstate);

maxspikes = 15;
for cc = 1:spikes.numcells
    MutInfo.Rate(cc) = mutualinfo(spkmat.data(spkmat.InWake,cc),spkmat.pos(spkmat.InWake));
    CONDXY(cc) = ConditionalHist(spkmat.pos(spkmat.InWake),spkmat.data(spkmat.InWake,cc),...
         'numXbins',20,'Xbounds',[0 2*pi],'numYbins',maxspikes+1,'Ybounds',[0 maxspikes],...
         'Xbinoverlap',2);
end

%% Spike count histogram
MIthresh = 0.3;
HDcells = find(MutInfo.Skaggs>MIthresh);

figure
for cc = 1:length(HDcells)
    subplot(4,3,cc)
        imagesc(CONDXY(1).Xbins,CONDXY(1).Ybins,CONDXY(HDcells(cc)).pYX')
        hold on
        %imagesc(CONDXY(1).Xbins+2.*pi,CONDXY(1).Ybins,CONDXY(HDcells(cc)).pYX')
        axis xy
        %xlim([0 4*pi])
        bz_piTickLabel('x')
        ColorbarWithAxis([0 0.4],'P[s|HD]')
        xlabel('Head Direction');ylabel('Spike Count')
        title(['UID: ',num2str(spikes.UID(HDcells(cc)))])

end
NiceSave('HDCells_SpikeCount',figfolder,baseName)
%%


testcell = 21;



FitEncodingModel_HD(spkmat.data(spkmat.InWake,testcell),spkmat.pos(spkmat.InWake),dt)






%%


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
    
    %Gamma Fit
    cellUID(cc) = spikes.UID(cc);
    GFIDX = find(GammaFit.WAKEstate.cellstats.UID==cellUID(cc));
    if isempty(GFIDX) %Cells with no gamma....
        continue
    end
    cellGamma = GammaFit.WAKEstate.singlecell(GFIDX);
    MutInfo.GSrate(cc) = GammaFit.WAKEstate.sharedfit.GSlogrates(GFIDX);
    MutInfo.GSweight(cc) = GammaFit.WAKEstate.sharedfit.GSweights(GFIDX);
%     [ISIbyHD_alignGam(cc)] = bz_ConditionalISI(spikes.times{cc},position_norm,...
%         'ints',SleepState.ints.WAKEstate,...
%         'showfig',false,'GammaFit',true,'GammaFitParms',cellGamma,...
%         'numXbins',2,'numISIbins',100,...
%         'normtype','none','Xwin',[-pi pi],'minX',10,'Xbinoverlap',2);
    
    
    %Mean ISI Dist and Return Map
    fieldrange = [-0.33*pi 0.33*pi];
    outfieldrange = [-1.33*pi -0.67*pi ; 0.67*pi 1.33*pi];
    infieldtimes = InIntervals(position_norm.data,fieldrange);
    outfieldtimes = InIntervals(position_norm.data,outfieldrange);
    nantimes = isnan(position_norm.data);
    iswake = InIntervals(position_norm.timestamps,SleepState.ints.WAKEstate);
    IDX.timestamps = position_norm.timestamps;
    IDX.states = infieldtimes*2 + outfieldtimes*1 + 2; %Will set infield = 4, outfield = 3;
    IDX.states(nantimes | ~iswake) = 1; %Set nonwake and no HD data to 1;
    IDX.statenames = {'NanTimes','Side','OutField','InField'};
    INT = bz_IDXtoINT(IDX);
    spikestemp.times = spikes.times(cc);
    spikestemp.UID = spikes.UID(cc);
    ISIstats(cc) = bz_ISIStats(spikestemp,'ints',INT,'showfig',false);
    
    [StateConditionalISI(cc)] = bz_ConditionalISI(spikes.times(cc),INT,...
        'ints','input','normtype','none',...
        'GammaFitParms',cellGamma,'GammaFit',true,...
        'showfig',false);
    %keyboard
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
%% ISI IN/OUT
cellISIStats.statenames = fieldnames(INT);

tempstruct = bz_CollapseStruct(ISIstats,3,'justcat');
meanISIhist = bz_CollapseStruct(tempstruct.ISIhist(MutInfo.Skaggs>MIthresh),3,'mean',true);

StateConditionalISI = bz_CollapseStruct(StateConditionalISI,3,'justcat');

cellISIStats.allISIhist = bz_CollapseStruct(tempstruct.ISIhist,3,'justcat',true);
cellISIStats.allJointhist = bz_CollapseStruct(tempstruct.Jointhist,3,'justcat',true);

cellISIStats.Dist  = bz_CollapseStruct(StateConditionalISI.Dist,3,'justcat',true);
cellISIStats.GammaModes = bz_CollapseStruct(StateConditionalISI.GammaModes,3,'justcat',true);
cellISIStats.GammaModes.states =  StateConditionalISI.states(1,:,1);
%MutInfo.cellclass = CellClass;
%%
figure
subplot(3,3,7)
hold on
for ss = 2:4
    plot(meanISIhist.logbins,meanISIhist.(cellISIStats.statenames{ss}).log)
end
legend(cellISIStats.statenames{2:3},'location','southoutside')

for ss = 3:4
subplot(3,3,(ss-3)*3+1)
hold on
    imagesc(meanISIhist.logbins,[0 1],squeeze(cellISIStats.allISIhist.(cellISIStats.statenames{ss}).log(:,:,MutInfo.Skaggs>MIthresh))')
end

for ss = 2:4
subplot(3,3,(ss-2)*3+2)
hold on
    imagesc(meanISIhist.logbins,meanISIhist.logbins,meanISIhist.(cellISIStats.statenames{ss}).return)
end
NiceSave('InOutField',figfolder,baseName)

%%
infield = 4;
outfield = 3;

figure
subplot(3,3,1)
plot(squeeze(1-cellISIStats.GammaModes.GSweights(1,outfield,(MutInfo.Skaggs>MIthresh))),...
    squeeze(1-cellISIStats.GammaModes.GSweights(1,infield,(MutInfo.Skaggs>MIthresh))),'.')
hold on
UnityLine
ylabel('In-Field AR');xlabel('Out-Field AR')

subplot(3,3,2)
plot(squeeze(cellISIStats.GammaModes.GSlogrates(1,outfield,(MutInfo.Skaggs>MIthresh))),...
    squeeze(cellISIStats.GammaModes.GSlogrates(1,infield,(MutInfo.Skaggs>MIthresh))),'.')
hold on
axis tight
UnityLine
ylabel('In-Field GS Rate');xlabel('Out-Field GS Rate')
LogScale('xy',10)

diffASweight = (cellISIStats.GammaModes.ASweights(infield,:,:)-...
    cellISIStats.GammaModes.ASweights(outfield,:,:));
subplot(3,3,3)
hold on
for aa = 1:5
scatter(-cellISIStats.GammaModes.ASlogrates(1,aa,(MutInfo.Skaggs>MIthresh)),...
    log10(cellISIStats.GammaModes.ASCVs(1,aa,(MutInfo.Skaggs>MIthresh))),...
    10*cellISIStats.GammaModes.ASweights(infield,aa,(MutInfo.Skaggs>MIthresh))+eps,...
    squeeze(diffASweight(1,aa,(MutInfo.Skaggs>MIthresh))))
end
colorbar
crameri('berlin','pivot',0)
caxis([-0.1 0.1])

diffAR = (1-cellISIStats.GammaModes.GSweights(1,infield,:))-(1-cellISIStats.GammaModes.GSweights(1,outfield,:));

subplot(3,3,4)
plot(squeeze(diffAR),log10(MutInfo.Skaggs),'.')
axis tight;box off
hold on
plot([0 0],ylim(gca),'k--')
xlabel('AR Modulation');ylabel('Skaggs Info')

subplot(3,3,6)
plot(squeeze(diffAR),log10(MutInfo.ISI),'.')
axis tight;box off
hold on
plot([0 0],ylim(gca),'k--')
xlabel('AR Modulation');ylabel('ISI Info')

subplot(3,3,5)
plot(squeeze(diffAR),log10(MutInfo.Rate),'.')
axis tight;box off
hold on
plot([0 0],ylim(gca),'k--')
xlabel('AR Modulation');ylabel('Rate Info')

subplot(3,3,7)
plot(squeeze(diffAR),log10(MutInfo.GSrate),'.')
axis tight;box off
hold on
plot([0 0],ylim(gca),'k--')
xlabel('AR Modulation');ylabel('GS Rate')
% subplot(2,2,4)
% plot(squeeze(cellISIStats.MIskaggs),squeeze(diffAR),'.')
NiceSave('GSASField',figfolder,baseName)
%%
ISIbyHD_align_mean = bz_CollapseStruct( ISIbyHD_align(MutInfo.Skaggs>MIthresh),3,'mean',true);
numcells = sum(MutInfo.Skaggs>MIthresh);

% try
% ISIbyHD_alignGam_mean = bz_CollapseStruct( ISIbyHD_alignGam(MutInfo.Skaggs>MIthresh),3,'mean',true);
% catch
%     display('Error Mean')
% end
%%
ISIbyHD_align = bz_CollapseStruct( ISIbyHD_align,3,'justcat',true);
%ISIbyHD_alignGam = bz_CollapseStruct( ISIbyHD_alignGam,3,'justcat',true);


%%
% lowWeight = squeeze(sum(ISIbyHD_alignGam.GammaModes.GSweights<0.05,2)<5);
% safeGSrates = MutInfo.Skaggs>MIthresh & lowWeight;
% meanGSrate = nanmean(ISIbyHD_alignGam.GammaModes.GSlogrates(:,:,safeGSrates),3)


%%

%% Get Peak Width
for cc = 1:spikes.numcells
[~,~,MutInfo.peakwidth(cc)] = findpeaks(ISIbyHD_align.Dist.SpikeRate(1,:,cc),...
    ISIbyHD_align.Dist.Xbins(1,:,cc),...
    'NPeaks',1,'SortStr','descend');
end

%%
% figure
% subplot(3,3,3)
%     polarplot(ISIbyHD_alignGam.Dist.Xbins(1,:,1),1-ISIbyHD_alignGam.GammaModes.GSweights,...
%         'k','linewidth',2)
% ax = gca;
% ax.ThetaAxisUnits = 'radians';
% %ax.ThetaDir = 'clockwise';
% ax.ThetaZeroLocation = 'top';
% rlim([0 1])
% 
% 
% subplot(3,3,6)
%     polarplot(ISIbyHD_alignGam.Dist.Xbins(1,:,1),ISIbyHD_alignGam_mean.GammaModes.GSlogrates,...
%         'k','linewidth',2)
% ax = gca;
% ax.ThetaAxisUnits = 'radians';
% %ax.ThetaDir = 'clockwise';
% ax.ThetaZeroLocation = 'top';
% rlim([-0.3 1])
% 
% NiceSave('PolarGSAR',figfolder,baseName)
%%

% MIforbest = ISIbyHD.MutInf;
% MIforbest(MutInfo.Rate<MIthresh) = nan;
% [~,bestcell] =max(MIforbest);
[~,sortMutInfo.ISI] = sort(MutInfo.ISI);
[~,sortMutInfo.Rate] = sort(MutInfo.Rate);
[~,sortMutInfo.Skaggs] = sort(MutInfo.Skaggs);
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
% subplot(3,3,6)
%     imagesc(ISIbyHD_alignGam.Dist.Xbins(1,:,1),[1 spikes.numcells],...
%         squeeze((ISIbyHD_alignGam.GammaModes.GSweights(:,:,sortMutInfo.Skaggs)))')
%     ylabel('Sort by MIRate')
%     hold on
%     plot(ISIbyHD_align_mean.Dist.Xbins([1 end]),spikes.numcells-numcells-0.5.*[1 1],'r')
    
    
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
 
% try
% subplot(3,3,4)
% plot(ISIbyHD_alignGam_mean.Dist.Xbins,ISIbyHD_alignGam_mean.GammaModes.GSweights,'k')
% hold on
% plot(ISIbyHD_alignGam_mean.Dist.Xbins+2*pi,ISIbyHD_alignGam_mean.GammaModes.GSweights,'k')
% box off 
% axis tight
% catch
%     display('Error Mean')
% end

    
NiceSave('HDCoding',figfolder,baseName)

%%
figure
subplot(2,2,2)
plot((MutInfo.GSrate),(MutInfo.Skaggs),'.')
xlabel('GS Rate');ylabel('MI')
LogScale('x',10)

subplot(2,2,4)
plot(MutInfo.GSweight,(MutInfo.Skaggs),'.')
xlabel('GS Weight');ylabel('MI')

% subplot(2,2,1)
% hold on
% for cc = 1:spikes.numcells
%     if MutInfo.ISI(cc)<MIthresh & MutInfo.Rate(cc)<MIthresh
%         continue
%     end
%     for mm = 1:5
%     scatter(ISIbyHD_alignGam.Dist.Xbins(1,:,1),...
%         ones(size(ISIbyHD_alignGam.Dist.Xbins(1,:,1))).*ISIbyHD_alignGam.GammaModes.ASlogrates(1,mm,cc),...
%         5,ISIbyHD_alignGam.GammaModes.ASweights(:,mm,cc))
%     end
% end
% LogScale('y',10)
%imagesc(ISIbyHD_alignGam(excell).Dist.pYX')      

%plot(ISIbyHD_alignGam(excell).GammaModes.GSweights)

% subplot(2,2,2)
% plot(ISIbyHD_alignGam_mean.GammaModes.ASweights)
% try
% subplot(2,2,3)
% plot(ISIbyHD_alignGam_mean.GammaModes.GSweights)
% catch
%     display('Error Mean')
% end
NiceSave('HDGSAS',figfolder,baseName)
%%
% [excell] = bz_ConditionalISI(spikes.times{bestcell},headdir,...
%     'ints',SleepState.ints.WAKEstate,...
%     'showfig',true,'GammaFit',false,'numXbins',25,'numISIbins',100,'minX',0,...
%     'normtype','none','Xwin',[0 2.*pi]);
