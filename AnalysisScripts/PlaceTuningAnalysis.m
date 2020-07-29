function [ISIbyPOS_norm,MutInfo,cellISIStats ] = PlaceTuningAnalysis(basePath,figfolder)
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
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/Cicero_09102014';
%basePath = [reporoot,'/Datasets/onProbox/AG_HPC/Achilles_10252013'];
%basePath = pwd;
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
statecolors = {'k','b','r',[0.6 0.6 0.6]};

[celltypes,~,typeidx] = unique(CellClass.label);
cellcolor = {'k','r'};


position = bz_LoadBehavior( basePath,'position' );

%%
nantimes = isnan(position.position.lin);
position.data = position.position.lin(~(nantimes));
%Possible here: remove drops for more than a... second?
minjump = 1;
numXbins = 100;
[ position.data ] = NanPadJumps( position.timestamps(~(nantimes)),position.data,minjump );
position.data = interp1(position.timestamps(~(nantimes)),position.data,position.timestamps);
[ISIbyPOS] = bz_ConditionalISI(spikes.times,position,...
    'ints',position.Epochs.MazeEpoch,...
    'showfig',false,'GammaFit',false,'numXbins',numXbins,'numISIbins',100,...
    'normtype','none','Xwin',[0 3],'Xbinoverlap',6,'minX',15);


%%
% excell = 44;
% figure
% subplot(2,2,1)
% plot(firingMaps.xbins{excell}{1},firingMaps.rateMaps{excell}{1})
% hold on
% plot(ISIbyPOS.Dist.Xbins(1,:,excell),ISIbyPOS.Dist.SpikeRate(1,:,excell),'k')
% 
% subplot(2,2,3)
%     imagesc(ISIbyPOS.Dist.Xbins(1,:,excell),ISIbyPOS.Dist.Ybins(1,:,excell),ISIbyPOS.Dist.pYX(:,:,excell)')
%     hold on
%     plot(ISIbyPOS.Dist.Xbins(1,:,excell),-log10(ISIbyPOS.Dist.SpikeRate(1,:,excell)),'r')
%     LogScale('y',10,'nohalf',true)
%     ylabel('ISI (s)')
%     bz_AddRightRateAxis
%     xlabel('Position relative to PF Peak (m)')
%%
MutInfo.numspks = squeeze(sum(ISIbyPOS.Dist.Xhist));
%%
meanrate = nansum(ISIbyPOS.Dist.SpikeRate.*ISIbyPOS.Dist.pX,2);
meanrate_full = repmat(meanrate,1,numXbins,1);
MutInfo.Skaggs = squeeze(nansum(ISIbyPOS.Dist.SpikeRate.*log2(ISIbyPOS.Dist.SpikeRate./meanrate_full).*ISIbyPOS.Dist.pX,2));
MutInfo.Skaggs_sec = MutInfo.Skaggs .*squeeze(meanrate);

MutInfo.meanrate = squeeze(meanrate);

%%
[~,ISIbyPOS.fieldpeak] = max(ISIbyPOS.Dist.SpikeRate,[],2);
ISIbyPOS.fieldpeak = ISIbyPOS.Dist.Xbins(ISIbyPOS.fieldpeak);


%% Compare ISI MI vs Rate MI


binsizes = logspace(-2.5,1.5,25);
for bb = 1:length(binsizes)
    bz_Counter(bb,length(binsizes),'Bin')
    spkmat = bz_SpktToSpkmat(spikes.times,'dt',binsizes(bb),'binsize',binsizes(bb),...
        'win',position.Epochs.MazeEpoch,'units','rate');
    spkmat.pos = interp1(position.timestamps,position.data,spkmat.timestamps);

    for cc = 1:spikes.numcells
        MutInfo.Rate_BinCompare(cc,bb) = mutualinfo(spkmat.data(:,cc),spkmat.pos);
    end
    MutInfo.Rate_SkaggsCompare(bb) = corr(MutInfo.Skaggs,MutInfo.Rate_BinCompare(:,bb));
end
MutInfo.binsizes = binsizes;



%% Compare ISI MI vs Rate MI

MutInfo.ISI = squeeze(ISIbyPOS.MutInf);
spkmat = bz_SpktToSpkmat(spikes.times,'dt',0.3,'binsize',0.3,'win',position.Epochs.MazeEpoch,'units','rate');
spkmat.pos = interp1(position.timestamps,position.data,spkmat.timestamps);

for cc = 1:spikes.numcells
    MutInfo.Rate(cc,1) = mutualinfo(spkmat.data(:,cc),spkmat.pos);
    

end

%%
figure
subplot(2,2,1)
plot(log10(MutInfo.numspks),log10(MutInfo.ISI),'.')
xlabel('numspks');ylabel('MI ISI')

subplot(2,2,2)
plot(log10(MutInfo.numspks),log10(MutInfo.Skaggs),'.')
xlabel('# Spks on Track');ylabel('Skaggs Info')

subplot(2,2,3)
scatter(log10(MutInfo.Skaggs),log10(MutInfo.ISI),5,log10(MutInfo.numspks))
%%
% figure
% subplot(2,2,1)
% plot(MutInfo.GSrate,log10(MutInfo.ISI),'.')
% 
% subplot(2,2,2)
% plot(MutInfo.GSweight,log10(MutInfo.ISI),'.')
%%
MIthresh = 0.03;
MIforbest = ISIbyPOS.MutInf;
MIforbest(MutInfo.Rate<MIthresh) = nan;
[~,bestcell] =max(MIforbest);
[~,sortMI_ISI] = sort(MutInfo.ISI);
[~,sortMutInfo.Rate] = sort(MutInfo.Rate);


%%

for cc = 1:spikes.numcells
    position_norm = position;
    position_norm.data = position.data-ISIbyPOS.fieldpeak(:,:,cc);
[ISIbyPOS_norm(cc)] = bz_ConditionalISI(spikes.times{cc},position_norm,...
    'ints',position_norm.Epochs.MazeEpoch,...
    'showfig',false,'GammaFit',false,'numXbins',120,'numISIbins',100,...
    'normtype','none','Xwin',[-1.5 1.5],'minX',10,'Xbinoverlap',5);
end



%MIthresh = 0.013;
ISIbyPOS_norm_mean = bz_CollapseStruct( ISIbyPOS_norm(squeeze(ISIbyPOS.MutInf)>MIthresh & MutInfo.Rate>MIthresh),3,'mean',true);
%ISIbyPOS_norm_mean = bz_CollapseStruct( ISIbyPOS_norm(squeeze(ISIbyPOS.MutInf)>MIthresh),3,'mean',true);


ISIbyPOS_norm = bz_CollapseStruct( ISIbyPOS_norm,3,'justcat',true);

%% Get gamma mode
GammaFit = bz_LoadCellinfo(basePath,'GammaFit');

MutInfo.GSrate = nan(1,spikes.numcells);
MutInfo.GSweight = nan(1,spikes.numcells);
for cc = 1:spikes.numcells
    cellUID(cc) = spikes.UID(cc);
    GFIDX = find(GammaFit.WAKEstate.cellstats.UID==cellUID(cc));
    if isempty(GFIDX)
        continue
    end
    %cellGamma = GammaFit.WAKEstate.singlecell(GFIDX);
    MutInfo.GSrate(cc) = GammaFit.WAKEstate.sharedfit.GSlogrates(GFIDX);
    MutInfo.GSweight(cc) = GammaFit.WAKEstate.sharedfit.GSweights(GFIDX);
end


%%
MutInfo.cellclass = CellClass;
%%
figure
subplot(3,3,1)
    imagesc(ISIbyPOS_norm_mean.Dist.Xbins,[1 spikes.numcells],squeeze(log10(ISIbyPOS_norm.Dist.SpikeRate(:,:,sortMI_ISI)))')
    ylabel('Sort by MIISI')
subplot(3,3,4)
    imagesc(ISIbyPOS_norm_mean.Dist.Xbins,[1 spikes.numcells],squeeze(log10(ISIbyPOS_norm.Dist.SpikeRate(:,:,sortMutInfo.Rate)))')
    ylabel('Sort by MIRate')
    
subplot(3,3,7)
    imagesc(ISIbyPOS_norm_mean.Dist.Xbins,ISIbyPOS_norm_mean.Dist.Ybins,ISIbyPOS_norm_mean.Dist.pYX')
    hold on
    plot(ISIbyPOS_norm_mean.Dist.Xbins,-log10(ISIbyPOS_norm_mean.Dist.SpikeRate),'r')
    LogScale('y',10,'nohalf',true)
    ylabel('ISI (s)')
    bz_AddRightRateAxis
    xlabel('Position relative to PF Peak (m)')
    

subplot(3,3,9)
    scatter(log10(MutInfo.Rate),squeeze(log10(ISIbyPOS.MutInf)),5,(MutInfo.GSrate))
    xlabel('MI - rate');ylabel('MI - ISI')
    
subplot(3,3,6)
imagesc(log10(binsizes),[1 spikes.numcells],log10(MutInfo.Rate_BinCompare(sortMutInfo.Rate,:)))
LogScale('x',10)

subplot(3,3,3)
imagesc(log10(binsizes),[1 spikes.numcells],log10(MutInfo.Rate_BinCompare(sortMI_ISI,:)))
LogScale('x',10)

NiceSave('PlaceCoding',figfolder,baseName)


%% Buzcode Placefield functions: Cells only with extracted fields
positions = [position.timestamps position.data];
firingMaps = bz_firingMapAvg(positions,spikes);
[placeFieldStats] = bz_findPlaceFields1D('firingMaps',firingMaps,'basepath',basePath);

%% 
%spikestemp
clear tempstruct
clear IDX stateIDX INT
for cc = 1:spikes.numcells
    bz_Counter(cc,spikes.numcells,'Cell')
    if isnan(placeFieldStats.mapStats{cc}{1}.fieldX)
        if cc==spikes.numcells
            tempstruct(cc).UID = [];
            passINT(cc).NanTimesstate = [];
        end
        %tempstruct(cc).GSrate = nan;
        %tempstruct(cc).GSweight = nan;
        continue
    end
    
    fieldrange = firingMaps.xbins{cc}{1}(placeFieldStats.mapStats{cc}{1}.fieldX);
    if length(placeFieldStats.mapStats{cc}{1}.x)==1
        infieldtimes = position.data>=fieldrange(1,1) & position.data<=fieldrange(1,2);
    elseif length(placeFieldStats.mapStats{cc}{1}.x)==2
        infieldtimes = (position.data>=fieldrange(1,1) & position.data<=fieldrange(1,2)) |...
            (position.data>=fieldrange(2,1) & position.data<=fieldrange(2,2));
    end


    IDX.timestamps = position.timestamps;
    IDX.states = infieldtimes+2;
    IDX.states(nantimes) = 1;
    IDX.statenames = {'NanTimes','OutField','InField'};
    INT = bz_IDXtoINT(IDX);
    INT.preWAKEstate = SleepState.ints.WAKEstate(SleepState.ints.WAKEstate(:,2)<position.Epochs.PREEpoch(2),:);
    INT.postWAKEstate = SleepState.ints.WAKEstate(SleepState.ints.WAKEstate(:,1)>position.Epochs.POSTEpoch(1),:);
    INT.preNREM = SleepState.ints.NREMstate(SleepState.ints.NREMstate(:,2)<position.Epochs.PREEpoch(2),:);
    INT.postNREM = SleepState.ints.NREMstate(SleepState.ints.NREMstate(:,1)>position.Epochs.POSTEpoch(1),:);
    passINT(cc) = INT;
    %PRE NREM POST NREM
    %PRE WAKE POST WAKE
    %RUN FIELD
    %RUN OUT FIELD
    %TRACK NORUN
    spikestemp.times = spikes.times(cc);
    spikestemp.UID = spikes.UID(cc);
    tempstruct(cc) = bz_ISIStats(spikestemp,'ints',INT,'showfig',false);
    tempstruct(cc).UID = spikes.UID(cc);
    %tempstruct(cc).GSrate = GammaFit.WAKEstate.sharedfit.GSlogrates(GFIDX);
    %tempstruct(cc).GSweight  = GammaFit.WAKEstate.sharedfit.GSweights(GFIDX);
    
    %stateIDX(cc) = bz_INTtoIDX(INT,'sf',10);
    %stateIDX(cc).data = stateIDX(cc).states;
end

clear StateConditionalISI
for cc = 1:spikes.numcells
    bz_Counter(cc,spikes.numcells,'Cell')
    
        cellUID(cc) = spikes.UID(cc);
    GFIDX = find(GammaFit.WAKEstate.cellstats.UID==cellUID(cc));
    if isempty(GFIDX) || isempty(passINT(cc).NanTimesstate)
        %tempstruct(cc).UID = spikes.UID(cc);
        tempstruct(cc).GSrate_all = nan;
        %tempstruct(cc).GSweight = nan;
        if cc==spikes.numcells
            StateConditionalISI(cc).states = [];
        end
        continue
    end
    tempstruct(cc).GSrate = GammaFit.WAKEstate.sharedfit.GSlogrates(GFIDX);
    tempstruct(cc).GSweight  = GammaFit.WAKEstate.sharedfit.GSweights(GFIDX);
    %StateConditionalISI(cc).GSrate = GammaFit.WAKEstate.sharedfit.GSlogrates(GFIDX);
    %StateConditionalISI(cc).GSweight  = GammaFit.WAKEstate.sharedfit.GSweights(GFIDX);
    tempstruct(cc).MIskaggs = MutInfo.Skaggs(cc);
    tempstruct(cc).MIISI = MutInfo.ISI(cc);
    
    cellGamma = GammaFit.WAKEstate.singlecell(GFIDX);
    %stateIDX(cc).data = stateIDX(cc).states;
    [StateConditionalISI(cc)] = bz_ConditionalISI(spikes.times(cc),passINT(cc),...
        'ints','input','normtype','none',...
        'GammaFitParms',cellGamma,'GammaFit',true,...
        'showfig',false);
    

    %keyboard
end

[~,sortGSrate] = sort([tempstruct(:).GSrate_all]);
tempstruct = bz_CollapseStruct(tempstruct(sortGSrate),3,'justcat');
StateConditionalISI = bz_CollapseStruct(StateConditionalISI(sortGSrate),3,'justcat');
%%
meanISIhist = bz_CollapseStruct(tempstruct.ISIhist,3,'mean',true);
cellISIStats.allISIhist = bz_CollapseStruct(tempstruct.ISIhist,3,'justcat',true);
cellISIStats.allJointhist = bz_CollapseStruct(tempstruct.Jointhist,3,'justcat',true);

cellISIStats.GSrate = bz_CollapseStruct(tempstruct.GSrate,3,'justcat',true);
cellISIStats.GSweight = bz_CollapseStruct(tempstruct.GSweight,3,'justcat',true);
cellISIStats.MIskaggs = bz_CollapseStruct(tempstruct.MIskaggs,3,'justcat',true);
cellISIStats.MIISI = bz_CollapseStruct(tempstruct.MIISI,3,'justcat',true);
cellISIStats.statenames = fieldnames(INT);

cellISIStats.GammaModes = bz_CollapseStruct(StateConditionalISI.GammaModes,3,'justcat',true);
cellISIStats.GammaModes.states =  StateConditionalISI.states(1,:,1);
%%
figure
subplot(2,2,1)
plot(squeeze(1-cellISIStats.GammaModes.GSweights(1,2,:)),squeeze(1-cellISIStats.GammaModes.GSweights(1,3,:)),'.')
hold on
UnityLine
ylabel('In-Field AR');xlabel('Out-Field AR')

subplot(2,2,2)
plot(squeeze(cellISIStats.GammaModes.GSlogrates(1,2,:)),squeeze(cellISIStats.GammaModes.GSlogrates(1,3,:)),'.')
hold on
axis tight
UnityLine
ylabel('In-Field GS Rate');xlabel('Out-Field GS Rate')
LogScale('xy',10)

diffASweight = (cellISIStats.GammaModes.ASweights(3,:,:)-cellISIStats.GammaModes.ASweights(2,:,:));
subplot(2,2,3)
hold on
for aa = 1:5
scatter(-cellISIStats.GammaModes.ASlogrates(1,aa,:),log10(cellISIStats.GammaModes.ASCVs(1,aa,:)),10*cellISIStats.GammaModes.ASweights(3,aa,:)+eps,...
    squeeze(diffASweight(1,aa,:)))
end
colorbar
crameri('berlin','pivot',0)
caxis([-0.1 0.1])

diffAR = (1-cellISIStats.GammaModes.GSweights(1,3,:))-(1-cellISIStats.GammaModes.GSweights(1,2,:));

% subplot(2,2,4)
% plot(squeeze(cellISIStats.MIskaggs),squeeze(diffAR),'.')
NiceSave('GSASField',figfolder,baseName)
%%


figure
subplot(3,3,7)
hold on
for ss = 2:3
    plot(meanISIhist.logbins,meanISIhist.(cellISIStats.statenames{ss}).log)
end
legend(cellISIStats.statenames{2:3},'location','southoutside')

subplot(3,3,8)
hold on
for ss = 4:5
    plot(meanISIhist.logbins,meanISIhist.(cellISIStats.statenames{ss}).log)
end
legend(cellISIStats.statenames{4:5},'location','southoutside')

subplot(3,3,9)
hold on
for ss = 6:7
    plot(meanISIhist.logbins,meanISIhist.(cellISIStats.statenames{ss}).log)
end
legend(cellISIStats.statenames{6:7},'location','southoutside')

for ss = 2:3
subplot(3,3,(ss-2)*3+1)
hold on
    imagesc(meanISIhist.logbins,[0 1],squeeze(cellISIStats.allISIhist.(cellISIStats.statenames{ss}).log)')
end
%legend(cellISIStats.statenames{1:3})

for ss = 4:5
subplot(3,3,(ss-4)*3+2)
hold on
    imagesc(meanISIhist.logbins,[0 1],squeeze(cellISIStats.allISIhist.(cellISIStats.statenames{ss}).log)')
end

for ss = 6:7
subplot(3,3,(ss-6)*3+3)
hold on
    imagesc(meanISIhist.logbins,[0 1],squeeze(cellISIStats.allISIhist.(cellISIStats.statenames{ss}).log)')
end
%For each cell, find the intervals in the field, calculate in-field runing,
%out-field runining, non-running, ISI distributions
NiceSave('InOutField',figfolder,baseName)
%%
% [excell] = bz_ConditionalISI(spikes.times{bestcell},position,...
%     'ints',position_norm.Epochs.MazeEpoch,...
%     'showfig',true,'GammaFit',false,'numXbins',10,'numISIbins',100,'minX',0);
