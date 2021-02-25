function [ ] = AnalysisXXXXXXXX(basePath,figfolder)
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
%reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
basePath = '/Users/dl2820/Dropbox/Research/Datasets/Cicero_09102014';
%basePath = pwd;
%basePath = fullfile(reporoot,'Datasets/onProbox/AG_HPC/Cicero_09102014');
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%Load Stuff
%sessionInfo = bz_getSessionInfo(basePath);
spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
%CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
states = fieldnames(SleepState.ints);
%states{4} = 'ALL';
%SleepState.ints.ALL = [0 Inf];
%statecolors = {'k','b','r',[0.6 0.6 0.6]};

% try
%     celltypes = CellClass.celltypes;
% catch
%     celltypes = unique(CellClass.label);
% end
% cellcolor = {'k','r'};

%%
GammaFit = bz_LoadCellinfo(basePath,'GammaFit');




%%
load([basePath,'/GammaProcessed1/hmm_out.mat'])



%% Prepare the mode structures
ModeHMM.WAKEstate = WAKEall;
ModeHMM.NREMstate = NREMall;

numcells = length(WAKEall);
spkthresh = 50;
MeanReturn.logbins = linspace(-3,2,50);
%get next ISI (nan for last one in the state)
%Cat all the cells
for ss = 1:2
for cc = 1:numcells

ModeHMM.(states{ss})(cc).next_isi = cellfun(@(prevISI) [prevISI(2:end) nan],ModeHMM.(states{ss})(cc).state_isi,'UniformOutput',false);
ModeHMM.(states{ss})(cc).next_state = cellfun(@(prevState) [prevState(2:end) nan],ModeHMM.(states{ss})(cc).decoded_mode,'UniformOutput',false);
ModeHMM.(states{ss})(cc).next_pMode = cellfun(@(prevProb) [prevProb(2:end,:);nan(1,6)],ModeHMM.(states{ss})(cc).prob_mode,'UniformOutput',false);


ModeHMM.(states{ss})(cc).prev_isi = cellfun(@(prevISI) prevISI(1:end),ModeHMM.(states{ss})(cc).state_isi,'UniformOutput',false);
ModeHMM.(states{ss})(cc).prev_state = cellfun(@(prevState) prevState(1:end),ModeHMM.(states{ss})(cc).decoded_mode,'UniformOutput',false);
ModeHMM.(states{ss})(cc).prev_pMode = cellfun(@(prevProb) prevProb(1:end,:),ModeHMM.(states{ss})(cc).prob_mode,'UniformOutput',false);


ModeHMM.(states{ss})(cc).prev_isi = cat(2,ModeHMM.(states{ss})(cc).prev_isi{:});
ModeHMM.(states{ss})(cc).next_isi = cat(2,ModeHMM.(states{ss})(cc).next_isi{:});

ModeHMM.(states{ss})(cc).prev_state = cat(2,ModeHMM.(states{ss})(cc).prev_state{:});
ModeHMM.(states{ss})(cc).next_state = cat(2,ModeHMM.(states{ss})(cc).next_state{:});
ModeHMM.(states{ss})(cc).state_spk = cat(2,ModeHMM.(states{ss})(cc).state_spk{:});

ModeHMM.(states{ss})(cc).prev_pMode = cat(1,ModeHMM.(states{ss})(cc).prev_pMode{:});
ModeHMM.(states{ss})(cc).next_pMode = cat(1,ModeHMM.(states{ss})(cc).next_pMode{:});

for sm = 1:6
    instate_both = ModeHMM.(states{ss})(cc).prev_state == sm & ModeHMM.(states{ss})(cc).next_state==sm;
    instate_either = ModeHMM.(states{ss})(cc).prev_state == sm | ModeHMM.(states{ss})(cc).next_state==sm;
    
    MeanReturn.(states{ss}).cells(cc).both(:,:,sm) = hist3([log10(ModeHMM.(states{ss})(cc).prev_isi(instate_both))',log10(ModeHMM.(states{ss})(cc).next_isi(instate_both))'],...
        {MeanReturn.logbins,MeanReturn.logbins});
    numspk = sum(sum(MeanReturn.(states{ss}).cells(cc).both(:,:,sm)));
    if numspk<spkthresh
        MeanReturn.(states{ss}).cells(cc).both(:,:,sm) = nan(length(MeanReturn.logbins));
    else
        MeanReturn.(states{ss}).cells(cc).both(:,:,sm) = MeanReturn.(states{ss}).cells(cc).both(:,:,sm)./numspk;
    end
    
    MeanReturn.(states{ss}).cells(cc).either(:,:,sm) = hist3([log10(ModeHMM.(states{ss})(cc).prev_isi(instate_either))',log10(ModeHMM.(states{ss})(cc).next_isi(instate_either))'],...
        {MeanReturn.logbins,MeanReturn.logbins});
    numspk = sum(sum(MeanReturn.(states{ss}).cells(cc).either(:,:,sm)));
    if numspk<spkthresh
        MeanReturn.(states{ss}).cells(cc).either(:,:,sm) = nan(length(MeanReturn.logbins));
    else
        MeanReturn.(states{ss}).cells(cc).either(:,:,sm) =  MeanReturn.(states{ss}).cells(cc).either(:,:,sm)./numspk;
    end
end

end
MeanReturn.(states{ss}).mean = bz_CollapseStruct(MeanReturn.(states{ss}),4,'mean',true);

MeanTransition.(states{ss}) = bz_CollapseStruct(ModeHMM.(states{ss}),3,'mean',false);
MeanTransition.(states{ss}) = MeanTransition.(states{ss}).trans_mat;
end


%Mode intervals
modenames = {'AS1','AS2','AS3','AS4','AS5','GS'};
for ss = 1:2
for cc = 1:numcells
        ModeInts.(states{ss}).cells(cc) = bz_IDXtoINT(ModeHMM.(states{ss})(cc).next_state',...
            'statenames',modenames,'jumptol',1);
%         nexttemp = bz_IDXtoINT(ModeHMM.(states{ss})(cc).next_state',...
%             'statenames',modenames);
        
        %Could just use the spike index above....
        for sm = 1:6
        ModeInts_time.(states{ss}).cells(cc).(modenames{sm}) = ...
            [ModeHMM.(states{ss})(cc).state_spk(ModeInts.(states{ss}).cells(cc).([modenames{sm},'state'])(:,1))' ...
            ModeHMM.(states{ss})(cc).state_spk(ModeInts.(states{ss}).cells(cc).([modenames{sm},'state'])(:,2)+1)'];
        end
end
end
%%
%Build this into bz_PlotModeRaster
GScolor = [0.6 0.4 0];
modecolors = crameri('bamako',5);
%modecolors = [modecolors;GScolor];
modecolors = {modecolors(1,:),modecolors(2,:),modecolors(3,:),modecolors(4,:),modecolors(5,:),...
    GScolor};

%%
checkUID = 5;
cellIDX = find(ismember([ModeHMM.WAKEstate(:).UID],checkUID));
%[ exwin ] = bz_RandomWindowInIntervals( SleepState.ints.WAKEstate,60,1 )
exwin = [33885 33935];
linethick = 2;

figure
subplot(2,3,1)

bz_PlotISIDistModes(GammaFit.WAKEstate,checkUID)
title(['Cell UID: ',num2str(checkUID)])

subplot(4,1,3)
bz_PlotModeRaster(ModeHMM.WAKEstate,ModeInts_time.WAKEstate.cells,cellIDX,modecolors,exwin,linethick);

subplot(4,1,4)
hold on
linewidths = [1 1 1 1 1 2];
for mm = 6:-1:1
bothspikes = [ModeHMM.WAKEstate(cellIDX).state_spk;ModeHMM.WAKEstate(cellIDX).state_spk];
bothpmodes = [ModeHMM.WAKEstate(cellIDX).prev_pMode(:,mm) ModeHMM.WAKEstate(cellIDX).next_pMode(:,mm)]';
%plot(ModeHMM.WAKEstate(checkUID).state_spk,ModeHMM.WAKEstate(checkUID).prev_pMode(:,6))
plot(bothspikes(:),bothpmodes(:),'color',modecolors{mm},'linewidth',linewidths(mm))
end
xlim(exwin)
%bz_ScaleBar('s')
set(gca,'xtick',[])
ylabel('p(state)')

GSints = ModeHMM.WAKEstate(cellIDX).prev_state==6;
subplot(3,3,3)
%hist(ModeHMM.WAKEstate(checkUID).prev_pMode(:,6))
hold on
plot(log10(ModeHMM.WAKEstate(cellIDX).prev_isi(GSints)),ModeHMM.WAKEstate(cellIDX).prev_pMode((GSints),6),...
    '.','color',GScolor,'markersize',1)

plot(log10(ModeHMM.WAKEstate(cellIDX).prev_isi(~GSints)),ModeHMM.WAKEstate(cellIDX).prev_pMode((~GSints),6),...
    'k.','markersize',1)
LogScale('x',10,'exp',true,'nohalf',true)
xlabel('ISI (s)');
ylabel('P(GS)')

subplot(4,3,2)
hold on
histogram(log10(ModeHMM.WAKEstate(cellIDX).prev_isi),'BinEdges', ISIStats.ISIhist.logbins,'facecolor','none')
histogram(log10(ModeHMM.WAKEstate(cellIDX).prev_isi(GSints)),'facecolor',GScolor,'BinEdges', ISIStats.ISIhist.logbins)
LogScale('x',10,'exp',true,'nohalf',true)
xlabel('ISI (s)');ylabel('# Intervals')

subplot(4,3,5)
hold on
histogram((ModeHMM.WAKEstate(cellIDX).prev_pMode(:,6)),'BinEdges',linspace(0,1,20),'facecolor','none')
histogram((ModeHMM.WAKEstate(cellIDX).prev_pMode(GSints,6)),'facecolor',GScolor,'BinEdges', linspace(0,1,20))
%LogScale('x',10,'exp',true,'nohalf',true)
xlabel('p(GS)');ylabel('# Intervals')


NiceSave(['GSPosterior_Cell',num2str(checkUID)],figfolder,baseName)

%%
ss =1
cc = 7;
figure
subplot(3,3,1)
plot(log10(ModeHMM.WAKEstate(cc).prev_isi),log10(ModeHMM.WAKEstate(cc).next_isi),'k.','markersize',1)
xlim([-3 2]);ylim([-3 2])
title('All Spikes')
xlabel('ISI_n');ylabel('ISI_n_+_1')
LogScale('xy',10,'exp',true)

subplot(3,3,7)
imagesc(MeanTransition.(states{ss})([2,4,1,5,3,6],[2,4,1,5,3,6],:))
xlabel('To State');ylabel('From State')
ColorbarWithAxis([0 0.6],'P(Transition)','inclusive',{'','>'})

for sm = 1:6
    %if ss==6
        instate_both = ModeHMM.WAKEstate(cc).prev_state == sm & ModeHMM.WAKEstate(cc).next_state==sm;
    %else
        instate_either = ModeHMM.WAKEstate(cc).prev_state == sm | ModeHMM.WAKEstate(cc).next_state==sm;
    %end
subplot(6,6,(sm-1)*6+3)
    plot(log10(ModeHMM.WAKEstate(cc).prev_isi(instate_either)),log10(ModeHMM.WAKEstate(cc).next_isi(instate_either)),...
        '.','color',[0.5 0.5 0.5],'markersize',0.5)
    hold on
    plot(log10(ModeHMM.WAKEstate(cc).prev_isi(instate_both)),log10(ModeHMM.WAKEstate(cc).next_isi(instate_both)),...
        'k.','markersize',0.5)
     set(gca,'yticklabel',[]);set(gca,'xticklabel',[])
    xlim([-3 2]);ylim([-3 2])
    if sm == 6
        title('GS')
    else
        title(['AS',num2str(sm)])
    end
end

for sm = 1:6
    subplot(6,6,(sm-1)*6+4)
        imagesc(MeanReturn.(states{ss}).mean.cells.both(:,:,sm))
        axis xy
        set(gca,'yticklabel',[]);set(gca,'xticklabel',[])
        
    subplot(6,6,(sm-1)*6+5)
        imagesc(MeanReturn.(states{ss}).mean.cells.either(:,:,sm))
        axis xy
        set(gca,'yticklabel',[]);set(gca,'xticklabel',[])

end
NiceSave('ModeReturnmaps',figfolder,baseName)






%%

%Find times at which the animal crosses the cell's field peak



%%
position = bz_LoadBehavior( basePath,'position' );
placeFieldStats = bz_LoadCellinfo(basePath,'placeFields');

nantimes = isnan(position.position.lin);
position.data = position.position.lin(~(nantimes));
%Possible here: remove drops for more than a... second?
minjump = 1;
numXbins = 100;
[ position.data ] = NanPadJumps( position.timestamps(~(nantimes)),position.data,minjump );
position.data = interp1(position.timestamps(~(nantimes)),position.data,position.timestamps);

positions = [position.timestamps position.data];
firingMaps = bz_firingMapAvg(positions,spikes);

%%
for cc = 1:length(placeFieldStats.UID)
    %bz_Counter(cc,spikes.numcells,'Cell')
    cellUID = placeFieldStats.UID(cc);
    whichcell = find(ismember([ModeHMM.WAKEstate(:).UID],cellUID));
    
    if isnan(placeFieldStats.mapStats{cc}{1}.fieldX)
        display('no field')
        fieldpeak(cc) = nan;
        continue
    end

    peakrate =placeFieldStats.mapStats{cc}{1}.peak;
    %fieldrate
    
    fieldrange = firingMaps.xbins{cc}{1}(placeFieldStats.mapStats{cc}{1}.fieldX);
    fieldpeak(cc) = firingMaps.xbins{cc}{1}(placeFieldStats.mapStats{cc}{1}.x(1));
    infieldtimes = InIntervals(position.data,fieldpeak(cc)+[0 0.1]);

    IDX.timestamps = position.timestamps;
    IDX.states = infieldtimes+2;
    IDX.states(nantimes) = 1;
    IDX.statenames = {'NanTimes','OutField','InField'};
    INT = bz_IDXtoINT(IDX);
    
    %longints = diff(INT.InFieldstate,1,2)<2;
    enterfield = INT.InFieldstate(:,1);
    enterfield([false;diff(enterfield)<2]) = [];




%%
GScolor = [0.6 0.4 0];
modecolors = crameri('bamako',5);
%modecolors = [modecolors;GScolor];
modecolors = {modecolors(1,:),modecolors(2,:),modecolors(3,:),modecolors(4,:),modecolors(5,:),...
    GScolor};

linethick = 2;

figure
subplot(2,3,1.25)
bz_PlotISIDistModes(GammaFit.WAKEstate,cellUID);

subplot(2,2,3)
for tt = 1:length(enterfield)
    temp_spikemodes(tt) = ModeHMM.WAKEstate(whichcell);
    temp_spikemodes(tt).state_spk =  temp_spikemodes(tt).state_spk - enterfield(tt);
    temp_modeintervals(tt) = ModeInts_time.WAKEstate.cells(whichcell);
    temp_modeintervals(tt) = structfun(@(X) X-enterfield(tt),temp_modeintervals(tt),'UniformOutput',false);
end
bz_PlotModeRaster(temp_spikemodes,temp_modeintervals,[1:length(enterfield)],modecolors,[-3 3],linethick);
plot([0 0],ylim(gca),'k--')
title(['Cell UID: ',num2str(cellUID),'. Peak Rate: ',num2str(peakrate)])
ylabel('Trial')


for xx = 1:4
tt = randsample(1:length(enterfield),1); %trial example
exwin = enterfield(tt)+[-1 1];
subplot(4,2,xx*2)
bz_MultiLFPPlot(thetalfp,'timewin',exwin,'LFPmidpoints',3,'scaleLFP',0.5e-3)
hold on
bz_PlotModeRaster(ModeHMM.WAKEstate,ModeInts_time.WAKEstate.cells,whichcell,modecolors,exwin,linethick);
hold on
ylim([-1 4])
%subplot(4,2,4)
end

NiceSave(['PlaceField_Cell',num2str(cellUID)],figfolder,baseName)

end

%%

%%
%% Get PFLocation for sorting
[~,ordermatch] = ismember([ModeHMM.WAKEstate(:).UID],placeFieldStats.UID);
PFlocation = fieldpeak(ordermatch);
% Many-neuron figure

%plotnumcells  = 30;
%plotrates = cellrates.WAKEstate(1:plotnumcells);
[~,cellorder] = sort(PFlocation);

%remove no field
nofield = find(isnan(PFlocation));
%%
%cellorder(ismember(cellorder,nofield))=[];
%%
firsttrial = position.timestamps(find(position.data>1,1));
exwin = firsttrial+[100 150];
figure
subplot(2,1,1)
bz_PlotModeRaster(ModeHMM.WAKEstate,ModeInts_time.WAKEstate.cells,cellorder,modecolors,exwin,linethick)
ylabel('Cell (30), sort by rate')

subplot(4,1,3)
plot(position.timestamps,position.data)
xlim(exwin)




%% Get rates for sorting
[~,ordermatch] = ismember([ModeHMM.NREMstate(:).UID],ISIStats.UID);
cellrates.NREMstate = ISIStats.summstats.NREMstate.meanrate(ordermatch);
%% Many-neuron figure

plotnumcells  = 75;
plotrates = cellrates.NREMstate(1:plotnumcells);
[~,cellorder] = sort(plotrates);
%[ exwin ] = bz_RandomWindowInIntervals( SleepState.ints.WAKEstate,60,1 )
exwin = randsample(SharpWaves.peaktimes,1) + [-2 2]
%exwin = [4980 5025]
linethick = 4
figure
subplot(3,1,2:3)
bz_PlotModeRaster(ModeHMM.NREMstate,ModeInts_time.NREMstate.cells,cellorder,modecolors,exwin,linethick)
hold on
plot(mean(exwin).*[1 1],ylim(gca),'k--')
ylabel('Cell (30), sort by rate')
%NiceSave('SWRExample',figfolder,baseName)

subplot(3,1,1)
bz_MultiLFPPlot(lfp,'timewin',exwin)
hold on
plot(mean(exwin).*[1 1],ylim(gca),'k--')
%plot(lfp.timestamps,lfp.data,'k')
%xlim(exwin)


NiceSave('SWRExample',figfolder,baseName)


%%





%%
ModeInts_all = bz_CollapseStruct(ModeInts,1,'justcat',true);

%%
for ss = 1:2
    for sm = 1:6
        ModeInts_all.(states{ss}).chainlength.(modenames{sm}) = diff(ModeInts_all.(states{ss}).cells.([modenames{sm},'state']),1,2);
        [ModeInts_all.(states{ss}).chainhist.(modenames{sm}),ModeInts_all.(states{ss}).chainbins.(modenames{sm})] = ...
            hist(ModeInts_all.(states{ss}).chainlength.(modenames{sm}),[0:10]);
        ModeInts_all.(states{ss}).chainhist.(modenames{sm}) = ModeInts_all.(states{ss}).chainhist.(modenames{sm})./sum(ModeInts_all.(states{ss}).chainhist.(modenames{sm}));
    end
end

%%
figure
for ss = 1:2
    for sm = 1:6
        subplot(6,2,(sm-1)*2+ss)
            bar(ModeInts_all.(states{ss}).chainbins.(modenames{sm})+1,ModeInts_all.(states{ss}).chainhist.(modenames{sm}))
            box off
            
            if ss == 1
                ylabel(modenames{sm})
            end
            if sm == 1
                title(states{ss})
            end
            if sm == 6
                xlabel('Chain Length (# ISIs)')
            end
    end
end
NiceSave('ModeChainLengths',figfolder,baseName)

%% Single cell example
GScolor = [0.6 0.4 0];
modecolors = crameri('bamako',5);
%modecolors = [modecolors;GScolor];
modecolors = {modecolors(1,:),modecolors(2,:),modecolors(3,:),modecolors(4,:),modecolors(5,:),...
    GScolor};

[ exwin ] = bz_RandomWindowInIntervals( SleepState.ints.WAKEstate,60,1 );

exwin = [3308 3329];

linethick = 10;

cc = 6;
figure
subplot(4,1,1)
bz_PlotModeRaster(ModeHMM.WAKEstate,ModeInts_time.WAKEstate.cells,6,modecolors,exwin,linethick)

subplot(3,3,4)
plot(log10(ModeHMM.WAKEstate(cc).prev_isi),log10(ModeHMM.WAKEstate(cc).next_isi),'k.','markersize',1)
xlim([-3 2]);ylim([-3 2])
title('All Spikes')
xlabel('ISI_n');ylabel('ISI_n_+_1')
LogScale('xy',10,'exp',true)

subplot(5,3,9)
hist(log10(ModeHMM.WAKEstate(cc).prev_isi),60)
box off
LogScale('x',10,'nohalf',true)

for sm = 1:6
    %if ss==6
        instate_both = ModeHMM.WAKEstate(cc).prev_state == sm & ModeHMM.WAKEstate(cc).next_state==sm;
    %else
        instate_either = ModeHMM.WAKEstate(cc).prev_state == sm | ModeHMM.WAKEstate(cc).next_state==sm;
    %end
subplot(6,6,5*6+(sm))
    plot(log10(ModeHMM.WAKEstate(cc).prev_isi(instate_either)),log10(ModeHMM.WAKEstate(cc).next_isi(instate_either)),...
        '.','color',[0.5 0.5 0.5],'markersize',0.5)
    hold on
    plot(log10(ModeHMM.WAKEstate(cc).prev_isi(instate_both)),log10(ModeHMM.WAKEstate(cc).next_isi(instate_both)),...
        'k.','markersize',0.5)
     set(gca,'yticklabel',[]);set(gca,'xticklabel',[])
    xlim([-3 2]);ylim([-3 2])
    if sm == 6
        title('GS','Color',modecolors{sm})
    else
        title(['AS',num2str(sm)],'Color',modecolors{sm})
    end
end




NiceSave(['CellExample_',num2str(cc)],figfolder,baseName)


%% Get rates for sorting
[~,ordermatch] = ismember([ModeHMM.WAKEstate(:).UID],ISIStats.UID);
cellrates.WAKEstate = ISIStats.summstats.WAKEstate.meanrate(ordermatch);
%% Many-neuron figure

plotnumcells  = 30;
plotrates = cellrates.WAKEstate(1:plotnumcells);
[~,cellorder] = sort(plotrates);
%[ exwin ] = bz_RandomWindowInIntervals( SleepState.ints.WAKEstate,60,1 )
exwin = [4980 5025]
linethick = 4
figure
subplot(2,1,1)
bz_PlotModeRaster(ModeHMM.WAKEstate,ModeInts_time.WAKEstate.cells,cellorder,modecolors,exwin,linethick)
ylabel('Cell (30), sort by rate')
NiceSave('MultiCellExample',figfolder,baseName)

%%

%[ exwin ] = bz_RandomWindowInIntervals( SleepState.ints.WAKEstate,60,1 )
exwin = [5005 5012]
linethick = 4
figure
subplot(2,1,1)
bz_PlotModeRaster(ModeHMM.WAKEstate,ModeInts_time.WAKEstate.cells,cellorder,modecolors,exwin,linethick)
ylabel('Cell (30), sort by rate')

exwin = [5009.5 5011]
subplot(2,1,2)
bz_PlotModeRaster(ModeHMM.WAKEstate,ModeInts_time.WAKEstate.cells,cellorder,modecolors,exwin,linethick+2)
ylabel('Cell (30), sort by rate')
ylim([15 31])
NiceSave('MultiCellExample_Zoom',figfolder,baseName)
%ylim(
%%

[ exwin ] = bz_RandomWindowInIntervals( SleepState.ints.WAKEstate,60,1 )

%spiketimes


plotnumcells  = 40;
plotrates = cellrates.WAKEstate(1:plotnumcells);
[~,cellorder] = sort(plotrates);
%rates = 
figure
subplot(2,1,1)
hold on
for cc = 1:plotnumcells
    whichcell = cellorder(cc);
    plotints =  structfun(@(modeints) RestrictInts(modeints,exwin,'inclusive',true),ModeInts_time.WAKEstate.cells(whichcell),'UniformOutput',false);
StateScorePlot( plotints,modecolors,'y',cc,'LineWidth',5)

for sm = 1:6
    %if ss==6
        instate_both = ModeHMM.WAKEstate(whichcell).prev_state == sm & ModeHMM.WAKEstate(whichcell).next_state==sm & ...
            InIntervals(ModeHMM.WAKEstate(whichcell).state_spk',exwin)';
    %else
        instate_either = (ModeHMM.WAKEstate(whichcell).prev_state == sm | ModeHMM.WAKEstate(whichcell).next_state==sm) & ...
            InIntervals(ModeHMM.WAKEstate(whichcell).state_spk',exwin)';
    %end
    plot([ModeHMM.WAKEstate(whichcell).state_spk(instate_either);ModeHMM.WAKEstate(whichcell).state_spk(instate_either)],...
        cc+[zeros(size(ModeHMM.WAKEstate(whichcell).state_spk(instate_either)))-0.4;0.4+zeros(size(ModeHMM.WAKEstate(whichcell).state_spk(instate_either)))],...
        'color',[0.5 0.5 0.5],'linewidth',0.5)
    plot([ModeHMM.WAKEstate(whichcell).state_spk(instate_both);ModeHMM.WAKEstate(whichcell).state_spk(instate_both)],...
        cc+[zeros(size(ModeHMM.WAKEstate(whichcell).state_spk(instate_both)))-0.4;0.4+zeros(size(ModeHMM.WAKEstate(whichcell).state_spk(instate_both)))],...
        'color','k','linewidth',1)
    
end
end
xlim(exwin)
ylim([0 plotnumcells+1])
box off
bz_ScaleBar('s')



% subplot(3,3,1)
% plot(log10(ModeHMM.WAKEstate(cc).prev_isi),log10(ModeHMM.WAKEstate(cc).next_isi),'k.','markersize',1)
% xlim([-3 2]);ylim([-3 2])
% title('All Spikes')
% xlabel('ISI_n');ylabel('ISI_n_+_1')
% LogScale('xy',10,'exp',true)
% 
% subplot(3,3,7)
% imagesc(MeanTransition.(states{ss})([2,4,1,5,3,6],[2,4,1,5,3,6],:))
% xlabel('To State');ylabel('From State')
% ColorbarWithAxis([0 0.6],'P(Transition)','inclusive',{'','>'})
cc = 1;
for sm = 1:6
    %if ss==6
        instate_both = ModeHMM.WAKEstate(cc).prev_state == sm & ModeHMM.WAKEstate(cc).next_state==sm;
    %else
        instate_either = ModeHMM.WAKEstate(cc).prev_state == sm | ModeHMM.WAKEstate(cc).next_state==sm;
    %end
subplot(6,6,3*6+(sm))
    plot(log10(ModeHMM.WAKEstate(cc).prev_isi(instate_either)),log10(ModeHMM.WAKEstate(cc).next_isi(instate_either)),...
        '.','color',[0.5 0.5 0.5],'markersize',0.5)
    hold on
    plot(log10(ModeHMM.WAKEstate(cc).prev_isi(instate_both)),log10(ModeHMM.WAKEstate(cc).next_isi(instate_both)),...
        'k.','markersize',0.5)
     set(gca,'yticklabel',[]);set(gca,'xticklabel',[])
    xlim([-3 2]);ylim([-3 2])
    if sm == 6
        title('GS','Color',modecolors{sm})
    else
        title(['AS',num2str(sm)],'Color',modecolors{sm})
    end
end

for sm = 1:6
    subplot(6,6,4*6+(sm))
        imagesc(MeanReturn.(states{ss}).mean.cells.both(:,:,sm))
        axis xy
        set(gca,'yticklabel',[]);set(gca,'xticklabel',[])
        
    subplot(6,6,5*6+(sm))
        imagesc(MeanReturn.(states{ss}).mean.cells.either(:,:,sm))
        axis xy
        set(gca,'yticklabel',[]);set(gca,'xticklabel',[])

end

NiceSave('ModeRaster',figfolder,baseName)


end
