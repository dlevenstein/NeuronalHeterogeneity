function [ ] = Analysis20190321(basePath,figfolder)
% Date 03/20/19
%
%Questions: Ripple coupling and modulation of the ISI distribution
%
%Plots
%-Ripple-band spike-phase coupling and I(power;ISI) for each cell/channel
%-
%
%% Load Header
%Initiate Paths
reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/Cicero_09102014';
basePath = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onDesktop/AG_HPC/Achilles_10252013';
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
%states{4} = 'ALL';
%SleepState.ints.ALL = [0 Inf];
statecolors = {'k','b','r',[0.6 0.6 0.6]};

[celltypes,~,typeidx] = unique(CellClass.label);
cellcolor = {'k','r'};

%% The necessary inputs
states = fieldnames(SleepState.ints);
numstates = length(states);
%% 
% loop spike groups, do each spike group in a row
%Get the region of each spike group. Only couple to cells in the region
spikeGroups = sessionInfo.spikeGroups;

%% get the region of the spikegroup
for gg = 1:length(spikeGroups.groups) 
    spikeGroups.region(gg) = unique(sessionInfo.region(ismember(sessionInfo.channels,spikeGroups.groups{gg})));
end



%%
regions = unique(spikeGroups.region);

AllCellClass = cell(1,spikes.numcells);
for rr = 1:length(regions)
    for cc = 1:length(celltypes)
        classname = [regions{rr},'_',celltypes{cc}];
        AllCellClass(CellClass.(celltypes{cc}) & strcmp(spikes.region,regions{rr})) = {classname};
    end
end
allclasses = unique(AllCellClass);
classcolors = {'k','r','k','r'};
classline = {'-','-','--','--'};
%%

rip.fband = [120 200];
%%
for gg = 1:length(spikeGroups.groups)
    display(['Mapping Spike Group ',num2str(gg)])
    %% Load the LFP in this spike group
    downsamplefactor = 2;
    lfp = bz_GetLFP(spikeGroups.groups{gg},...
        'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);

    %% Restrict to state
    for ss = 1:numstates
        state = states{ss};


        %% Spike coupling with wavelets

        %Take only subset of time (random intervals) so wavelets doesn't break
        %computer (total 625s)
        usetime = 2000;%2500
        winsize = 25;
        if sum(diff(SleepState.ints.(state),1,2))>usetime
            nwin = round(usetime./winsize);
            %winsize = 30; %s
            try
                windows = bz_RandomWindowInIntervals( SleepState.ints.(state),winsize,nwin );
            catch
                windows = SleepState.ints.(state);
            end
        else
            windows = SleepState.ints.(state);
        end

        %Calculate pop-phase coupling for all channels
        [SpikeRippleCoupling(gg).(state)] = ...
            bz_GenSpikeLFPCoupling(spikes,lfp,'channel',spikeGroups.groups{gg},...
            'int',windows,'frange',rip.fband,'ncyc',4,...
            'cellclass',AllCellClass,'synchwin',0.002,'synchdt',0.002,...
            'nfreqs',1,'ISIpower',true,'spikeLim',10000);
            close all

    end
    
end
%Pick the channel with the best coupling to the pE population
% [~,gamma.pEchanIDX] = max(synchgammaphasecoupling.pE);
% gamma.pEchan = lfp.channels(gamma.pEchanIDX);
%%
SpikRPCollapsed = bz_CollapseStruct(SpikeRippleCoupling,'match','justcat',true);

%%
thisstate = 'NREMstate';
figure
subplot(2,2,1)
    hold on
    for cc = 1:length(allclasses)
        plot(SpikRPCollapsed.(thisstate).pop.(allclasses{cc}).phasemag,...
            SpikRPCollapsed.(thisstate).detectorinfo.detectionchannel,...
            classline{cc},'color',classcolors{cc})
    end
    legend(allclasses)

subplot(2,3,4)
    imagesc(squeeze(SpikRPCollapsed.NREMstate.cell.ISIpowermodulation)')
    
subplot(2,3,5)
    imagesc(squeeze(SpikRPCollapsed.NREMstate.cell.spikephasemag)')
    
subplot(2,3,6)
    imagesc(squeeze(SpikRPCollapsed.NREMstate.cell.ratepowercorr)')
%% Load the figure for Wavelet coupling map

load('/home/dlevenstein/Desktop/LFPcoupling.mat')

%% Figure: pop-phase coupling by channel
figure
    for tt = 1:length(celltypes)
        subplot(2,3,tt)
        imagesc(log2(freqs),[1 length(ctxchan)],synchphasecoupling.(celltypes{tt})(spikegroupsort,:));
        hold on
        plot(log2([gamma.fband ; gamma.fband]),[1 length(ctxchan);1 length(ctxchan)]',...
            'w--')
        plot(bz_NormToRange(synchgammaphasecoupling.(celltypes{tt})(spikegroupsort),log2(gamma.fband)),...
            1:length(ctxchan),'color',cellcolor{tt})
        plot(log2(gamma.fband(2)),find(ismember(spikegroupordering,gamma.pEchan)),'+')
        title([celltypes{tt},' Pop. Rate'])
        LogScale('x',2)
        colorbar('northoutside')

        xlabel('freq (Hz)');
        if tt==1
        ylabel({'Phase Coupling','Channel'})
        end 
        %xlim(log2([16 192]))
        
        
        subplot(2,3,tt+3)
        imagesc(log2(freqs),[1 length(ctxchan)],synchpowercorr.(celltypes{tt})(spikegroupsort,:));
        hold on
        plot(log2([gamma.fband ; gamma.fband]),[1 length(ctxchan);1 length(ctxchan)]',...
            'w--')
        plot(bz_NormToRange(synchgammapowercorr.(celltypes{tt})(spikegroupsort),log2(gamma.fband)),...
            1:length(ctxchan),'color',cellcolor{tt})
        title([celltypes{tt},' Pop. Rate'])
        LogScale('x',2)
        colorbar('northoutside')
        xlabel('freq (Hz)');
        if tt==1
        ylabel({'Rate-Power Corr.','Channel'})
        end 

    end
    
    
    subplot(2,3,3)
    imagesc([1 spikes.numcells],[1 length(ctxchan)],spikephasemag(spikegroupsort,cellchanneltypesort));
    hold on
    idxvec = 0;
    for tt = 1:length(celltypes)
        idxvec=idxvec(end)+[1:sum(CellClass.(celltypes{tt}))];
         plot(idxvec,cellchannelIDX(CellClass.(celltypes{tt})),...
             '.','color',cellcolor{tt})
    end
    title('All Cells')
            colorbar('northoutside')

   % LogScale('x',2)
    xlabel('Cell');%ylabel('Channel')
    subplot(2,3,6)
    imagesc([1 spikes.numcells],[1 length(ctxchan)],ratepowercorr(spikegroupsort,cellchanneltypesort));
    hold on
    idxvec = 0;
    for tt = 1:length(celltypes)
        idxvec=idxvec(end)+[1:sum(CellClass.(celltypes{tt}))];
         plot(idxvec,cellchannelIDX(CellClass.(celltypes{tt})),...
             '.','color',cellcolor{tt})
    end
    title('All Cells')
            colorbar('northoutside')

   % LogScale('x',2)
    xlabel('Cell');%ylabel('Channel')
NiceSave(['GammaCoupling',state],figfolder,baseName,'includeDate',true)

%%
figure
plot(log10(ISIStats.summstats.(state).meanrate),spikephasemag,'.')
