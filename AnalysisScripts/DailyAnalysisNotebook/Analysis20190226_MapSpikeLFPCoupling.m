function [ ] = Analysis20190227(basePath,figfolder)
% Date XX/XX/20XX
%
%Goal: create a function to map pop-LFP coupling for every channel in a
%recording
%
%Plots
%-
%-
%
%% Load Header
%Initiate Paths
%reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/Cicero_09102014';
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
%%

for gg = 8
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
            windows = bz_RandomWindowInIntervals( SleepState.ints.(state),winsize,nwin );
        else
            windows = SleepState.ints.(state);
        end

        %Calculate pop-phase coupling for all channels
        [SpikeLFPCoupling(gg).(state)] = ...
            bz_GenSpikeLFPCoupling(spikes,lfp,'channel',spikeGroups.groups{gg},...
            'int',windows,'frange',[1 312],'ncyc',15,...
            'subpop',CellClass.pE+2.*CellClass.pI,'synchwin',0.002,'synchdt',0.002,...
            'nfreqs',150,'ISIpower',true,'spikeLim',20000);
%             for tt = 1:length(celltypes)
%                 synchphasecoupling.(state).(celltypes{tt})(:,:,gg) = synchcoupling(tt).phasemag;
%                 synchphaseangle.(state).(celltypes{tt})(:,:,gg) = synchcoupling(tt).phaseangle;
%                 synchpowercorr.(state).(celltypes{tt})(:,:,gg) = synchcoupling(tt).powercorr;
%             end
            close all

    end
    
end

%%
%Need: region for each spike group. Split populations by region and cell
%type
%Histogram of peak channel for cells (by population)

%Histogram of proportion of cells coupled to each frequency? (need
%significance....)

%Use this to pick gamma(/ripple?) channel for each cell
for gg = 1:length(spikeGroups.groups)
    spikeGroups.ingroupcells{gg} = ismember(spikes.maxWaveformCh,spikeGroups.groups{gg});
end

%%
%Mean in group cells
%Mean one group away cells
%Mean opposite hemisphere cells
%Mean same hemisphere cells
%Figure: example in group cell, out of group cell
excell = randsample(find(spikeGroups.ingroupcells{11}),1);
cellchan = find(ismember(spikeGroups.groups{8},spikes.maxWaveformCh(excell)));
if isempty(cellchan); cellchan=nan;end
figure
for ss = 1:3
        subplot(3,3,1+(ss-1)*3)
        imagesc(log2(SpikeLFPCoupling(8).(states{ss}).freqs),...
            [1 length(spikeGroups.groups{8})],...
            squeeze(SpikeLFPCoupling(8).(states{ss}).cell.spikephasemag(excell,:,:))');
        hold on
        plot(log2(SpikeLFPCoupling(8).(states{ss}).freqs(end)),cellchan,'r+')
       % plot(log2(rip.fband(2)),find(ismember(spikegroupordering,rip.pEchan)),'+')
        title('Phase Coupling')
        LogScale('x',2)
        colorbar
        
        %caxis([0 0.12])
        xlabel('freq (Hz)');
        if ss==1
        ylabel({states{ss},'Channel'})
        end 
        %xlim(log2([16 192]))
        
        
        subplot(3,3,2+(ss-1)*3)
        imagesc(log2(SpikeLFPCoupling(8).(states{ss}).freqs),...
            [1 length(spikeGroups.groups{8})],...
            squeeze(SpikeLFPCoupling(8).(states{ss}).cell.ratepowercorr(excell,:,:))');
        hold on
        plot(log2(SpikeLFPCoupling(8).(states{ss}).freqs(end)),cellchan,'r+')
        title('Rate-Power')
        LogScale('x',2)
        %caxis([0 0.2])
        colorbar
        xlabel('freq (Hz)');
        
        subplot(3,3,3+(ss-1)*3)
        imagesc(log2(SpikeLFPCoupling(8).(states{ss}).freqs),...
            [1 length(spikeGroups.groups{8})],...
            squeeze(SpikeLFPCoupling(8).(states{ss}).cell.ISIpowermodulation(:,excell,:))');
        hold on
        plot(log2(SpikeLFPCoupling(8).(states{ss}).freqs(end)),cellchan,'r+')
        title('Power-ISI')
        LogScale('x',2)
        %caxis([0 0.2])
        colorbar
        xlabel('freq (Hz)');

end    

%% Figure: pop-phase coupling by channel
for gg =8

figure
for ss = 1:3
    for tt = 1:length(celltypes)
        subplot(3,4,tt+(ss-1)*4)
        imagesc(log2(SpikeLFPCoupling(gg).(states{ss}).freqs),...
            [1 length(spikeGroups.groups{gg})],...
            SpikeLFPCoupling(gg).(states{ss}).pop(tt).phasemag');
        hold on
       % plot(log2(rip.fband(2)),find(ismember(spikegroupordering,rip.pEchan)),'+')
       % title([celltypes{tt},' Pop. Rate'])
        LogScale('x',2)
        colorbar('northoutside')
        %caxis([0 0.12])
        xlabel('freq (Hz)');
        if tt==1
        ylabel({states{ss},'Channel'})
        end 
        %xlim(log2([16 192]))
        
        
        subplot(3,4,tt+2+(ss-1)*4)
        imagesc(log2(SpikeLFPCoupling(gg).(states{ss}).freqs),...
            [1 length(spikeGroups.groups{gg})],...
            SpikeLFPCoupling(gg).(states{ss}).pop(tt).powercorr');
        hold on
       % title([celltypes{tt},' Pop. Rate'])
        LogScale('x',2)
        %caxis([0 0.2])
        colorbar('northoutside')
        xlabel('freq (Hz)');
        if tt==1
        %ylabel({'Rate-Power Corr.','Channel'})
        end 

    end
end    
    %NiceSave(['SpikeLFPCoupling_SpikeGroup',num2str(gg)],figfolder,baseName,'includeDate',true)


end
end
