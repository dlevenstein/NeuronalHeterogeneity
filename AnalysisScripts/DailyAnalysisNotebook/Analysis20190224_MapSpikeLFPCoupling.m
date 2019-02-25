function [ ] = Analysis20190224(basePath,figfolder)
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
clear synchphasecoupling
clear synchpowercorr
clear synchphaseangle
for gg = 1:2
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
        usetime = 1500;%2500
        winsize = 25;
        if sum(diff(SleepState.ints.(state),1,2))>usetime
            nwin = round(usetime./winsize);
            %winsize = 30; %s
            windows = bz_RandomWindowInIntervals( SleepState.ints.(state),winsize,nwin );
        else
            windows = SleepState.ints.(state);
        end

        %Calculate pop-phase coupling for all channels
        [freqs,synchcoupling] = ...
            bz_GenSpikeLFPCoupling(spikes.times,lfp,'channel',spikeGroups.groups{gg}(cc),...
            'int',windows,'DOWNSAMPLE',1,'frange',[1 312],'ncyc',15,...
            'subpop',CellClass.pE+2.*CellClass.pI,'synchwin',0.002,'synchdt',0.002,...
            'nfreqs',150);
            for tt = 1:length(celltypes)
                synchphasecoupling.(state).(celltypes{tt})(:,:,gg) = synchcoupling(tt).phasemag;
                synchphaseangle.(state).(celltypes{tt})(:,:,gg) = synchcoupling(tt).phaseangle;
                synchpowercorr.(state).(celltypes{tt})(:,:,gg) = synchcoupling(tt).powercorr;
            end
            close all

    end
    
end

%%%% PUT THE CAPABILITY FOR GENSPIKELFPCOUPLING TO DO MULTIPLE CHANNELS...
%% Figure: pop-phase coupling by channel
for gg =1:2;

figure
for ss = 1:3
    for tt = 1:length(celltypes)
        subplot(3,4,tt+(ss-1)*4)
        imagesc(log2(freqs),[1 length(spikeGroups.groups{gg})],synchphasecoupling.(states{ss}).(celltypes{tt})(:,:,gg));
        hold on
       % plot(log2(rip.fband(2)),find(ismember(spikegroupordering,rip.pEchan)),'+')
       % title([celltypes{tt},' Pop. Rate'])
        LogScale('x',2)
        colorbar('northoutside')
        caxis([0 0.12])
        xlabel('freq (Hz)');
        if tt==1
        ylabel({states{ss},'Channel'})
        end 
        %xlim(log2([16 192]))
        
        
        subplot(3,4,tt+2+(ss-1)*4)
        imagesc(log2(freqs),[1 length(spikeGroups.groups{gg})],synchpowercorr.(states{ss}).(celltypes{tt})(:,:,gg));
        hold on
       % title([celltypes{tt},' Pop. Rate'])
        LogScale('x',2)
        caxis([0 0.2])
        colorbar('northoutside')
        xlabel('freq (Hz)');
        if tt==1
        %ylabel({'Rate-Power Corr.','Channel'})
        end 

    end
end    
    NiceSave(['SpikeLFPCoupling_SpikeGroup',num2str(gg)],figfolder,baseName,'includeDate',true)


end
end
