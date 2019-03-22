function [ ] = Analysis20190322(basePath,figfolder)
% Date 03/20/19
%
%Questions: Ripple coupling and modulation of the ISI distribution - by
%channel. How local is phase coupling and ISI modulation?
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

%%
SpikRPCollapsed = bz_CollapseStruct(SpikeRippleCoupling,'match','justcat',true);

%%
thisstate = 'NREMstate';
figure
subplot(2,2,2)
    hold on
    for cc = 1:length(allclasses)
        plot(SpikRPCollapsed.(thisstate).pop.(allclasses{cc}).phasemag,...
            SpikRPCollapsed.(thisstate).detectorinfo.detectionchannel,...
            classline{cc},'color',classcolors{cc})
    end
    legend(allclasses,'location','eastoutside')
    axis tight
    title('Spike-Phase')

subplot(2,4,1)
    imagesc(squeeze(SpikRPCollapsed.NREMstate.cell.ISIpowermodulation)')
    hold on
    plot([1:spikes.numcells],spikes.maxWaveformCh,'.r','markersize',0.5)
    axis xy
    xlabel('Cell');ylabel('Channel')
    title('Power-ISI')
    
subplot(2,4,2)
    imagesc(squeeze(SpikRPCollapsed.NREMstate.cell.spikephasemag)')
    hold on
    %check indexinghere...
    plot([1:spikes.numcells],spikes.maxWaveformCh,'.r','markersize',0.5)

    axis xy
    xlabel('Cell');ylabel('Channel')
    title('Spike-Phase')

    
subplot(2,3,6)
    imagesc(squeeze(SpikRPCollapsed.NREMstate.cell.ratepowercorr)')
    xlabel('Cell');ylabel('Channel')
    title('Rate-Power')
    
NiceSave(['NREMRippleCouplingByChannel',num2str(gg)],figfolder,baseName,'includeDate',true)
%% Load the figure for Wavelet coupling map

load('/home/dlevenstein/Desktop/LFPcoupling.mat')


%% Ripple coupling by on/off group/region

for gg = 1:length(spikeGroups.groups)
    %Cells in same spike group
    spikeGroups.ingroupcellIDX{gg} = ismember(spikes.maxWaveformCh,spikeGroups.groups{gg});
    spikeGroups.inregioncellIDX{gg} = ismember(spikes.region,spikeGroups.region{gg});
    
    for ss = 1:numstates
        state = states{ss};
        
        %Rate-power
       SpikeRippleCoupling(gg).(state).ingroupmean.ratepowercorr = ...
           squeeze(mean(SpikeRippleCoupling(gg).(state).cell.ratepowercorr(spikeGroups.ingroupcellIDX{gg},:,:),1));
       SpikeRippleCoupling(gg).(state).outgroupmean.ratepowercorr = ...
           squeeze(mean(SpikeRippleCoupling(gg).(state).cell.ratepowercorr(...
           spikeGroups.inregioncellIDX{gg}&~spikeGroups.ingroupcellIDX{gg},:,:),1));
       SpikeRippleCoupling(gg).(state).outregionmean.ratepowercorr = ...
           squeeze(mean(SpikeRippleCoupling(gg).(state).cell.ratepowercorr(...
           ~spikeGroups.inregioncellIDX{gg},:,:),1));

       %Spike-Phase
       SpikeRippleCoupling(gg).(state).ingroupmean.spikephasemag = ...
           squeeze(mean(SpikeRippleCoupling(gg).(state).cell.spikephasemag(spikeGroups.ingroupcellIDX{gg},:,:),1));
       SpikeRippleCoupling(gg).(state).outgroupmean.spikephasemag = ...
           squeeze(mean(SpikeRippleCoupling(gg).(state).cell.spikephasemag(...
           spikeGroups.inregioncellIDX{gg}&~spikeGroups.ingroupcellIDX{gg},:,:),1));
       SpikeRippleCoupling(gg).(state).outregionmean.spikephasemag = ...
           squeeze(mean(SpikeRippleCoupling(gg).(state).cell.spikephasemag(...
           ~spikeGroups.inregioncellIDX{gg},:,:),1));

       %ISI-Power
       SpikeRippleCoupling(gg).(state).ingroupmean.ISIpowermodulation = ...
           squeeze(mean(SpikeRippleCoupling(gg).(state).cell.ISIpowermodulation(spikeGroups.ingroupcellIDX{gg},:,:),1));
       SpikeRippleCoupling(gg).(state).outgroupmean.ISIpowermodulation = ...
           squeeze(mean(SpikeRippleCoupling(gg).(state).cell.ISIpowermodulation(...
           spikeGroups.inregioncellIDX{gg}&~spikeGroups.ingroupcellIDX{gg},:,:),1));
       SpikeRippleCoupling(gg).(state).outregionmean.ISIpowermodulation = ...
           squeeze(mean(SpikeRippleCoupling(gg).(state).cell.ISIpowermodulation(...
           ~spikeGroups.inregioncellIDX{gg},:,:),1)); 

       for tt = 1:length(allclasses)
           SpikeRippleCoupling(gg).(state).pop.(allclasses{tt}).ISIpowermodulation = ...
               squeeze(mean(SpikeRippleCoupling(gg).(state).cell.ISIpowermodulation(...
               ismember(AllCellClass,allclasses{tt}),:,:),1));
       end
    end
end

%%
%Pick the channel with the best coupling to the in-region cells for each
%spikeGroup during NREM
for gg = 1:length(spikeGroups.groups)
    [~,rip.pEchanIDX(gg)] = max(SpikeRippleCoupling(gg).NREMstate.outgroupmean.spikephasemag);
    rip.pEchan(gg) = spikeGroups.groups{gg}(rip.pEchanIDX(gg));
end
sessionInfo = bz_tagChannel(basePath,rip.pEchan,'ripchans');
%% In/Out Group Cells: Ripple Coupling
inout = {'ingroupmean','outgroupmean','outregionmean'};

for gg = 1:length(spikeGroups.groups)

figure
    
for ss = 1:3
        subplot(3,3,ss*3-2)
            hold on
            for io = 1:3
            plot(SpikeRippleCoupling(gg).(states{ss}).(inout{io}).spikephasemag,...
                -[1:length(spikeGroups.groups{gg})]);
            end
            pickchan = spikeGroups.groups{gg}==rip.pEchan(gg);
            plot(SpikeRippleCoupling(gg).(states{ss}).outgroupmean.spikephasemag(pickchan),...
                -find(pickchan),'*')
           % plot(log2(rip.fband(2)),find(ismember(spikegroupordering,rip.pEchan)),'+')
           if ss ==1
            title('Phase Coupling')
           end
            %xlabel('freq (Hz)');
            ylabel({states{ss},'Channel'})
            
            
        subplot(3,3,ss*3-1)
            hold on
            for io = 1:3
            plot(SpikeRippleCoupling(gg).(states{ss}).(inout{io}).ratepowercorr,...
                -[1:length(spikeGroups.groups{gg})]);
            end
           % plot(log2(rip.fband(2)),find(ismember(spikegroupordering,rip.pEchan)),'+')
           if ss ==1
            title('Rate Modulation')
           end
            %xlabel('freq (Hz)');
            %ylabel({states{ss},'Channel'})

        subplot(3,3,ss*3)
            hold on
            for io = 1:3
            plot(SpikeRippleCoupling(gg).(states{ss}).(inout{io}).ISIpowermodulation,...
                -[1:length(spikeGroups.groups{gg})]);
            end
           % plot(log2(rip.fband(2)),find(ismember(spikegroupordering,rip.pEchan)),'+')
           if ss ==1
            title('ISI Modulation')
           end
            %xlabel('freq (Hz)');
            %ylabel({states{ss},'Channel'})
end 


NiceSave(['SpikeRippleCoupling_InOutGroup_SpikeGroup',num2str(gg)],figfolder,baseName,'includeDate',true)

end





%%
for gg = 1:length(spikeGroups.groups)
    
    for ss = 1:numstates
        state = states{ss};
       SpikeLFPCoupling(gg).(state).ingroupmean.ratepowercorr = ...
           squeeze(mean(SpikeLFPCoupling(gg).(state).cell.ratepowercorr(spikeGroups.ingroupcellIDX{gg},:,:),1));
       SpikeLFPCoupling(gg).(state).outgroupmean.ratepowercorr = ...
           squeeze(mean(SpikeLFPCoupling(gg).(state).cell.ratepowercorr(...
           spikeGroups.inregioncellIDX{gg}&~spikeGroups.ingroupcellIDX{gg},:,:),1));
       SpikeLFPCoupling(gg).(state).outregionmean.ratepowercorr = ...
           squeeze(mean(SpikeLFPCoupling(gg).(state).cell.ratepowercorr(...
           ~spikeGroups.inregioncellIDX{gg},:,:),1));   

       SpikeLFPCoupling(gg).(state).ingroupmean.spikephasemag = ...
           squeeze(mean(SpikeLFPCoupling(gg).(state).cell.spikephasemag(spikeGroups.ingroupcellIDX{gg},:,:),1));
       SpikeLFPCoupling(gg).(state).outgroupmean.spikephasemag = ...
           squeeze(mean(SpikeLFPCoupling(gg).(state).cell.spikephasemag(...
           spikeGroups.inregioncellIDX{gg}&~spikeGroups.ingroupcellIDX{gg},:,:),1));
       SpikeLFPCoupling(gg).(state).outregionmean.spikephasemag = ...
           squeeze(mean(SpikeLFPCoupling(gg).(state).cell.spikephasemag(...
           ~spikeGroups.inregioncellIDX{gg},:,:),1));

       SpikeLFPCoupling(gg).(state).ingroupmean.ISIpowermodulation = ...
           squeeze(mean(SpikeLFPCoupling(gg).(state).cell.ISIpowermodulation(spikeGroups.ingroupcellIDX{gg},:,:),1));
       SpikeLFPCoupling(gg).(state).outgroupmean.ISIpowermodulation = ...
           squeeze(mean(SpikeLFPCoupling(gg).(state).cell.ISIpowermodulation(...
           spikeGroups.inregioncellIDX{gg}&~spikeGroups.ingroupcellIDX{gg},:,:),1));
       SpikeLFPCoupling(gg).(state).outregionmean.ISIpowermodulation = ...
           squeeze(mean(SpikeLFPCoupling(gg).(state).cell.ISIpowermodulation(...
           ~spikeGroups.inregioncellIDX{gg},:,:),1));

%        for tt = 1:2
%            SpikeLFPCoupling(gg).(state).pop(tt).ISIpowermodulation = ...
%                squeeze(mean(SpikeLFPCoupling(gg).(state).cell.ISIpowermodulation(CellClass.(celltypes{tt}),:,:),1));
%        end
    end
end




%% In/Out Group Cells

for gg = 1:length(spikeGroups.groups)
figure
for io = 1:3
    
for ss = 1:2
        subplot(6,3,1+(io-1)*3+(ss-1)*9)
        imagesc(log2(SpikeLFPCoupling(gg).(states{ss}).freqs),...
            [1 length(spikeGroups.groups{gg})],...
            squeeze(SpikeLFPCoupling(gg).(states{ss}).(inout{io}).spikephasemag)');
        hold on
        plot(log2([rip.fband ; rip.fband]),[1 length(spikeGroups.groups{gg});1 length(spikeGroups.groups{gg})]',...
            'w--')
        plot(bz_NormToRange(SpikeRippleCoupling(gg).(states{ss}).(inout{io}).spikephasemag,log2(rip.fband)),...
            1:length(spikeGroups.groups{gg}),'color','w')

       % plot(log2(rip.fband(2)),find(ismember(spikegroupordering,rip.pEchan)),'+')
       if io==1 & ss==1
        title('Phase Coupling')
       end
        LogScale('x',2)
        colorbar
        
        %caxis([0 0.12])
        if io ==3
        xlabel('freq (Hz)');
        end
        %if ss==1
        ylabel({states{ss},'Channel'})
       % end 
        %xlim(log2([16 192]))
        
        
        subplot(6,3,2+(io-1)*3+(ss-1)*9)
        imagesc(log2(SpikeLFPCoupling(gg).(states{ss}).freqs),...
            [1 length(spikeGroups.groups{gg})],...
            squeeze(SpikeLFPCoupling(gg).(states{ss}).(inout{io}).ratepowercorr)');
        hold on
        plot(log2([rip.fband ; rip.fband]),[1 length(spikeGroups.groups{gg});1 length(spikeGroups.groups{gg})]',...
            'w--')
        plot(bz_NormToRange(SpikeRippleCoupling(gg).(states{ss}).(inout{io}).ratepowercorr,log2(rip.fband)),...
            1:length(spikeGroups.groups{gg}),'color','w')
        if io==1 & ss==1
        title('Rate-Power')
        end
        LogScale('x',2)
        %caxis([0 0.2])
        colorbar
        if io ==3
        xlabel('freq (Hz)');
        end
        
        subplot(6,3,3+(io-1)*3+(ss-1)*9)
        imagesc(log2(SpikeLFPCoupling(gg).(states{ss}).freqs),...
            [1 length(spikeGroups.groups{gg})],...
            squeeze(SpikeLFPCoupling(gg).(states{ss}).(inout{io}).ISIpowermodulation)');
        hold on
        plot(log2([rip.fband ; rip.fband]),[1 length(spikeGroups.groups{gg});1 length(spikeGroups.groups{gg})]',...
            'w--')
        plot(bz_NormToRange(SpikeRippleCoupling(gg).(states{ss}).(inout{io}).ISIpowermodulation,log2(rip.fband)),...
            1:length(spikeGroups.groups{gg}),'color','w')
        if io==1 & ss==1
        title('Power-ISI')
        end
        LogScale('x',2)
        %caxis([0 0.2])
        colorbar
        if io ==3
        xlabel('freq (Hz)');
        end

end 
end

NiceSave(['SpikeLFPCoupling_InOutGroup_SpikeGroup',num2str(gg)],figfolder,baseName,'includeDate',true)

end












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
