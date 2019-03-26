function [ ] = Analysis20190325(basePath,figfolder)
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
        usetime = 3000;%2500
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
        [SpikeLFPCoupling(gg).(state)] = ...
            bz_GenSpikeLFPCoupling(spikes,lfp,'channel',spikeGroups.groups{gg},...
            'int',windows,'frange',[1 312],'ncyc',15,...
            'cellclass',AllCellClass,'synchwin',0.002,'synchdt',0.002,...
            'nfreqs',150,'ISIpower',true,'spikeLim',20000);
            close all

    end
    
end

%% Saving...
saveswitch = false;

savename = fullfile(basePath,[baseName,'.SpikeLFPCoupling.chaninfo.mat']);

if saveswitch
    save(savename,'SpikeLFPCoupling');
end

%%
%Need: region for each spike group. Split populations by region and cell
%type
%Histogram of peak channel for cells (by population)

%Histogram of proportion of cells coupled to each frequency? (need
%significance....)

%Use this to pick gamma(/ripple?) channel for each cell
for gg = 1:length(spikeGroups.groups)
	spikeGroups.ingroupcellIDX{gg} = ismember(spikes.maxWaveformCh,spikeGroups.groups{gg});
    spikeGroups.inregioncellIDX{gg} = ismember(spikes.region,spikeGroups.region{gg});    
    
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

       for tt = 1:2
           SpikeLFPCoupling(gg).(state).pop.(celltypes{tt}).ISIpowermodulation = ...
               squeeze(mean(SpikeLFPCoupling(gg).(state).cell.ISIpowermodulation(CellClass.(celltypes{tt}),:,:),1));
       end
       
       for tt = 1:length(allclasses)
           SpikeLFPCoupling(gg).(state).pop.(allclasses{tt}).ISIpowermodulation = ...
               squeeze(mean(SpikeLFPCoupling(gg).(state).cell.ISIpowermodulation(...
               ismember(AllCellClass,allclasses{tt}),:,:),1));
       end
    end
end


%% Try Collapsing

%SpikLFPCollapsed = bz_CollapseStruct(SpikeLFPCoupling,'match','justcat',true);


%% R vs L

order = [1 3 2 4];
for gg = 1:length(spikeGroups.groups)
figure
suptitle(['SpikeGroup ',num2str(gg),': ',spikeGroups.region{gg}])
for ss = 1:length(states)
for tt = 1:length(allclasses)
    subplot(4,3,ss+(order(tt)*3)-3)
        imagesc(log2(SpikeLFPCoupling(gg).(states{ss}).freqs),[1 length(spikeGroups.groups{gg})],...
            SpikeLFPCoupling(gg).(states{ss}).pop.(allclasses{tt}).ISIpowermodulation')
        if ss == 1
            ylabel(allclasses{tt})
        end
        LogScale('x',2)
        if tt ==4
            xlabel('f (Hz)')
        end
        colorbar
        if tt==1
            title(states{ss})
        end
end
end

NiceSave(['ISIModulation_LR_',num2str(gg)],figfolder,baseName,'includeDate',true)


end

%%
%state = 'WAKEstate';
for gg = 1:length(spikeGroups.groups)
figure

specchan = find(ismember(spikeGroups.groups{gg},lamina.RAD.chans) | ismember(spikeGroups.groups{gg},lamina.PYR.chans));

for ss = 1:length(states)
for tt = 1:length(CellClass.celltypes)
    subplot(3,2,tt+(ss*2)-2)
        imagesc(log2(SpikeLFPCoupling(gg).(states{ss}).freqs),[1 length(spikeGroups.groups{gg})],...
            SpikeLFPCoupling(gg).(states{ss}).pop.(celltypes{tt}).ISIpowermodulation')
        hold on
        plot(zeros(size(specchan)),specchan,'r*')
        if ss == 1
            title(celltypes{tt})
        end
        LogScale('x',2)
        switch tt
            case 1
                if ss == 2
                    caxis([0.0125 0.02])
                else
                    caxis([0.0125 0.045])
                end
            case 2
                caxis([0.005 0.02])
        end
        colorbar
        if tt==1
            ylabel(states{ss})
        end
        set(gca,'ytick',[])
        if ss ==3
            xlabel('f (Hz)')
        end
        
            
end
end

NiceSave(['ISIModulation_',num2str(gg)],figfolder,baseName,'includeDate',true)

end

%% Tag LM channel, PYR channel - two for each hemisphere.
lamina.RAD.chans = [19 39 103 113];
lamina.PYR.chans = [10 42 77 116];
[ sessionInfo ] = bz_tagChannel( basePath,[],'LMChan','overwrite',true );
[ sessionInfo ] = bz_tagChannel( basePath,lamina.RAD.chans,'RADChan','overwrite',true );
[ sessionInfo ] = bz_tagChannel( basePath,lamina.PYR.chans,'PYRChan','overwrite',true  );

%% Calcluate conditional spike phase coupling on LM chan and PYR chans

ISIStats.allspikes.logISIs = cellfun(@(X) log10(X),ISIStats.allspikes.ISIs,'UniformOutput',false);
ISIStats.allspikes.logISIs_next = cellfun(@(X) log10([X(2:end);nan]),ISIStats.allspikes.ISIs,'UniformOutput',false);

doubleISIs.times = cellfun(@(X) [X;X],ISIStats.allspikes.times,'UniformOutput',false);
doubleISIs.ISIs = cellfun(@(X,Y) [X;Y],ISIStats.allspikes.logISIs,ISIStats.allspikes.logISIs_next,'UniformOutput',false);

%% Load the LFP for spike-phase

downsamplefactor = 2;
lfp = bz_GetLFP([sessionInfo.channelTags.RADChan sessionInfo.channelTags.PYRChan],...
    'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);

%%
 thchan = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.THchanID;
 downsamplefactor = 5;
 th_lfp = bz_GetLFP(thchan,...
     'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);


thetalfp = bz_Filter(th_lfp,'passband',[6 10]);
thetalfp.amp = NormToInt((thetalfp.amp),'mean',SleepState.ints.WAKEstate);

ISIStats.allspikes.thetapower =cellfun(@(X) ...
    interp1(thetalfp.timestamps,thetalfp.amp,X,'nearest'),...
    ISIStats.allspikes.times,'UniformOutput',false);
%% Restrict to state

for ss = 1:3
    state = states{ss};
    %ints = SleepState.ints.(state);

    %Take only subset of time (random intervals) so wavelets doesn't break
    %computer (total 625s)
    usetime = 3000;%2500
    winsize = 25;
    if sum(diff(SleepState.ints.(state),1,2))>usetime
        nwin = round(usetime./winsize);
        %winsize = 30; %s
        windows = bz_RandomWindowInIntervals( SleepState.ints.(state),winsize,nwin );
    else
        windows = SleepState.ints.(state);
    end


    %% Get complex valued wavelet transform at each timestamp
    for cc = 1:length(lfp.channels)
    wavespec = bz_WaveSpec(lfp,'intervals',windows,'showprogress',true,'ncyc',15,...
        'nfreqs',150,'frange',[1 312],'chanID',lfp.channels(cc));


    %% Coupling COnditioned on ISI
    LFPCoupling(ss,cc) = bz_ConditionalLFPCoupling( doubleISIs,doubleISIs.ISIs,wavespec,...
        'Xbounds',[-2.6 1],'intervals',windows,'showFig',true,... %true
    'minX',25,'CellClass',CellClass,'spikeLim',20000,...
    'showFig',true);

    %% Coupling Conditioned on theta
        if ss~=2
            LFPCoupling_theta(ss,cc) = bz_ConditionalLFPCoupling( ISIStats.allspikes,ISIStats.allspikes.thetapower,wavespec,...
                'Xbounds',[0.1 2.5],'intervals',windows,'showFig',true,'numXbins',30,...
            'minX',25,'CellClass',CellClass,'spikeLim',20000,...
            'showFig',true);
        end
    end

end

%% For each cell - pick its LM and PYR channel. 
%Combine LFPCoupling results for each cell from the appropriate channels
lamina.names = {'PYR','RAD'};

for ll = 1:length(lamina.names)
    lamina.(lamina.names{ll}).spikegroup = zeros(size(lamina.(lamina.names{ll}).chans));

    for gg = 1:length(spikeGroups.groups)
        lamina.(lamina.names{ll}).spikegroup = ...
            lamina.(lamina.names{ll}).spikegroup+gg.*ismember(lamina.(lamina.names{ll}).chans,spikeGroups.groups{gg});
    end

    lamina.(lamina.names{ll}).chanIDX = find(ismember(sessionInfo.channels,lamina.(lamina.names{ll}).chans));

    lamina.(lamina.names{ll}).region = sessionInfo.region(lamina.(lamina.names{ll}).chanIDX);
end
%%
for cc = 1:spikes.numcells
    for ll = 1:length(lamina.names)
    %In same hemisphere, not in same spikegroup
    sameregion = strcmp(spikes.region(cc),lamina.(lamina.names{ll}).region);
    samespikegroup = spikes.shankID(cc) == lamina.(lamina.names{ll}).spikegroup;

    lamina.(lamina.names{ll}).cellchan(cc) = randsample(find(sameregion&~samespikegroup),1);
    end
end

%% Get the spike-lfp coupling for each cell from its appropriate spot
clear laminarLFPCoupling
ss = 1;
for ll =1:2
    %ll = 2;
    offset = 0;
    if ll == 1
        offset = length(lamina.RAD.chans);
    end
    fields = fieldnames(LFPCoupling);
    laminarLFPCoupling.(lamina.names{ll}).Xbins = LFPCoupling(1).Xbins;
    laminarLFPCoupling.(lamina.names{ll}).freqs	 = LFPCoupling(1).freqs;
    laminarLFPCoupling_theta.(lamina.names{ll}).Xbins = LFPCoupling_theta(1).Xbins;
    laminarLFPCoupling_theta.(lamina.names{ll}).freqs = LFPCoupling_theta(1).freqs;

    for cc = 1:spikes.numcells
    laminarLFPCoupling.(lamina.names{ll}).mrl(:,:,cc) = ...
        LFPCoupling(ss,lamina.(lamina.names{ll}).cellchan(cc)+offset).mrl(:,:,cc);
    laminarLFPCoupling.(lamina.names{ll}).meanpower(:,:,cc) = ...
        LFPCoupling(ss,lamina.(lamina.names{ll}).cellchan(cc)+offset).meanpower(:,:,cc);
    laminarLFPCoupling.(lamina.names{ll}).mutInfoXPower(cc,:) = ...
        LFPCoupling(ss,lamina.(lamina.names{ll}).cellchan(cc)+offset).mutInfoXPower(cc,:);
    
    laminarLFPCoupling_theta.(lamina.names{ll}).mrl(:,:,cc) = ...
        LFPCoupling_theta(ss,lamina.(lamina.names{ll}).cellchan(cc)+offset).mrl(:,:,cc);
    laminarLFPCoupling_theta.(lamina.names{ll}).meanpower(:,:,cc) = ...
        LFPCoupling_theta(ss,lamina.(lamina.names{ll}).cellchan(cc)+offset).meanpower(:,:,cc);
    laminarLFPCoupling_theta.(lamina.names{ll}).mutInfoXPower(cc,:) = ...
        LFPCoupling_theta(ss,lamina.(lamina.names{ll}).cellchan(cc)+offset).mutInfoXPower(cc,:);
    end


    %Get Cell Type Average
    for tt = 1:length(celltypes)
        laminarLFPCoupling.(lamina.names{ll}).allmeanpower.(celltypes{tt}) = ...
            nanmean(laminarLFPCoupling.(lamina.names{ll}).meanpower(:,:,CellClass.(celltypes{tt})),3);
        laminarLFPCoupling.(lamina.names{ll}).almeanpMRL.(celltypes{tt}) = ...
            nanmean(laminarLFPCoupling.(lamina.names{ll}).mrl(:,:,CellClass.(celltypes{tt})),3);
        laminarLFPCoupling.(lamina.names{ll}).groupmutinf.(celltypes{tt}) = ...
            nanmean(laminarLFPCoupling.(lamina.names{ll}).mutInfoXPower(CellClass.(celltypes{tt}),:),1);

        laminarLFPCoupling_theta.(lamina.names{ll}).allmeanpower.(celltypes{tt}) = ...
            nanmean(laminarLFPCoupling_theta.(lamina.names{ll}).meanpower(:,:,CellClass.(celltypes{tt})),3);
        laminarLFPCoupling_theta.(lamina.names{ll}).almeanpMRL.(celltypes{tt}) = ...
            nanmean(laminarLFPCoupling_theta.(lamina.names{ll}).mrl(:,:,CellClass.(celltypes{tt})),3);
        laminarLFPCoupling_theta.(lamina.names{ll}).groupmutinf.(celltypes{tt}) = ...
            nanmean(laminarLFPCoupling_theta.(lamina.names{ll}).mutInfoXPower(CellClass.(celltypes{tt}),:),1);
    end

end

%%

    powermap = makeColorMap([0 0 0.8],[1 1 1],[0.8 0 0]);
    figure
    
    subplot(3,3,1)
        hold on
        for tt = 1:length(celltypes)
            plot(log2(laminarLFPCoupling.PYR.freqs),laminarLFPCoupling.PYR.groupmutinf.(celltypes{tt}),...
                'linewidth',1,'color',cellcolor{tt})
        end
        for tt = 1:length(celltypes)
            plot(log2(laminarLFPCoupling.RAD.freqs),laminarLFPCoupling.RAD.groupmutinf.(celltypes{tt}),...
                '--','linewidth',1,'color',cellcolor{tt})
        end
        box off
        axis tight
        xlabel('f (Hz)');ylabel('I(Power;ISI)')
            LogScale('x',2)    
            
    
    for ll=1:2
    for tt = 1:length(celltypes)
    subplot(4,3,tt+(ll*3)+4)
    colormap(gca,powermap)
        imagesc(laminarLFPCoupling.(lamina.names{ll}).Xbins,log2(laminarLFPCoupling.(lamina.names{ll}).freqs),...
            log2(laminarLFPCoupling.(lamina.names{ll}).allmeanpower.(celltypes{tt}))')
        colorbar
        ColorbarWithAxis([-1.25 1.25],'Power (mean^-^1)')
        LogScale('x',10);
        LogScale('y',2)
        LogScale('c',2)
        axis xy
        xlabel('ISI (s)');ylabel('freq (Hz)')
        title((celltypes{tt}))
    end  
    
    
    for tt = 1:length(celltypes)
    subplot(4,3,tt+(ll*3)-2)
        imagesc(laminarLFPCoupling.(lamina.names{ll}).Xbins,log2(laminarLFPCoupling.(lamina.names{ll}).freqs),...
            laminarLFPCoupling.(lamina.names{ll}).almeanpMRL.(celltypes{tt})')
        colorbar
        hold on
        %caxis([0.5 1.5])
        LogScale('x',10);
        LogScale('y',2)
        ColorbarWithAxis([0 0.6],'Phase Coupling (pMRL)')

        axis xy
        xlabel('ISI (s)');ylabel('freq (Hz)')
        title((celltypes{tt}))
    end 
    end
    
    NiceSave('CouplingbyISI',figfolder,baseName,'includeDate',true)

    

%% Figure Theta
    figure
    
    subplot(3,3,1)
        hold on
        for tt = 1:length(celltypes)
            plot(log2(laminarLFPCoupling_theta.PYR.freqs),laminarLFPCoupling_theta.PYR.groupmutinf.(celltypes{tt}),...
                'linewidth',1,'color',cellcolor{tt})
        end
        for tt = 1:length(celltypes)
            plot(log2(laminarLFPCoupling_theta.RAD.freqs),laminarLFPCoupling_theta.RAD.groupmutinf.(celltypes{tt}),...
                '--','linewidth',1,'color',cellcolor{tt})
        end
        box off
        axis tight
        xlabel('f (Hz)');ylabel('I(Power;Theta)')
            LogScale('x',2)    
            
    
    for ll=1:2
    for tt = 1:length(celltypes)
    subplot(4,3,tt+(ll*3)+4)
    colormap(gca,powermap)
        imagesc(laminarLFPCoupling_theta.(lamina.names{ll}).Xbins,log2(laminarLFPCoupling_theta.(lamina.names{ll}).freqs),...
            log2(laminarLFPCoupling_theta.(lamina.names{ll}).allmeanpower.(celltypes{tt}))')
        colorbar
        ColorbarWithAxis([-1.25 1.25],'Power (mean^-^1)')
        %LogScale('x',10);
        LogScale('y',2)
        LogScale('c',2)
        axis xy
        xlabel('Theta Power (mean^-1)');ylabel('freq (Hz)')
        title((celltypes{tt}))
    end  
    
    
    for tt = 1:length(celltypes)
    subplot(4,3,tt+(ll*3)-2)
        imagesc(laminarLFPCoupling_theta.(lamina.names{ll}).Xbins,log2(laminarLFPCoupling_theta.(lamina.names{ll}).freqs),...
            laminarLFPCoupling_theta.(lamina.names{ll}).almeanpMRL.(celltypes{tt})')
        colorbar
        hold on
        %caxis([0.5 1.5])
        %LogScale('x',10);
        LogScale('y',2)
        ColorbarWithAxis([0 0.5],'Phase Coupling (pMRL)')

        axis xy
        xlabel('Theta Power (mean^-1)');ylabel('freq (Hz)')
        title((celltypes{tt}))
    end 
    end
    
    NiceSave('CouplingbyTheta',figfolder,baseName,'includeDate',true)

end
