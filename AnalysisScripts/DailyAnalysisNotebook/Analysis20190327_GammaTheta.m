function [ ] = Analysis20190327(basePath,figfolder)
% Date 03/27/2019
%
%Goal: invesitgate the relationship between theta and high gamma/ripple
%oscillation, and specifically how it effects the ISI distribution
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
spikeGroups = sessionInfo.spikeGroups;
% get the region of the spikegroup
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
%% Get theta power at each time/spike

 thchan = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.THchanID;
 downsamplefactor = 5;
 th_lfp = bz_GetLFP(thchan,...
     'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);

%%
thetalfp = bz_Filter(th_lfp,'passband',[6 10]);
deltalfp = bz_Filter(th_lfp,'passband',[2 20]);

thetalfp.thetadelta = NormToInt(thetalfp.amp./deltalfp.amp,'mean',SleepState.ints.WAKEstate);
thetalfp.amp = NormToInt((thetalfp.amp),'mean',SleepState.ints.WAKEstate);


ISIStats.allspikes.thetapower =cellfun(@(X) ...
    interp1(thetalfp.timestamps,thetalfp.amp,X,'nearest'),...
    ISIStats.allspikes.times,'UniformOutput',false);
ISIStats.allspikes.thetarat =cellfun(@(X) ...
    interp1(thetalfp.timestamps,thetalfp.thetadelta,X,'nearest'),...
    ISIStats.allspikes.times,'UniformOutput',false);
ISIStats.allspikes.thetarat_log = cellfun(@(X) log2(X),ISIStats.allspikes.thetarat,'UniformOutput',false);

%% Get gamma power at each time(/spike)

 gammachan = sessionInfo.channelTags.PYRChan;
 downsamplefactor = 2;
 ga_lfp = bz_GetLFP(gammachan,...
     'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);


gammalfp = bz_Filter(ga_lfp,'passband',[120 200]);
gammalfp.amp = NormToInt((gammalfp.amp),'mean',SleepState.ints.WAKEstate);

%% Calculate Joint Probability

thetalfp.gammaamp = interp1(gammalfp.timestamps,gammalfp.amp,thetalfp.timestamps,'nearest');

%%
thetalfp.wakeidx = InIntervals(thetalfp.timestamps,SleepState.ints.WAKEstate);
%%
[ga_bytheta] = ConditionalHist(log2(thetalfp.amp(thetalfp.wakeidx)),log10(thetalfp.gammaamp(thetalfp.wakeidx)),...
    'Xbounds',[-3 2],'numXbins',30,'Ybounds',[-2 1.5],'numYbins',150,'minX',400);
[ga_bythetarat] = ConditionalHist(log2(thetalfp.thetadelta(thetalfp.wakeidx)),log10(thetalfp.gammaamp(thetalfp.wakeidx)),...
    'Xbounds',[-3 2.5],'numXbins',30,'Ybounds',[-2 1.5],'numYbins',150,'minX',400);

%%
figure
subplot(2,2,1)
imagesc(ga_bytheta.Xbins,ga_bytheta.Ybins,log10(ga_bytheta.pYX)')
axis xy
LogScale('y',10)
LogScale('x',2)

subplot(2,2,3)
bar(ga_bytheta.Xbins,ga_bytheta.Xhist)
axis tight
box off

subplot(2,2,2)
imagesc(ga_bythetarat.Xbins,ga_bythetarat.Ybins,log10(ga_bythetarat.pYX)')
axis xy
LogScale('y',10)
LogScale('x',2)

subplot(2,2,4)
bar(ga_bythetarat.Xbins,ga_bythetarat.Xhist)
axis tight
box off

    NiceSave('GammaPowerByThetaRat',figfolder,baseName,'includeDate',true)

%%

%% ISI Coupling conditioned on theta, theta/delta

state = states{1};
%ints = SleepState.ints.(state);

%Take only subset of time (random intervals) so wavelets doesn't break
%computer (total 625s)
usetime = 4000;%2500
winsize = 25;
if sum(diff(SleepState.ints.(state),1,2))>usetime
    nwin = round(usetime./winsize);
    %winsize = 30; %s
    windows = bz_RandomWindowInIntervals( SleepState.ints.(state),winsize,nwin );
else
    windows = SleepState.ints.(state);
end

downsamplefactor = 2;
lfp_laminar = bz_GetLFP([sessionInfo.channelTags.PYRChan sessionInfo.channelTags.RADChan],...
    'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);

%% Get complex valued wavelet transform at each timestamp
for cc = 1:length(lfp_laminar.channels)
wavespec = bz_WaveSpec(lfp_laminar,'intervals',windows,'showprogress',true,'ncyc',15,...
    'nfreqs',150,'frange',[1 312],'chanID',lfp_laminar.channels(cc));

%% Coupling Conditioned on theta
    LFPCoupling_theta(cc) = bz_ConditionalLFPCoupling( ISIStats.allspikes,ISIStats.allspikes.thetapower,wavespec,...
        'Xbounds',[0.1 2.5],'intervals',windows,'showFig',true,'numXbins',30,...
    'minX',25,'CellClass',CellClass,'spikeLim',20000,...
    'showFig',true);
%%
    LFPCoupling_thetarat(cc) = bz_ConditionalLFPCoupling( ISIStats.allspikes,ISIStats.allspikes.thetarat,wavespec,...
        'Xbounds',[0.1 2],'intervals',windows,'showFig',true,'numXbins',30,...
    'minX',25,'CellClass',CellClass,'spikeLim',20000,...
    'showFig',true);

end


%% For each cell - pick its LM and PYR channel. 
%Combine LFPCoupling results for each cell from the appropriate channels
lamina.names = {'PYR','RAD'};

for ll = 1:length(lamina.names)
    lamina.(lamina.names{ll}).chans = sessionInfo.channelTags.([(lamina.names{ll}),'Chan']);
    lamina.(lamina.names{ll}).spikegroup = zeros(size(lamina.(lamina.names{ll}).chans));

    for gg = 1:length(spikeGroups.groups)
        lamina.(lamina.names{ll}).spikegroup = ...
            lamina.(lamina.names{ll}).spikegroup+gg.*ismember(lamina.(lamina.names{ll}).chans,spikeGroups.groups{gg});
    end

    lamina.(lamina.names{ll}).chanIDX = find(ismember(sessionInfo.channels,lamina.(lamina.names{ll}).chans));

    lamina.(lamina.names{ll}).region = sessionInfo.region(lamina.(lamina.names{ll}).chanIDX);
end

for cc = 1:spikes.numcells
    for ll = 1:length(lamina.names)
    %In same hemisphere, not in same spikegroup
    sameregion = strcmp(spikes.region(cc),lamina.(lamina.names{ll}).region);
    samespikegroup = spikes.shankID(cc) == lamina.(lamina.names{ll}).spikegroup;

    lamina.(lamina.names{ll}).cellchan(cc) = randsample(find(sameregion&~samespikegroup),1);
    end
end


%% Get the spike-lfp coupling for each cell from its appropriate spot
clear laminarLFPCoupling_theta
clear laminarLFPCoupling_thetarat
ss = 1;
for ll =1:2
    %ll = 2;
    offset = 0;
    if ll == 2
        offset = length(lamina.PYR.chans);
    end
    %fields = fieldnames(LFPCoupling);
    laminarLFPCoupling_thetarat.(lamina.names{ll}).Xbins = LFPCoupling_thetarat(1).Xbins;
    laminarLFPCoupling_thetarat.(lamina.names{ll}).freqs	 = LFPCoupling_thetarat(1).freqs;
    laminarLFPCoupling_theta.(lamina.names{ll}).Xbins = LFPCoupling_theta(1).Xbins;
    laminarLFPCoupling_theta.(lamina.names{ll}).freqs = LFPCoupling_theta(1).freqs;

    for cc = 1:spikes.numcells
    laminarLFPCoupling_thetarat.(lamina.names{ll}).mrl(:,:,cc) = ...
        LFPCoupling_thetarat(ss,lamina.(lamina.names{ll}).cellchan(cc)+offset).mrl(:,:,cc);
    laminarLFPCoupling_thetarat.(lamina.names{ll}).meanpower(:,:,cc) = ...
        LFPCoupling_thetarat(ss,lamina.(lamina.names{ll}).cellchan(cc)+offset).meanpower(:,:,cc);
    laminarLFPCoupling_thetarat.(lamina.names{ll}).mutInfoXPower(cc,:) = ...
        LFPCoupling_thetarat(ss,lamina.(lamina.names{ll}).cellchan(cc)+offset).mutInfoXPower(cc,:);
    
    laminarLFPCoupling_theta.(lamina.names{ll}).mrl(:,:,cc) = ...
        LFPCoupling_theta(ss,lamina.(lamina.names{ll}).cellchan(cc)+offset).mrl(:,:,cc);
    laminarLFPCoupling_theta.(lamina.names{ll}).meanpower(:,:,cc) = ...
        LFPCoupling_theta(ss,lamina.(lamina.names{ll}).cellchan(cc)+offset).meanpower(:,:,cc);
    laminarLFPCoupling_theta.(lamina.names{ll}).mutInfoXPower(cc,:) = ...
        LFPCoupling_theta(ss,lamina.(lamina.names{ll}).cellchan(cc)+offset).mutInfoXPower(cc,:);
    end


    %Get Cell Type Average
    for tt = 1:length(celltypes)
        laminarLFPCoupling_thetarat.(lamina.names{ll}).allmeanpower.(celltypes{tt}) = ...
            nanmean(laminarLFPCoupling_thetarat.(lamina.names{ll}).meanpower(:,:,CellClass.(celltypes{tt})),3);
        laminarLFPCoupling_thetarat.(lamina.names{ll}).almeanpMRL.(celltypes{tt}) = ...
            nanmean(laminarLFPCoupling_thetarat.(lamina.names{ll}).mrl(:,:,CellClass.(celltypes{tt})),3);
        laminarLFPCoupling_thetarat.(lamina.names{ll}).groupmutinf.(celltypes{tt}) = ...
            nanmean(laminarLFPCoupling_thetarat.(lamina.names{ll}).mutInfoXPower(CellClass.(celltypes{tt}),:),1);

        laminarLFPCoupling_theta.(lamina.names{ll}).allmeanpower.(celltypes{tt}) = ...
            nanmean(laminarLFPCoupling_theta.(lamina.names{ll}).meanpower(:,:,CellClass.(celltypes{tt})),3);
        laminarLFPCoupling_theta.(lamina.names{ll}).almeanpMRL.(celltypes{tt}) = ...
            nanmean(laminarLFPCoupling_theta.(lamina.names{ll}).mrl(:,:,CellClass.(celltypes{tt})),3);
        laminarLFPCoupling_theta.(lamina.names{ll}).groupmutinf.(celltypes{tt}) = ...
            nanmean(laminarLFPCoupling_theta.(lamina.names{ll}).mutInfoXPower(CellClass.(celltypes{tt}),:),1);
    end

end


%% Figure Theta

powermap = makeColorMap([0 0 0.8],[1 1 1],[0.8 0 0]);
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
        %title((celltypes{tt}))
                if ll==1
        title((celltypes{tt}))
        end
    end 
    end
    
    NiceSave('CouplingbyTheta',figfolder,baseName,'includeDate',true)
    
    
%% Figure Theta Ratio
    figure
    
    subplot(3,3,1)
        hold on
        for tt = 1:length(celltypes)
            plot(log2(laminarLFPCoupling_thetarat.PYR.freqs),laminarLFPCoupling_thetarat.PYR.groupmutinf.(celltypes{tt}),...
                'linewidth',1,'color',cellcolor{tt})
        end
        for tt = 1:length(celltypes)
            plot(log2(laminarLFPCoupling_thetarat.RAD.freqs),laminarLFPCoupling_thetarat.RAD.groupmutinf.(celltypes{tt}),...
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
        imagesc(laminarLFPCoupling_thetarat.(lamina.names{ll}).Xbins,log2(laminarLFPCoupling_thetarat.(lamina.names{ll}).freqs),...
            log2(laminarLFPCoupling_thetarat.(lamina.names{ll}).allmeanpower.(celltypes{tt}))')
        colorbar
        ColorbarWithAxis([-1.25 1.25],'Power (mean^-^1)')
        %LogScale('x',10);
        LogScale('y',2)
        LogScale('c',2)
        axis xy
        xlabel('Theta Power (mean^-1)');ylabel('freq (Hz)')
    end  
    
    
    for tt = 1:length(celltypes)
    subplot(4,3,tt+(ll*3)-2)
        imagesc(laminarLFPCoupling_thetarat.(lamina.names{ll}).Xbins,log2(laminarLFPCoupling_thetarat.(lamina.names{ll}).freqs),...
            laminarLFPCoupling_thetarat.(lamina.names{ll}).almeanpMRL.(celltypes{tt})')
        colorbar
        hold on
        %caxis([0.5 1.5])
        %LogScale('x',10);
        LogScale('y',2)
        ColorbarWithAxis([0 0.5],'Phase Coupling (pMRL)')

        axis xy
        xlabel('Theta Power (mean^-1)');ylabel('freq (Hz)')
        %title((celltypes{tt}))
        if ll==1
        title((celltypes{tt}))
        end
    end 
    end
    
    NiceSave('CouplingbyThetarat',figfolder,baseName,'includeDate',true)
end
