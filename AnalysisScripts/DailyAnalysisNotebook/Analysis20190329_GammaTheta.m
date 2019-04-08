function [ ] = Analysis20190329(basePath,figfolder)
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
thetalfp = bz_Filter(th_lfp,'passband',[6 9]);
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
[ga_bythetarat] = ConditionalHist((thetalfp.thetadelta(thetalfp.wakeidx)),(thetalfp.gammaamp(thetalfp.wakeidx)),...
    'Xbounds',[0.1 2],'numXbins',20,'Ybounds',[0 15],'numYbins',200,'minX',400);

[ga_bythetarat_log] = ConditionalHist(log2(thetalfp.thetadelta(thetalfp.wakeidx)),log2(thetalfp.gammaamp(thetalfp.wakeidx)),...
    'Xbounds',[-2.5 2],'numXbins',20,'Ybounds',[-5 5],'numYbins',150,'minX',400);


%%
figure
subplot(3,3,1)
imagesc(ga_bytheta.Xbins,ga_bytheta.Ybins,log10(ga_bytheta.pYX)')
axis xy
LogScale('y',10)
LogScale('x',2)

subplot(3,3,4)
bar(ga_bytheta.Xbins,ga_bytheta.Xhist)
axis tight
box off

subplot(3,3,2)
a = imagesc(ga_bythetarat.Xbins,ga_bythetarat.Ybins,log10(ga_bythetarat.pYX)');
alpha(a,single((ga_bythetarat.XYhist')>5))
box off
axis xy
%LogScale('y',10)
%LogScale('x',2)
xlabel('Theta Ratio');ylabel('Ga/Rp Power (120-200Hz)')

subplot(3,3,5)
bar(ga_bythetarat.Xbins,ga_bythetarat.Xhist./thetalfp.samplingRate)
axis tight
box off



subplot(3,3,3)
a = imagesc(ga_bythetarat_log.Xbins,ga_bythetarat_log.Ybins,log10(ga_bythetarat_log.pYX)');
alpha(a,single((ga_bythetarat_log.XYhist')>5))
box off
axis xy
%LogScale('y',10)
LogScale('xy',2)
xlabel('Theta Ratio');ylabel('Ga/Rp Power (120-200Hz)')

subplot(3,3,6)
bar(ga_bythetarat_log.Xbins,ga_bythetarat_log.Xhist./thetalfp.samplingRate)
axis tight
box off

    NiceSave('GammaPowerByThetaRat',figfolder,baseName,'includeDate',true)

%%

%% ISI Coupling conditioned on theta, theta/delta

state = states{1};
%ints = SleepState.ints.(state);

%Take only subset of time (random intervals) so wavelets doesn't break
%computer (total 625s)
usetime = 6000;%2500
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

%%

    LFPCoupling_thetarat(cc) = bz_ConditionalLFPCoupling( ISIStats.allspikes,ISIStats.allspikes.thetarat,wavespec,...
        'Xbounds',[0.1 2],'intervals',windows,'showFig',true,'numXbins',20,...
    'minX',25,'CellClass',CellClass,'spikeLim',20000,...
    'showFig',true,'binNorm',true);
%%
    LFPCoupling_thetarat_log(cc) = bz_ConditionalLFPCoupling( ISIStats.allspikes,ISIStats.allspikes.thetarat_log,wavespec,...
        'Xbounds',[-2.5 2],'intervals',windows,'showFig',true,'numXbins',20,...
    'minX',25,'CellClass',CellClass,'spikeLim',20000,...
    'showFig',true,'binNorm',true);

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
clear laminarLFPCoupling_thetarat_log
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
    
	laminarLFPCoupling_thetarat_log.(lamina.names{ll}).Xbins = LFPCoupling_thetarat_log(1).Xbins;
    laminarLFPCoupling_thetarat_log.(lamina.names{ll}).freqs	 = LFPCoupling_thetarat_log(1).freqs;


    for cc = 1:spikes.numcells
    laminarLFPCoupling_thetarat.(lamina.names{ll}).mrl(:,:,cc) = ...
        LFPCoupling_thetarat(ss,lamina.(lamina.names{ll}).cellchan(cc)+offset).mrl(:,:,cc);
    laminarLFPCoupling_thetarat.(lamina.names{ll}).meanpower(:,:,cc) = ...
        LFPCoupling_thetarat(ss,lamina.(lamina.names{ll}).cellchan(cc)+offset).meanpower(:,:,cc);
    laminarLFPCoupling_thetarat.(lamina.names{ll}).mutInfoXPower(cc,:) = ...
        LFPCoupling_thetarat(ss,lamina.(lamina.names{ll}).cellchan(cc)+offset).mutInfoXPower(cc,:);
    
    laminarLFPCoupling_thetarat_log.(lamina.names{ll}).mrl(:,:,cc) = ...
        LFPCoupling_thetarat_log(ss,lamina.(lamina.names{ll}).cellchan(cc)+offset).mrl(:,:,cc);
    laminarLFPCoupling_thetarat_log.(lamina.names{ll}).meanpower(:,:,cc) = ...
        LFPCoupling_thetarat_log(ss,lamina.(lamina.names{ll}).cellchan(cc)+offset).meanpower(:,:,cc);
    laminarLFPCoupling_thetarat_log.(lamina.names{ll}).mutInfoXPower(cc,:) = ...
        LFPCoupling_thetarat_log(ss,lamina.(lamina.names{ll}).cellchan(cc)+offset).mutInfoXPower(cc,:);

    end


    %Get Cell Type Average
    for tt = 1:length(celltypes)
        laminarLFPCoupling_thetarat.(lamina.names{ll}).allmeanpower.(celltypes{tt}) = ...
            nanmean(laminarLFPCoupling_thetarat.(lamina.names{ll}).meanpower(:,:,CellClass.(celltypes{tt})),3);
        laminarLFPCoupling_thetarat.(lamina.names{ll}).almeanpMRL.(celltypes{tt}) = ...
            nanmean(laminarLFPCoupling_thetarat.(lamina.names{ll}).mrl(:,:,CellClass.(celltypes{tt})),3);
        laminarLFPCoupling_thetarat.(lamina.names{ll}).groupmutinf.(celltypes{tt}) = ...
            nanmean(laminarLFPCoupling_thetarat.(lamina.names{ll}).mutInfoXPower(CellClass.(celltypes{tt}),:),1);

        laminarLFPCoupling_thetarat_log.(lamina.names{ll}).allmeanpower.(celltypes{tt}) = ...
            nanmean(laminarLFPCoupling_thetarat_log.(lamina.names{ll}).meanpower(:,:,CellClass.(celltypes{tt})),3);
        laminarLFPCoupling_thetarat_log.(lamina.names{ll}).almeanpMRL.(celltypes{tt}) = ...
            nanmean(laminarLFPCoupling_thetarat_log.(lamina.names{ll}).mrl(:,:,CellClass.(celltypes{tt})),3);
        laminarLFPCoupling_thetarat_log.(lamina.names{ll}).groupmutinf.(celltypes{tt}) = ...
            nanmean(laminarLFPCoupling_thetarat_log.(lamina.names{ll}).mutInfoXPower(CellClass.(celltypes{tt}),:),1);
    end

end



    
%% Figure Theta Ratio
powermap = makeColorMap([0 0 0.8],[1 1 1],[0.8 0 0]);

    figure
       
    
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
        xlabel('Theta Ratio (mean^-^1)');ylabel('freq (Hz)')
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
        ColorbarWithAxis([0 0.3],'Phase Coupling (pMRL)')

        axis xy
        xlabel('Theta Ratio (mean^-^1)');ylabel('freq (Hz)')
        %title((celltypes{tt}))
        if ll==1
        title((celltypes{tt}))
        end
    end 
    end
    
    
    subplot(4,3,1)
        a = imagesc(ga_bythetarat.Xbins,ga_bythetarat.Ybins,log10(ga_bythetarat.pYX)');
        alpha(a,single((ga_bythetarat.XYhist')>5))
        box off
        axis xy
        %LogScale('y',10)
        %LogScale('x',2)
        xlabel('Theta Ratio');ylabel('Ga/Rp Power (120-200Hz)')

	subplot(4,3,4)
        bar(ga_bythetarat.Xbins,ga_bythetarat.Xhist./thetalfp.samplingRate)
        axis tight
        box off
        xlabel('Theta Ratio');ylabel('Occupancy (s)')





    
    NiceSave('CouplingbyThetarat',figfolder,baseName,'includeDate',true)
    
    
    
%% Figure Theta Ratio Log
    figure
       
    
    for ll=1:2
    for tt = 1:length(celltypes)
    subplot(4,3,tt+(ll*3)+4)
    colormap(gca,powermap)
        imagesc(laminarLFPCoupling_thetarat_log.(lamina.names{ll}).Xbins,log2(laminarLFPCoupling_thetarat_log.(lamina.names{ll}).freqs),...
            log2(laminarLFPCoupling_thetarat_log.(lamina.names{ll}).allmeanpower.(celltypes{tt}))')
        colorbar
        ColorbarWithAxis([-1.25 1.25],'Power (mean^-^1)')
        %LogScale('x',10);
        LogScale('xy',2)
        LogScale('c',2)
        axis xy
        xlabel('Theta Ratio (mean^-^1)');ylabel('freq (Hz)')
    end  
    
    
    for tt = 1:length(celltypes)
    subplot(4,3,tt+(ll*3)-2)
        imagesc(laminarLFPCoupling_thetarat_log.(lamina.names{ll}).Xbins,log2(laminarLFPCoupling_thetarat_log.(lamina.names{ll}).freqs),...
            laminarLFPCoupling_thetarat_log.(lamina.names{ll}).almeanpMRL.(celltypes{tt})')
        colorbar
        hold on
        %caxis([0.5 1.5])
        %LogScale('x',10);
        LogScale('xy',2)
        ColorbarWithAxis([0 0.3],'Phase Coupling (pMRL)')

        axis xy
        xlabel('Theta Ratio (mean^-^1)');ylabel('freq (Hz)')
        %title((celltypes{tt}))
        if ll==1
        title((celltypes{tt}))
        end
    end 
    end
    
   


forcolor = ga_bythetarat_log.pYX;forcolor((ga_bythetarat_log.XYhist)<=20)=0;
    subplot(4,3,1)
        a = imagesc(ga_bythetarat_log.Xbins,ga_bythetarat_log.Ybins,log10(forcolor)');
        alpha(a,single((ga_bythetarat_log.XYhist')>20))
        box off
        axis xy
        %LogScale('y',10)
        LogScale('xy',2)
        xlabel('Theta Ratio');ylabel('Ga/Rp Power (120-200Hz)')

    subplot(4,3,4)
        bar(ga_bythetarat_log.Xbins,ga_bythetarat_log.Xhist./thetalfp.samplingRate)
        axis tight
        box off

    
    NiceSave('CouplingbyThetarat_log',figfolder,baseName,'includeDate',true)
end
