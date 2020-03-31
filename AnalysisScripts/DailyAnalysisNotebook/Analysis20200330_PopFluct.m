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
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
basePath = '/Users/dl2820/Dropbox/Research/Datasets/Cicero_09102014';
%basePath = pwd;
%basePath = fullfile(reporoot,'Datasets/onProbox/AG_HPC/Achilles_11012013');
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
statecolors = {'k','b','r'};

try
    celltypes = CellClass.celltypes;
catch
    celltypes = unique(CellClass.label);
end
cellcolor = {'k','r'};

%%
ISIStats.allspikes.ISInp1 = cellfun(@(X) [X(2:end);nan],ISIStats.allspikes.ISIs,...
    'UniformOutput',false);
%% Restrict to state

state = states{1};
ISIStats.allspikes.instate = cellfun(@(X) InIntervals(X,double(SleepState.ints.(state))),...
    ISIStats.allspikes.times,'UniformOutput',false);


%% CCGs
ccgspikes = cellfun(@(X,Y) X(Y),ISIStats.allspikes.times,ISIStats.allspikes.instate,'UniformOutput',false);
[allccg,t_ccg] = CCG(ccgspikes,[],'binSize',0.001,'duration',0.4,'norm','rate'); 

popspikes = {cat(1,ccgspikes{CellClass.pE}),cat(1,ccgspikes{CellClass.pI})};
[popccg,t_ccg] = CCG(popspikes,[],'binSize',0.001,'duration',0.4,'norm','rate'); 

%%
clear meanCCG
for cc = 1:spikes.numcells
    for tt = 1:length(celltypes)
        popothercells = CellClass.(celltypes{tt});
        popothercells(cc) = false;
        meanCCG.(celltypes{tt})(:,cc) = mean(allccg(:,cc,popothercells),3);
    end
end
%%
figure
for tt = 1:2
    subplot(2,2,tt)
        imagesc(t_ccg,[0 spikes.numcells],log10(meanCCG.(celltypes{tt})(:,ISIStats.sorts.(state).ratebyclass))')
        hold on
        plot(t_ccg,bz_NormToRange(-popccg(:,tt,1),[0 sum(CellClass.pE)]),'color',cellcolor{tt})
        plot(t_ccg,bz_NormToRange(-popccg(:,tt,2),sum(CellClass.pE)+[1 sum(CellClass.pI)]),'color',cellcolor{tt})
        plot(xlim(gca),sum(CellClass.pE).*[1 1],'w')
        colorbar
end

%%
nspkthresh = 10;
clear ISICCG
for cc = 1:spikes.numcells
bz_Counter(cc,spikes.numcells,'Cell')
clear binccg
cellISIspikes = {};
for tt = 1:length(celltypes)
    popothercells = CellClass.(celltypes{tt});
    popothercells(cc) = false;
    cellISIspikes{tt} = cat(1,ccgspikes{popothercells});
end

for ii = 2:length(ISIStats.ISIhist.logbins)
   inbinspikes = ((log10(ISIStats.allspikes.ISIs{cc}) >= ISIStats.ISIhist.logbins(ii-1)) & ...
       (log10(ISIStats.allspikes.ISIs{cc}) <= ISIStats.ISIhist.logbins(ii))) | ...
       ((log10(ISIStats.allspikes.ISInp1{cc}) >= ISIStats.ISIhist.logbins(ii-1)) & ...
       (log10(ISIStats.allspikes.ISInp1{cc}) <= ISIStats.ISIhist.logbins(ii)));
   
    inbinspiketimes = ISIStats.allspikes.times{cc}(inbinspikes & ISIStats.allspikes.instate{cc});
    if length(inbinspiketimes)<nspkthresh
        inbinspiketimes = [];
    end
    cellISIspikes{ii+1} = inbinspiketimes;
          
end

%%
[binccg,ISICCG.t_ccg] = CCG(cellISIspikes,[],'binSize',0.002,'duration',0.5,'norm','rate'); 
for tt = 1:2

    ISICCG.(celltypes{tt})(:,:,cc) = binccg(:,:,tt)./sum(~ismember(find(CellClass.(celltypes{tt})),cc));
end


end
%%
for tt = 1:2
    for tt2 = 1:2
        ISICCG.popmean.(celltypes{tt}).(celltypes{tt2}) = ...
            nanmean(ISICCG.(celltypes{tt})(:,:,CellClass.(celltypes{tt2})),3);
    end
end
%%
figure
for tt = 1:2
    for tt2 = 1:2
        subplot(2,2,(tt-1)*2+tt2)
            imagesc(t_ccg,ISIStats.ISIhist.logbins,ISICCG.popmean.(celltypes{tt}).(celltypes{tt2})')
            hold on
            %plot(ISIStats.ISIhist.(state).log(cc,:),ISIStats.ISIhist.logbins,'k')
            LogScale('y',10)
            colorbar
            if tt == 1
                clim([0 2.5])
            else 
                clim([0 30])
            end
    end
    
end
%%
figure
for tt = 1:2
    subplot(2,2,tt)
        imagesc(t_ccg,ISIStats.ISIhist.logbins,binccg(:,:,tt)')
        hold on
        plot(ISIStats.ISIhist.(state).log(cc,:),ISIStats.ISIhist.logbins,'k')
        LogScale('y',10)
        colorbar
        %caxis([0 5])
end
