function [popCCG, ISICCG, CellClass,ISIStats ] = PopCCGAnalysis(basePath,figfolder)
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
%%
for ss = 1:3
    %%
  % ss = 2;
    
[popCCG.(states{ss})] = PopCCG(spikes,'showfig',true,'cellclass',CellClass.label,...
    'ints',SleepState.ints.(states{ss}),'sort',ISIStats.sorts.(states{ss}).ratebyclass);

subplot(2,2,1)
caxis([0 2])
subplot(2,2,2)
caxis([0 25])
NiceSave(['PopCCG_',states{ss}],figfolder,baseName)

%%
ISIStats.allspikes.instate = cellfun(@(X) InIntervals(X,double(SleepState.ints.(states{ss}))),...
    ISIStats.allspikes.times,'UniformOutput',false);

%%
 ccgspikes = cellfun(@(X,Y) X(Y),ISIStats.allspikes.times,ISIStats.allspikes.instate,'UniformOutput',false);


nspkthresh = 25;
%clear ISICCG
for cc = 1:spikes.numcells
bz_Counter(cc,spikes.numcells,'Cell')
clear binccg
cellISIspikes = {};
for tt = 1:length(celltypes)
    popothercells = CellClass.(celltypes{tt});
    popothercells(cc) = false;
    cellISIspikes{tt} = cat(1,ccgspikes{popothercells});
end

ISICCG.logISIbins  = ISIStats.ISIhist.logbins;
for ii = 2:length(ISICCG.logISIbins)
   inbinspikes = ((log10(ISIStats.allspikes.ISIs{cc}) >= ISICCG.logISIbins(ii-1)) & ...
       (log10(ISIStats.allspikes.ISIs{cc}) <= ISICCG.logISIbins(ii))) | ...
       ((log10(ISIStats.allspikes.ISInp1{cc}) >= ISICCG.logISIbins(ii-1)) & ...
       (log10(ISIStats.allspikes.ISInp1{cc}) <= ISICCG.logISIbins(ii)));
   
    inbinspiketimes = ISIStats.allspikes.times{cc}(inbinspikes & ISIStats.allspikes.instate{cc});
    if length(inbinspiketimes)<nspkthresh
        inbinspiketimes = [];
    end
    cellISIspikes{ii+1} = inbinspiketimes;
          
end

%%
[binccg,ISICCG.t_ccg] = CCG(cellISIspikes,[],'binSize',0.002,'duration',0.5,'norm','rate'); 
for tt = 1:2
    ISICCG.(states{ss}).(celltypes{tt})(:,:,cc) = binccg(:,:,tt)./sum(~ismember(find(CellClass.(celltypes{tt})),cc));
end


end
%%
for tt = 1:2
    for tt2 = 1:2
        ISICCG.(states{ss}).popmean.(celltypes{tt}).(celltypes{tt2}) = ...
            nanmean(ISICCG.(states{ss}).(celltypes{tt})(:,:,CellClass.(celltypes{tt2})),3);
    end
end

%%

%%
figure
for tt = 1:2
    for tt2 = 1:2
        subplot(2,2,(tt-1)*2+tt2)
            imagesc(ISICCG.t_ccg,ISICCG.logISIbins,ISICCG.(states{ss}).popmean.(celltypes{tt}).(celltypes{tt2})')
            hold on
            %plot(ISIStats.ISIhist.(state).log(cc,:),ISICCG.logISIbins,'k')
            LogScale('y',10)
            
            if tt == 1
                ColorbarWithAxis([0 2.5],'E Pop Rate')
            else 
                ColorbarWithAxis([0 30],'I Pop Rate')
            end
            xlabel(['t lag (s) - ',(celltypes{tt2})]);ylabel('ISI (s)')
            
    end
    
end
NiceSave(['CCGbyISI_',(states{ss})],figfolder,baseName)

end

ISIStats = rmfield(ISIStats,'allspikes');
if isfield(ISIStats,'cellinfo')
    ISIStats = rmfield(ISIStats,'cellinfo');
end

try
    popCCG.cellinfo.regions = spikes.region;
catch
    popCCG.cellinfo.regions = nan(size(spikes.UID));
end


%% example cell
% figure
% for tt = 1:2
%     subplot(2,2,tt)
%         imagesc(t_ccg,ISIStats.ISIhist.logbins,binccg(:,:,tt)')
%         hold on
%         plot(ISIStats.ISIhist.(state).log(cc,:),ISIStats.ISIhist.logbins,'k')
%         LogScale('y',10)
%         colorbar
%         %caxis([0 5])
% end
