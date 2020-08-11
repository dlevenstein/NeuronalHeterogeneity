function [GammaFit] = SharedGammaModeFitAnalysis(basePath,figfolder)
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
%reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/';
%basePath = pwd;
%basePath = '/Users/dl2820/Dropbox/research/Datasets/Cicero_09102014';
%basePath = '/Users/dlevenstein/Dropbox/research/Datasets/20140526_277um';
%figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);
SAVECELLINFO = true;

%Load Stuff
sessionInfo = bz_getSessionInfo(basePath);
spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
%ISIStats = bz_LoadCellinfo(basePath,'ISIStats');

try
    celltypes = CellClass.celltypes;
catch
    celltypes = unique(CellClass.label);
end
cellcolor = {'k','r'};
statecolor = {'b','k','r'};

%%
statenames = {'NREMstate','WAKEstate','REMstate'};
%%

[ ISIstats ] = bz_ISIStats( spikes,'ints',SleepState.ints,...
    'cellclass',CellClass.label,'shuffleCV2',false,...
    'savecellinfo',false,'basePath',basePath,'forceRedetect',true,...
    'numISIbins',150,'logISIbounds',[0.0001 500]);

cellinfofilename = fullfile(basePath,[baseName,'.GammaFit.cellinfo.mat']);
%%
numAS.NREMstate = 5;
numAS.WAKEstate = 5;
numAS.REMstate = 5;

%%
spkthresh = 250;
for ss = 1:3
    numspks = cellfun(@sum,ISIstats.allspikes.instate.(statenames{ss}));
    logtimebins = ISIstats.ISIhist.logbins;
    logISIhist = ISIstats.ISIhist.(statenames{ss}).log;
    usecells{ss} = find(CellClass.pE & numspks>spkthresh);
    logISIhist = logISIhist(usecells{ss},:)';
    logISIhist = logISIhist./mode(diff(logtimebins));
    GammaFit.(statenames{ss}) = bz_FitISISharedGammaModes(logISIhist,logtimebins,...
        'numAS',numAS.(statenames{ss}),...
        'figfolder','detectionfigures','basePath',basePath,...
        'AScost_lambda',0.13,'AScost_p',1/2,'ASguess',true,'MScost',3,'figname',(statenames{ss}));
    
    
    GammaFit.(statenames{ss}).cellstats.meanrate = ...
        ISIstats.summstats.(statenames{ss}).meanrate(usecells{ss});
    GammaFit.(statenames{ss}).cellstats.UID = spikes.UID(usecells{ss});
    if isfield(spikes,'region')
        GammaFit.(statenames{ss}).cellstats.region = spikes.region(usecells{ss});
    end
    if length(GammaFit.(statenames{ss}).cellstats.meanrate) ~= ...
            length(GammaFit.(statenames{ss}).singlecell)
        error('bad number of cells')
    end
    


end


if SAVECELLINFO
    save(cellinfofilename,'GammaFit')
end



%%



end
