function [ ] = MultiRecordingGammaModes(basePaths)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% DEV
%Note - should be able to load GammaFit from basepath OR be given filenames
GFfilenames = {'20140526_277um.AnalysisResults.SharedGammaModeFitAnalysis.mat', ...
    '20140527_421um.AnalysisResults.SharedGammaModeFitAnalysis.mat'};

%ISSUE: regions!

%Need to keep track of.... basePath/baseName for each cell
clear LoadGF
for ff = 1:length(GFfilenames)
    LoadGF(ff) = load(GFfilenames{ff});
    %savefilename{ff} here: figure out the filename to re-save this GammaFit
    statenames = fieldnames(LoadGF(ff).GammaFit);
    for ss = 1:length(statenames)
        LoadGF(ff).GammaFit.(statenames{ss}).recordingIDX = ff.*ones(size(LoadGF(ff).GammaFit.(statenames{ss}).sharedfit.GSlogrates));
    end
end
LoadGF = bz_CollapseStruct(LoadGF,'match','justcat',true);

%%
statenames = fieldnames(LoadGF.GammaFit);
%%
%Consider parfor to run in parallel on cluster.
for ss = 1:length(statenames)

AScost = LoadGF.GammaFit.(statenames{ss}).parms.AScost_lambda(1);
MScost = LoadGF.GammaFit.(statenames{ss}).parms.MScost(1);
% MScost = 10;
% AScost = 0.05;
keepAS = 2;
GammaFit_all.(statenames{ss}) = bz_FitISISharedGammaModes_new(LoadGF.GammaFit.(statenames{ss}).ISIdists,'logtimebins',LoadGF.GammaFit.(statenames{ss}).logtimebins(1,:),...
    'maxAS',keepAS,'numAS',keepAS,...
    'AScost_lambda',AScost,'AScost_p',1,'ASguess',false,'MScost',MScost,'figname',(statenames{ss}),...
    'savecellinfo',false,'forceRedetect',true,'singlefit',true,...
    'display_results','iter','meanFR',LoadGF.GammaFit.(statenames{ss}).cellstats.meanrate);
end
%% FIgure here showing results of full fits. 
%Consider, moving fit figure into bz_FitISISharedGammaModes_new function
%Compare - shared fits each recording...
%%
for ff = 1:length(GFfilenames)
    for ss = 1:2
        reccells = LoadGF.GammaFit.(statenames{ss}).recordingIDX==ff;
        thisrecfit = GammaFit_all.(statenames{ss});
        
        thisrecfit.singlecell = thisrecfit.singlecell(reccells);
        thisrecfit.ISIdists = thisrecfit.ISIdists(:,reccells);
        thisrecfit.numcells = sum(reccells);
        thisrecfit.costval = thisrecfit.costval(:,reccells);
        thisrecfit.cellstats.meanrate = thisrecfit.cellstats.meanrate(reccells);
        
        numsharedfits = length(thisrecfit.sharedfit);
        for sf = 1:numsharedfits
            thisrecfit.sharedfit(sf).GSlogrates = thisrecfit.sharedfit(sf).GSlogrates(reccells);
            thisrecfit.sharedfit(sf).GSCVs = thisrecfit.sharedfit(sf).GSCVs(reccells);
            thisrecfit.sharedfit(sf).GSweights = thisrecfit.sharedfit(sf).GSweights(reccells);
            thisrecfit.sharedfit(sf).ASweights = thisrecfit.sharedfit(sf).ASweights(reccells,:);
        end
        GammaFit_full(ff).(statenames{ss}) = thisrecfit;
    end
end

%%
%Save each GammaFit_full in a baseName.GammaFit_full.cellinfo.mat file
%And the GammaFit_all file somewhere...
end

