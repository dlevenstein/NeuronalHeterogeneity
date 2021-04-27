function [ ] = MultiRecordingGammaModes(basePaths,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%
% parse args
p = inputParser;
addParameter(p,'saveName',[])
addParameter(p,'saveFolder',[])

parse(p,varargin{:})
saveName_full = p.Results.saveName;
saveFolder = p.Results.saveFolder;

%% DEV
%Note - should be able to load GammaFit from basepath OR be given filenames
%Here: Check if basePaths are folders. 
%If yes. make GFfilenames: basePath/baseName.GammaFit.cellinfo.mat
%If no: GFfilenames = basePaths, set basepaths to be the folder in which
%each file is in...


GFfilenames = {'20140526_277um.AnalysisResults.SharedGammaModeFitAnalysis.mat', ...
    '20140527_421um.AnalysisResults.SharedGammaModeFitAnalysis.mat'};

GFfilenames = {'Achilles_11012013.AnalysisResults.SharedGammaModeFitAnalysis.mat',...
	'Achilles_10252013.AnalysisResults.SharedGammaModeFitAnalysis.mat',...
    'Buddy_06272013.AnalysisResults.SharedGammaModeFitAnalysis.mat',...
    'Cicero_09012014.AnalysisResults.SharedGammaModeFitAnalysis.mat',...
    'Cicero_09102014.AnalysisResults.SharedGammaModeFitAnalysis.mat',...
    'Cicero_09172014.AnalysisResults.SharedGammaModeFitAnalysis.mat',...
    'Gatsby_08022013.AnalysisResults.SharedGammaModeFitAnalysis.mat',...
    'Gatsby_08282013.AnalysisResults.SharedGammaModeFitAnalysis.mat'};

if ~iscell(basePaths)
    saveFolder = basePaths;
    [temp{1:length(GFfilenames)}] = deal(basePaths);
    basePaths = temp;
end


%ISSUE: regions! some recordings have cells from different regions...

%Need to keep track of.... basePath/baseName for each cell
clear LoadGF
for ff = 1:length(GFfilenames)
    LoadGF(ff) = load(GFfilenames{ff});
    baseName{ff} = bz_BasenameFromBasepath(LoadGF(ff).GammaFit.WAKEstate.detectorinfo.detectionparms.basePath);
    saveName{ff} = [baseName{ff},'.GammaFit_full.cellinfo.mat'];
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
pc = parcluster('local');
% % store temporary files in the 'scratch' drive on the cluster, labeled by job ID
pc.JobStorageLocation = strcat(getenv('SCRATCH'), '/', getenv('SLURM_JOB_ID'));
% % enable MATLAB to utilize the multiple cores allocated in the job script
% % SLURM_NTASKS_PER_NODE is a variable set in the job script by the flag --tasks-per-node
% % we use SLURM_NTASKS_PER_NODE - 1, because one of these tasks is the original MATLAB script itself
parpool(pc, str2num(getenv('SLURM_NTASKS_PER_NODE'))-1);

%Consider parfor to run in parallel on cluster.
parfor ss = 1:length(statenames)
    AScost = LoadGF.GammaFit.(statenames{ss}).detectorinfo.detectionparms.AScost_lambda(1);
    MScost = LoadGF.GammaFit.(statenames{ss}).detectorinfo.detectionparms.MScost(1);
    % MScost = 10;
    % AScost = 0.05;
    keepAS = 6;
    temp(ss) = bz_FitISISharedGammaModes_new(LoadGF.GammaFit.(statenames{ss}).ISIdists,...
        'logtimebins',LoadGF.GammaFit.(statenames{ss}).logtimebins(1,:),...
        'maxAS',keepAS,'numAS',keepAS,'figfolder',saveFolder,...
        'AScost_lambda',AScost,'AScost_p',1,'ASguess',false,'MScost',MScost,'figname',[saveName_full,(statenames{ss})],...
        'savecellinfo',false,'forceRedetect',true,'singlefit',true,...
        'display_results','iter','meanFR',LoadGF.GammaFit.(statenames{ss}).cellstats.meanrate,...
        'basePath',saveFolder);
end

for ss = 1:length(statenames)
    GammaFit_all.(statenames{ss}) = temp(ss);
end 
%% FIgure here showing results of full fits. 
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
        GammaFit_eachrec(ff).(statenames{ss}) = thisrecfit;
    end
end

%%
%Save the GammaFit_all file somewhere...
cellinfofilename = fullfile(saveFolder,[saveName_full,'.GammaFit_all.cellinfo.mat']); 
save(cellinfofilename,'GammaFit_all')

%Each recordings cellinfo file
for ff = 1:length(GFfilenames)
    %GammaFit_full = GammaFit_full(ff)
    cellinfofilename = fullfile(basePaths{ff},saveName{ff}); 
    GammaFit_full = GammaFit_eachrec(ff);
    save(cellinfofilename,'GammaFit_full')
end

end

