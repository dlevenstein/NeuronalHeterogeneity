function [ ] = MultiRecordingGammaModes(basePaths,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%   INPUTS
%       basePaths   cell array of basePaths to load from, must have
%                   GammaFit.cellinfo.mat file
%
%       (optional)
%       'saveName'
%       'saveFolder'
%       'region'
%       'clusterpar'
%
%%
% parse args
p = inputParser;
addParameter(p,'saveName',[])
addParameter(p,'saveFolder',[])
addParameter(p,'region',[])
addParameter(p,'clusterpar',false)

parse(p,varargin{:})
saveName_full = p.Results.saveName;
saveFolder = p.Results.saveFolder;
region = p.Results.region;
clusterpar = p.Results.clusterpar;

%% DEV
%Note - should be able to load GammaFit from basepath OR be given filenames
%Here: Check if basePaths are folders. 
%If yes. make GFfilenames: basePath/baseName.GammaFit.cellinfo.mat
%If no: GFfilenames = basePaths, set basepaths to be the folder in which
%each file is in...


% GFfilenames = {'20140526_277um.AnalysisResults.SharedGammaModeFitAnalysis.mat', ...
%     '20140527_421um.AnalysisResults.SharedGammaModeFitAnalysis.mat'};
% 
% GFfilenames = {'Achilles_11012013.AnalysisResults.SharedGammaModeFitAnalysis.mat',...
% 	'Achilles_10252013.AnalysisResults.SharedGammaModeFitAnalysis.mat',...
%     'Buddy_06272013.AnalysisResults.SharedGammaModeFitAnalysis.mat',...
%     'Cicero_09012014.AnalysisResults.SharedGammaModeFitAnalysis.mat',...
%     'Cicero_09102014.AnalysisResults.SharedGammaModeFitAnalysis.mat',...
%     'Cicero_09172014.AnalysisResults.SharedGammaModeFitAnalysis.mat',...
%     'Gatsby_08022013.AnalysisResults.SharedGammaModeFitAnalysis.mat',...
%     'Gatsby_08282013.AnalysisResults.SharedGammaModeFitAnalysis.mat'};

% GFfilenames = {'Achilles_11012013.AnalysisResults.SharedGammaModeFitAnalysis.mat',...
% 	'Achilles_10252013.AnalysisResults.SharedGammaModeFitAnalysis.mat'};
% saveName_full = 'CA1_test';
% basePaths = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/AnalysisScripts/AnalysisFigs/SharedGammaModeFitAnalysis';
% clusterpar = false;
% region = [];
%%
display(['Running Gamma Fit: ',saveName_full])

if ~iscell(basePaths)
    saveFolder = basePaths;
    [temp{1:length(GFfilenames)}] = deal(basePaths);
    basePaths = temp;
else
    baseNames = cellfun(@bz_BasenameFromBasepath,basePaths,'UniformOutput',false);
    GFfilenames = cellfun(@(X,Y) fullfile(X,[Y,'.GammaFit.cellinfo.mat']),basePaths,baseNames,'UniformOutput',false);
end


%ISSUE: regions! some recordings have cells from different regions...
%Keep only cells with proper region tag

%Need to keep track of.... basePath/baseName for each cell
clear LoadGF
success = true(size(GFfilenames));
for ff = 1:length(GFfilenames)
    try
        LoadGF(ff) = load(GFfilenames{ff});
    catch
        display(['Failed to load ',GFfilenames{ff}])
        success(ff) = false;
        continue
    end
    baseName{ff} = bz_BasenameFromBasepath(LoadGF(ff).GammaFit.WAKEstate.detectorinfo.detectionparms.basePath);
    saveName{ff} = [baseName{ff},'.GammaFit_full.cellinfo.mat']; %Note: add an option here for a tag (e.g. region...)
    %savefilename{ff} here: figure out the filename to re-save this GammaFit
    statenames = fieldnames(LoadGF(ff).GammaFit);
    for ss = 1:length(statenames)
        LoadGF(ff).GammaFit.(statenames{ss}).recordingIDX = ff.*ones(size(LoadGF(ff).GammaFit.(statenames{ss}).sharedfit(1).GSlogrates));
    end
end
display(['Loaded ',num2str(sum(success)),' recordings'])
LoadGF = bz_CollapseStruct(LoadGF(success),'match','justcat',true);

%%
statenames = fieldnames(LoadGF.GammaFit);
statecolors = {'b','k'};
%%
if clusterpar
    pc = parcluster('local');
    % % store temporary files in the 'scratch' drive on the cluster, labeled by job ID
    pc.JobStorageLocation = strcat(getenv('SCRATCH'), '/', getenv('SLURM_JOB_ID'));
    % % enable MATLAB to utilize the multiple cores allocated in the job script
    % % SLURM_NTASKS_PER_NODE is a variable set in the job script by the flag --tasks-per-node
    % % we use SLURM_NTASKS_PER_NODE - 1, because one of these tasks is the original MATLAB script itself
    parpool(pc, str2num(getenv('SLURM_NTASKS_PER_NODE'))-1,'IdleTimeout', Inf);
end
%%
%Consider parfor to run in parallel on cluster.
for ss = 1:find(success)
    
    %Select only the cells in the proper region
    if isempty(region)
        keepcells = true(size(LoadGF.GammaFit.(statenames{ss}).cellstats.meanrate));
    else
        keepcells = strcmp(LoadGF.GammaFit.(statenames{ss}).cellstats.region,region);
    end
    ISIdists4fit = LoadGF.GammaFit.(statenames{ss}).ISIdists(:,keepcells);
    meanFR = LoadGF.GammaFit.(statenames{ss}).cellstats.meanrate(keepcells);
    recordingIDX.(statenames{ss}) = LoadGF.GammaFit.(statenames{ss}).recordingIDX(keepcells);
    
    AScost = LoadGF.GammaFit.(statenames{ss}).detectorinfo.detectionparms.AScost_lambda(1);
    MScost = LoadGF.GammaFit.(statenames{ss}).detectorinfo.detectionparms.MScost(1);
    % MScost = 10;
    % AScost = 0.05;
    keepAS = 5;
    GammaFit_all.(statenames{ss}) = bz_FitISISharedGammaModes_new(ISIdists4fit,...
        'logtimebins',LoadGF.GammaFit.(statenames{ss}).logtimebins(1,:),...
        'maxAS',keepAS,'numAS',keepAS,'figfolder',saveFolder,...
        'AScost_lambda',AScost,'AScost_p',1,'ASguess',false,'MScost',MScost,'figname',[saveName_full,(statenames{ss})],...
        'savecellinfo',false,'forceRedetect',true,'singlefit',true,...
        'display_results','iter','meanFR',meanFR,...
        'basePath',saveFolder,'UseParallel',true);
    
    GammaFit_all.(statenames{ss}).cellstats.NW = LoadGF.GammaFit.(statenames{ss}).cellstats.NW;
end

% for ss = 1:length(statenames)
%     GammaFit_all.(statenames{ss}) = temp(ss);
% end 
%% FIgure here showing results of full fits. 
%Compare - shared fits each recording...
%%
for ff = 1:find(success)
    for ss = 1:2
        reccells = recordingIDX.(statenames{ss})==ff;
        thisrecfit = GammaFit_all.(statenames{ss});
        
        thisrecfit.singlecell = thisrecfit.singlecell(reccells);
        thisrecfit.ISIdists = thisrecfit.ISIdists(:,reccells);
        thisrecfit.numcells = sum(reccells);
        thisrecfit.costval = thisrecfit.costval(:,reccells);
        thisrecfit.cellstats.meanrate = thisrecfit.cellstats.meanrate(reccells);
        thisrecfit.cellstats.NW = thisrecfit.cellstats.NW(reccells);
        
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

%%
%load('CA1.GammaFit_all.cellinfo.mat')
%%
   whichAS.WAKEstate = 5;
   whichAS.NREMstate = 5;
   
% weightthresh = 0.01;
% clear modeweightcorr allweights numsigAS
% for pp = 1:keepAS+1
%     allweights{pp} = GammaFit_all.(statenames{ss}).sharedfit(pp).ASweights;
%     allweights{pp}(log10(allweights{pp})<-4) = 1e-4;
%     numsigAS{pp} = sum(allweights{pp}>weightthresh,2);
%     
%     modeweightcorr{pp} = corr([allweights{pp} GammaFit_all.(statenames{ss}).sharedfit(pp).GSweights'] ,...
%         'type','spearman');
% end
%%
GScolor = [0.6 0.4 0];
lowthreshcolor = [0.95 0.95 0.95];
numrepeats = 3;
%excell = excells;
histcolors = [repmat([1 1 1],numrepeats,1);makeColorMap(lowthreshcolor,[0 0 0])];

for ss = 1:2
 
   
    [~,sortGSrate] = sort(GammaFit_all.(statenames{ss}).sharedfit(whichAS.(statenames{ss})+1).GSlogrates);
    logtimebins = GammaFit_all.(statenames{ss}).logtimebins;
    numcells = GammaFit_all.(statenames{ss}).numcells;
    logISIhist = GammaFit_all.(statenames{ss}).ISIdists;
    meanFR = GammaFit_all.(statenames{ss}).cellstats.meanrate;
    costval = GammaFit_all.(statenames{ss}).costval;
    
    weightthresh = 0.01;
    for pp = 1:keepAS+1
        allweights.(statenames{ss}){pp} = GammaFit_all.(statenames{ss}).sharedfit(pp).ASweights;
        allweights.(statenames{ss}){pp}(log10(allweights.(statenames{ss}){pp})<-4) = 1e-4;
        numsigAS{pp} = sum(allweights.(statenames{ss}){pp}>weightthresh,2);

        modeweightcorr{pp} = corr([allweights.(statenames{ss}){pp} GammaFit_all.(statenames{ss}).sharedfit(pp).GSweights'] ,...
            'type','spearman');
    end
    
    
    figure
        subplot(3,3,1)
            imagesc(logtimebins,[1 numcells],logISIhist(:,sortGSrate)')
            hold on
            plot(-log10(meanFR(sortGSrate)),[1:numcells],'k.')
            plot(-GammaFit_all.(statenames{ss}).sharedfit(whichAS.(statenames{ss})+1).GSlogrates(sortGSrate),[1:numcells],'.','color',GScolor)
            %colorbar
            colormap(gca,histcolors)
            title([saveName_full,(statenames{ss})])
            xlim([-3 2])
            ylabel(['Cells (',num2str(numcells),')'])
            set(gca,'yticklabel',[])
            LogScale('x',10,'exp',true)

        subplot(3,3,2)
            plot(linspace(0,size(costval,1)-1,size(costval,1)),mean(log10(costval),2),'ko-')
            hold on
            plot(linspace(0,size(costval,1)-1,size(costval,1)),log10(costval),'k.')
            box off
            axis tight
            LogScale('y',10,'nohalf',true)
            xlabel('# Modes');ylabel('Error')

        for pp = 1:keepAS+1
        subplot(3,7,pp+7)
            bz_PlotISIDistModes(GammaFit_all.(statenames{ss}),'all','showSingleFits',true,...
                'whichShare',pp,'dotscale',10,'dotscaleAS',150)
            ylim([-1.5 1.6])
            LogScale('y',10,'nohalf',true)
            if pp>1
                set(gca,'yticklabels',[])
                ylabel('')
            end
            box off

        subplot(4,7,pp+21)
            hist(numsigAS{pp},linspace(0,pp-1,pp))
            hold on
            box off
            axis tight
            xlabel('# Modes');ylabel('# Cells') 

    %         subplot(4,7,pp)
    %             imagesc(modeweightcorr{pp})
    %             colorbar
    %             crameri('vik','pivot',0)
    %             xlabel('Mode');ylabel('Mode')
    %             set(gca,'ytick',[0:pp])
    %             set(gca,'xtick',[0:pp])
        end
        
        pp = whichAS.(statenames{ss})+1;
        subplot(3,3,3)
            imagesc(modeweightcorr{pp})
            colorbar
            crameri('vik','pivot',0)
            xlabel('Mode');ylabel('Mode')
            set(gca,'ytick',[0:pp])
            set(gca,'xtick',[0:pp])

NiceSave(['CompareNAS_',(statenames{ss})],saveFolder,saveName_full);
end

%% Relate WAKE and NREM
        modeweightcorr_NW = corr([allweights.(statenames{1}){whichAS.(statenames{1})+1} GammaFit_all.(statenames{1}).sharedfit(whichAS.(statenames{1})+1).GSweights'] ,...
            [allweights.(statenames{2}){whichAS.(statenames{2})+1} GammaFit_all.(statenames{2}).sharedfit(whichAS.(statenames{2})+1).GSweights'],...
            'type','spearman');

figure
    subplot(2,4,1)
    for ss = 1:2
        bz_PlotISIDistModes(GammaFit_all.(statenames{ss}),'all','showSingleFits',false,...
            'whichShare',whichAS.(statenames{ss})+1,'dotscale',10,'dotscaleAS',150)
    end
        ylim([-2 1.9])
        box off
        LogScale('y',10,'nohalf',true)
        
    subplot(3,3,3)
        imagesc(modeweightcorr_NW)
        colorbar
        crameri('vik','pivot',0)
        xlabel('Mode WAKE');ylabel('Mode NREM')
        set(gca,'ytick',[0:whichAS.(statenames{1})+1])
        set(gca,'xtick',[0:whichAS.(statenames{2})+1])
        
    subplot(3,3,2)
    hold on
    for ss = 1:2
        plot(linspace(0,size(GammaFit_all.(statenames{ss}).costval,1)-1,size(GammaFit_all.(statenames{ss}).costval,1)),...
            mean(log10(GammaFit_all.(statenames{ss}).costval),2),'o-','color',statecolors{ss})
    end
        box off
        axis tight
        LogScale('y',10,'nohalf',true)
        xlabel('# Modes');ylabel('Error')
        
NiceSave(['CompareStates'],saveFolder,saveName_full);
end

