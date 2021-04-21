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
%basePath = '/Users/dl2820/Dropbox/research/Datasets/20140526_277um';
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

for ss = 1:2
%ss = 1;
AScost = 0.05; %Formerly 0.05
MScost = 10;  %Formerly 10
%Here: fit with all the stuff (final parms)
keepAS = 6;
GammaFit.(statenames{ss}) = bz_FitISISharedGammaModes_new(spikes,...
    'figfolder',figfolder,'basePath',basePath,'ints',SleepState.ints.(statenames{ss}),...
    'usecells',CellClass.pE,'maxAS',keepAS,'numAS',keepAS,...
    'AScost_lambda',AScost,'AScost_p',1,'ASguess',false,'MScost',MScost,'figname',(statenames{ss}),...
    'savecellinfo',false,'forceRedetect',true,'singlefit',true,...
    'display_results','final');

%     %Here: 'holdweights', false
%     GammaFit_noweights = bz_FitISISharedGammaModes_new(spikes,...
%         'figfolder',figfolder,'basePath',basePath,'ints',SleepState.ints.(statenames{ss}),...
%         'usecells',CellClass.pE,'maxAS',6,'holdweights',false,...
%         'AScost_lambda',0.05,'AScost_p',1,'ASguess',false,'MScost',10,'figname',[(statenames{ss}),'noweights']);
%%
%For Testing weight holdover
TESTHOLD = false;
if TESTHOLD
    %Here: new hold old weights
    GammaFit_holdweights = bz_FitISISharedGammaModes_new(spikes,...
        'figfolder',figfolder,'basePath',basePath,'ints',SleepState.ints.(statenames{ss}),...
        'usecells',CellClass.pE,'maxAS',6,'holdweights',true,...
        'AScost_lambda',0.05,'AScost_p',1,'ASguess',false,'MScost',10,'figname',[(statenames{ss}),'holdweights']);

    %Compare: compute time, which/number of modes, error
    %%
    figure
    subplot(3,4,1)
        plot(GammaFit.(statenames{ss}).computetime)
        hold on
        plot(GammaFit_noweights.computetime)
        plot(GammaFit_holdweights.computetime)
        xlabel('Iteration');ylabel('Compute Time')
        legend('Hold Weights','Reset Weights','RenormWeights','location','northwest');

    subplot(3,4,5)
        plot(mean(log10(GammaFit.(statenames{ss}).costval),2))
        hold on
        plot(mean(log10(GammaFit_noweights.costval),2))
        plot(mean(log10(GammaFit_holdweights.costval),2))
        xlabel('Iteration');ylabel('Mean Loss')
        legend('Hold Weights','Reset Weights','RenormWeights','location','northwest');

        panels = [4 5 7];
    for pp = 1:3        
         subplot(3,4,1+pp)
            bz_PlotISIDistModes(GammaFit.(statenames{ss}),GammaFit.(statenames{ss}).cellstats.UID,...
                'whichShare',panels(pp),'dotscale',50)

         subplot(3,4,5+pp)
            bz_PlotISIDistModes(GammaFit_noweights,GammaFit_noweights.cellstats.UID,...
                'whichShare',panels(pp),'dotscale',50)

         subplot(3,4,9+pp)
            bz_PlotISIDistModes(GammaFit_holdweights,GammaFit_holdweights.cellstats.UID,...
                'whichShare',panels(pp),'dotscale',50)
    end
    NiceSave(['CompareInit_',(statenames{ss})],figfolder,baseName);
end
%% Example figures for process
exUID = 10; %5, 27, 95 problem children
exUID = 79; %good cicero cell 

exUID = randsample(GammaFit.(statenames{ss}).cellstats.UID,1);

% bz_PlotISIDistModes(GammaFit.(statenames{ss}),GammaFit.(statenames{ss}).cellstats.UID,'whichShare',pp)
 lowthreshcolor = [0.95 0.95 0.95];
numrepeats = 3;
%excell = excells;
histcolors = [repmat([1 1 1],numrepeats,1);makeColorMap(lowthreshcolor,[0 0 0])];   

figure
for pp = 1:keepAS+1
    fitISI = GSASmodel(GammaFit.(statenames{ss}).sharedfit(pp),...
        GammaFit.(statenames{ss}).taubins,GammaFit.(statenames{ss}).numcells,pp-1);
    [~,sortGSrate] = sort(GammaFit.(statenames{ss}).sharedfit(pp).GSlogrates);

    subplot(3,7,pp)
        imagesc(GammaFit.(statenames{ss}).logtimebins,[1 GammaFit.(statenames{ss}).numcells],fitISI(:,sortGSrate)')
        hold on
        %plot(log10(MSthresh).*[1 1],ylim(gca),'r')
        %plot(logtimebins,-bz_NormToRange(meanISIdist_fit,0.3)+numcells,'k','linewidth',2)
        %colorbar
        colormap(gca,histcolors)
        xlim([-3 2])
        set(gca,'xtick',[]);set(gca,'ytick',[])
        title([num2str(pp-1),' AS Modes'])

    subplot(3,7,pp+7)
        bz_PlotISIDistModes(GammaFit.(statenames{ss}),GammaFit.(statenames{ss}).cellstats.UID,...
            'whichShare',pp,'dotscale',10,'dotscaleAS',150)
        ylim([-1.5 1.6])
        LogScale('y',10,'nohalf',true)
        if pp>1
            set(gca,'yticklabels',[])
            ylabel('')
        end
        box off

    subplot(3,7,pp+14)
        bz_PlotISIDistModes(GammaFit.(statenames{ss}),exUID,'whichShare',pp)
        ylim([-1.5 1.6])
        LogScale('y',10,'nohalf',true)
        if pp>1
            set(gca,'yticklabels',[])
            ylabel('')
        end
        box off
end
NiceSave(['ModeExamples_',(statenames{ss})],figfolder,baseName);

%%


%%
    weightthresh = 0.01;
    clear modeweightcorr allweights numsigAS
for pp = 1:keepAS+1
    allweights{pp} = GammaFit.(statenames{ss}).sharedfit(pp).ASweights;
    allweights{pp}(log10(allweights{pp})<-4) = 1e-4;
    numsigAS{pp} = sum(allweights{pp}>weightthresh,2);
    
    modeweightcorr{pp} = corr([allweights{pp} GammaFit.(statenames{ss}).sharedfit(pp).GSweights'] ,...
        'type','spearman');
end

%%

figure
for pp = 1:keepAS+1       

    subplot(4,7,pp)
        imagesc(modeweightcorr{pp})
        colorbar
        crameri('vik','pivot',0)
        xlabel('Mode');ylabel('Mode')
        set(gca,'ytick',[0:pp])
        set(gca,'xtick',[0:pp])

    subplot(4,7,pp+7)
        hist(log10(allweights{pp}(:)))
        hold on;box off; axis tight
        plot(log10(weightthresh).*[1 1],ylim(gca),'k--')
        LogScale('x',10)
        xlabel('Weight');ylabel('# Modes (All Cells)')

    subplot(4,7,pp+14)
        hist(numsigAS{pp},linspace(0,pp-1,pp))
        hold on
        box off
        axis tight
        xlabel('# Modes');ylabel('# Cells') 
end
NiceSave(['ModeWeightProperties_',(statenames{ss})],figfolder,baseName);

%%
TESTCONSTRAINTS = false;
if TESTCONSTRAINTS
    %Examples with no constraints and example constraints? Input initial
    %For determining best constraint range
    %conditions from Fit...
    GammaFit_nocon = bz_FitISISharedGammaModes_new(spikes,...
        'figfolder',figfolder,'basePath',basePath,'ints',SleepState.ints.(statenames{ss}),...
        'usecells',CellClass.pE,...
        'AScost_lambda',0,'AScost_p',1,'ASguess',false,'MScost',0,'figname',[(statenames{ss}),'_noconst'],...
        'init_struct',GammaFit.(statenames{ss}).sharedfit(4));
 
    %%
    c_ref = [0.3 1 3 10 30 100 300];
    numAS = 3;
    %GammaFit_ref(1) = GammaFit_nocon;
    for cc = 1:7
    GammaFit_ref(cc) = bz_FitISISharedGammaModes_new(spikes,...
        'figfolder',figfolder,'basePath',basePath,'ints',SleepState.ints.(statenames{ss}),...
        'usecells',CellClass.pE,...
        'AScost_lambda',0.1,'AScost_p',1,'ASguess',false,'MScost',c_ref(cc),...
        'init_struct',GammaFit.(statenames{ss}).sharedfit(numAS+1),...
        'forceRedetect',true);
    end
    
    %% Figure - Refractory Constraint
    
    GammaFit_ref_all = bz_CollapseStruct(GammaFit_ref,1,'justcat',true);
    figure
    subplot(3,3,1)
        plot(log10(c_ref),log10(GammaFit_ref_all.costval),'k.')
        hold on
        plot(log10(c_ref),mean(log10(GammaFit_ref_all.costval),2),'ko-')
        box off
        xlabel('Refractory Cost');ylabel('Loss')
        LogScale('xy',10)
    subplot(3,3,2) 
        plot(log10(c_ref),(GammaFit_ref_all.sharedfit.GSCVs),'k.')
        hold on
        plot(log10(c_ref),mean((GammaFit_ref_all.sharedfit.GSCVs),2),'ko-')
        box off
        xlabel('Refractory Cost');ylabel('GS CV')
        LogScale('x',10)
        
    subplot(3,3,3) 
        plot(log10(c_ref),(GammaFit_ref_all.sharedfit.GSweights),'k.')
        hold on
        plot(log10(c_ref),mean((GammaFit_ref_all.sharedfit.GSweights),2),'ko-')
        box off
        xlabel('Refractory Cost');ylabel('GS Weight')
        LogScale('x',10)
        
    subplot(3,3,4) 
        plot(log10(c_ref),(GammaFit_ref_all.sharedfit.GSlogrates),'k.')
        hold on
        plot(log10(c_ref),mean((GammaFit_ref_all.sharedfit.GSlogrates),2),'ko-')
        box off
        xlabel('Refractory Cost');ylabel('GS Rate')
        LogScale('x',10)
        
    subplot(3,3,5)    
        plot(-GammaFit_ref_all.sharedfit.ASlogrates',GammaFit_ref_all.sharedfit.ASCVs','.')
        
    subplot(3,3,7)
        bz_PlotISIDistModes(GammaFit_ref(1),GammaFit_ref(1).cellstats.UID,...
            'dotscale',10,'dotscaleAS',150)
        xlim([-3 2])
        
    subplot(3,3,8)
        bz_PlotISIDistModes(GammaFit_ref(5),GammaFit_ref(5).cellstats.UID,...
            'dotscale',10,'dotscaleAS',150)
        xlim([-3 2])
        
    subplot(3,3,9)
        bz_PlotISIDistModes(GammaFit_ref(7),GammaFit_ref(7).cellstats.UID,...
            'dotscale',10,'dotscaleAS',150)
        xlim([-3 2])

    NiceSave('RefractoryCost',figfolder,baseName);
    %%
    c_AS = [0 0.05 0.1 0.15 0.2 0.25];
    for cc = 1:6
    GammaFit_AS(cc) = bz_FitISISharedGammaModes_new(spikes,...
        'figfolder',figfolder,'basePath',basePath,'ints',SleepState.ints.(statenames{ss}),...
        'usecells',CellClass.pE,...
        'AScost_lambda',c_AS(cc),'AScost_p',1,'ASguess',false,'MScost',30,...
        'init_struct',GammaFit.(statenames{ss}).sharedfit(numAS+1),...
        'forceRedetect',true);
    end
    
    %% Figure - AS Constraint
    
    GammaFit_AS_all = bz_CollapseStruct(GammaFit_AS,1,'justcat',true);
    figure
    subplot(3,3,1)
        plot(c_AS,log10(GammaFit_AS_all.costval),'k.')
        hold on
        plot(c_AS,mean(log10(GammaFit_AS_all.costval),2),'ko-')
        box off
        xlabel('AS Weight Cost');ylabel('Loss')
       % LogScale('xy',10)
    subplot(3,3,2) 
        plot(c_AS,(GammaFit_AS_all.sharedfit.GSCVs),'k.')
        hold on
        plot(c_AS,mean((GammaFit_AS_all.sharedfit.GSCVs),2),'ko-')
        box off
        xlabel('AS Weight Cost');ylabel('GS CV')
      %  LogScale('x',10)
        
    subplot(3,3,3) 
        plot(c_AS,(GammaFit_AS_all.sharedfit.GSweights),'k.')
        hold on
        plot(c_AS,mean((GammaFit_AS_all.sharedfit.GSweights),2),'ko-')
        box off
        xlabel('AS Weight Cost');ylabel('GS Weight')
       % LogScale('x',10)
        
    subplot(3,3,4) 
        plot(c_AS,(GammaFit_AS_all.sharedfit.GSlogrates),'k.')
        hold on
        plot(c_AS,mean((GammaFit_AS_all.sharedfit.GSlogrates),2),'ko-')
        box off
        xlabel('AS Weight Cost');ylabel('GS Rate')
       % LogScale('x',10)
        
    subplot(3,3,5)    
        plot(-GammaFit_AS_all.sharedfit.ASlogrates',GammaFit_AS_all.sharedfit.ASCVs','.')
        
    subplot(3,3,7)
        bz_PlotISIDistModes(GammaFit_AS(1),GammaFit_AS(1).cellstats.UID,...
            'dotscale',10,'dotscaleAS',150)
        xlim([-3 2])
        
    subplot(3,3,8)
        bz_PlotISIDistModes(GammaFit_AS(2),GammaFit_AS(2).cellstats.UID,...
            'dotscale',10,'dotscaleAS',150)
        xlim([-3 2])
        
    subplot(3,3,9)
        bz_PlotISIDistModes(GammaFit_AS(5),GammaFit_AS(5).cellstats.UID,...
            'dotscale',10,'dotscaleAS',150)
        xlim([-3 2])

    NiceSave('ASWeightCost',figfolder,baseName);

%%
    weightthresh = 0.01;
    clear modeweightcorr allweights numsigAS
    for pp = 1:6
        allweights{pp} = GammaFit_AS(pp).sharedfit.ASweights;
        allweights{pp}(log10(allweights{pp})<-4) = 1e-4;
        numsigAS{pp} = sum(allweights{pp}>weightthresh,2);

        modeweightcorr{pp} = corr([allweights{pp} GammaFit_AS(pp).sharedfit.GSweights'] ,...
            'type','spearman');
    end
%%
    figure
    for pp = 1:6      
       
        subplot(4,6,pp)
            imagesc(modeweightcorr{pp})
            colorbar
            crameri('vik','pivot',0)
            xlabel('Mode');ylabel('Mode')
            set(gca,'ytick',[0:pp])
            set(gca,'xtick',[0:pp])
            
        subplot(4,6,pp+6)
            hist(log10(allweights{pp}(:)))
            hold on;box off; axis tight
            plot(log10(weightthresh).*[1 1],ylim(gca),'k--')
            LogScale('x',10)
            xlabel('Weight');ylabel('# Modes (All Cells)')

        subplot(4,6,pp+12)
            hist(numsigAS{pp},linspace(0,numAS,numAS+1))
            hold on
            box off
            axis tight
            xlabel('# Modes');ylabel('# Cells') 
    end
    NiceSave('ModeWeightProperties_AS',figfolder,baseName);
end
end

%% Keep only the selected iteration. Change this when move to function

for ss = 1:2
    GammaFit.(statenames{ss}).sharedfit = GammaFit.(statenames{ss}).sharedfit(keepAS);
end
%%
%joint indexing across different sets of intervals
GammaFit.(statenames{1}).cellstats.NW = false(size(GammaFit.(statenames{1}).cellstats.meanrate));
GammaFit.(statenames{2}).cellstats.NW = false(size(GammaFit.(statenames{2}).cellstats.meanrate));
[~,NW1,NW2]=intersect(GammaFit.(statenames{1}).cellstats.UID,GammaFit.(statenames{2}).cellstats.UID);
GammaFit.(statenames{1}).cellstats.NW(NW1)=true;
GammaFit.(statenames{2}).cellstats.NW(NW2)=true;


cellinfofilename = fullfile(basePath,[baseName,'.GammaFit.cellinfo.mat']); %Update in a bit and below
if SAVECELLINFO
    save(cellinfofilename,'GammaFit')
end




GScolor = [0.6 0.4 0];

%%


%%
diffrate = log10(GammaFit.WAKEstate.cellstats.meanrate(GammaFit.WAKEstate.cellstats.NW))-...
    log10(GammaFit.NREMstate.cellstats.meanrate(GammaFit.NREMstate.cellstats.NW));
diffGS = [GammaFit.WAKEstate.singlecell(GammaFit.WAKEstate.cellstats.NW).GSlogrates]-...
     [GammaFit.NREMstate.singlecell(GammaFit.NREMstate.cellstats.NW).GSlogrates];
 diffAR = (1-[GammaFit.WAKEstate.singlecell(GammaFit.WAKEstate.cellstats.NW).GSweights]) - ...
      (1-[GammaFit.NREMstate.singlecell(GammaFit.NREMstate.cellstats.NW).GSweights]);
%%
figure
subplot(3,3,1)
    plot(log10(GammaFit.WAKEstate.cellstats.meanrate(GammaFit.WAKEstate.cellstats.NW)),...
        log10(GammaFit.NREMstate.cellstats.meanrate(GammaFit.NREMstate.cellstats.NW)),'.');
    hold on
    UnityLine
    xlabel('WAKE');ylabel('NREM')
    title('Mean Rate')

subplot(3,3,2)
    plot([GammaFit.WAKEstate.singlecell(GammaFit.WAKEstate.cellstats.NW).GSlogrates],...
        [GammaFit.NREMstate.singlecell(GammaFit.NREMstate.cellstats.NW).GSlogrates],'.');
    hold on
    UnityLine
    xlabel('WAKE ');ylabel('NREM')
    title('GS Rate')

subplot(3,3,3)
    plot(1-[GammaFit.WAKEstate.singlecell(GammaFit.WAKEstate.cellstats.NW).GSweights],...
        1-[GammaFit.NREMstate.singlecell(GammaFit.NREMstate.cellstats.NW).GSweights],'.');
    hold on
    UnityLine
    xlabel('WAKE');ylabel('NREM')
    title('AS Ratio')

subplot(3,3,4)
    plot(diffrate,diffGS,'.')
    hold on
    UnityLine
    plot([0 0],ylim(gca),'k')
    plot(xlim(gca),[0 0],'k')
    xlabel('NW Rate Diff (W-N)');ylabel('NW GSRate Diff (W-N)')
subplot(3,3,5)
    plot(diffrate,diffAR,'.')
    hold on
    UnityLine
    plot([0 0],ylim(gca),'k')
    plot(xlim(gca),[0 0],'k')
    xlabel('NW Rate Diff (W-N)');ylabel('NW ASWeight Diff (W-N)')
subplot(3,3,6)
    plot(diffGS,diffAR,'.')
    hold on
    UnityLine
    plot([0 0],ylim(gca),'k')
    plot(xlim(gca),[0 0],'k')
    xlabel('NW GSRate Diff (W-N)');ylabel('NW ASWeight Diff (W-N)')
    
    
subplot(3,4,9)
    plot(log10(GammaFit.NREMstate.cellstats.meanrate(GammaFit.NREMstate.cellstats.NW)),diffGS,'.')
    hold on
    plot(xlim(gca),[0 0],'k')
    xlabel('N Rate');ylabel('NW GSRate Diff (W-N)')
subplot(3,4,10)
    plot(log10(GammaFit.NREMstate.cellstats.meanrate(GammaFit.NREMstate.cellstats.NW)),diffAR,'.')
    hold on
    plot(xlim(gca),[0 0],'k')
    xlabel('N Rate');ylabel('NW ASWeight Diff (W-N)')
    
subplot(3,4,11)
    plot(log10(GammaFit.WAKEstate.cellstats.meanrate(GammaFit.WAKEstate.cellstats.NW)),diffGS,'.')
    hold on
    plot(xlim(gca),[0 0],'k')
    xlabel('N Rate');ylabel('NW GSRate Diff (W-N)')
subplot(3,4,12)
    plot(log10(GammaFit.WAKEstate.cellstats.meanrate(GammaFit.WAKEstate.cellstats.NW)),diffAR,'.')
    hold on
    plot(xlim(gca),[0 0],'k')
    xlabel('N Rate');ylabel('NW ASWeight Diff (W-N)')
    
        NiceSave('NWDiff',figfolder,baseName);
    
    

        
        

%% Example cell: 3 states
for ff = 1:3
numex=2;
excells = randi(GammaFit.(statenames{ss}).numcells,numex);
figure
for ee = 1:2
for ss = 1:2
    excell = excells(ee);
subplot(6,3,ss+(ee-1)*9)
    plot(GammaFit.(statenames{ss}).logtimebins,...
        GammaFit.(statenames{ss}).ISIdists(:,excell),...
        'color',[0.5 0.5 0.5],'linewidth',2)
    hold on
    plot(GammaFit.(statenames{ss}).logtimebins,...
        GSASmodel(GammaFit.(statenames{ss}).singlecell(excell),...
        GammaFit.(statenames{ss}).taubins),...
        statecolor{ss},'linewidth',1)
    hold on
    plot(GammaFit.(statenames{ss}).logtimebins,...
        LogGamma(GammaFit.(statenames{ss}).singlecell(excell).GSlogrates,...
        GammaFit.(statenames{ss}).singlecell(excell).GSCVs,...
        GammaFit.(statenames{ss}).singlecell(excell).GSweights',...
        GammaFit.(statenames{ss}).taubins'),'color',GScolor,'linewidth',0.25);
    for aa = 1:keepAS
        plot(GammaFit.(statenames{ss}).logtimebins,...
            LogGamma(GammaFit.(statenames{ss}).singlecell(excell).ASlogrates(aa),...
            GammaFit.(statenames{ss}).singlecell(excell).ASCVs(aa),...
            GammaFit.(statenames{ss}).singlecell(excell).ASweights(aa)',...
            GammaFit.(statenames{ss}).taubins'),'k','linewidth',0.25);
    end
    box off
    axis tight
    if ee == 1
        title(statenames{ss})
    end
    if ss == 1
        ylabel(['UID: ',num2str(GammaFit.(statenames{ss}).cellstats.UID(excell))])
    end
    xlim([-3 2])

subplot(6,3,[3 6]+ss+(ee-1)*9)
    scatter(-GammaFit.(statenames{ss}).singlecell(excell).ASlogrates(:),...
        log10(GammaFit.(statenames{ss}).singlecell(excell).ASCVs(:)),...
        100*GammaFit.(statenames{ss}).singlecell(excell).ASweights(:)+0.00001,'k','filled')
    hold on
    scatter(-GammaFit.(statenames{ss}).singlecell(excell).GSlogrates,...
        log10(GammaFit.(statenames{ss}).singlecell(excell).GSCVs),...
        100*GammaFit.(statenames{ss}).singlecell(excell).GSweights+0.00001,GScolor,'filled')
    plot(GammaFit.(statenames{ss}).logtimebins([1 end]),[0 0],'k--')
    ylabel('CV');xlabel('mean ISI (s)')
    xlim([-3 1.9])
    ylim([-2 0.75])
    LogScale('x',10,'exp',true)
    LogScale('y',10)
    box on

end
end
%if figfolder
    NiceSave(['CellExample_states',num2str(ff)],figfolder,baseName);
%end
end
%% Mean and all points (single cell and group)
%Mean dist with group AS. use mean weight.

%% Mean rate and GS rate
% figure
% plot(sharedfit.GSlogrates,log10(ISIStats.summstats.WAKEstate.meanrate(ISIStats.sorts.WAKEstate.ratepE)),'.')
% hold on
% %UnityLine
% xlabel('GS Rate');
% ylabel('Mean Rate')









%%
%cc = 1
% Nmodes = 5;
% maxNmodes = 12;
% numcells = length(ISIStats.summstats.WAKEstate.meanrate);
% clear ISIfits
% for ss = 1:3
%     for cc = 1:numcells
%      
%     bz_Counter(cc,numcells,'Cell')
%     fitISIs = InIntervals(ISIStats.allspikes.times{cc},SleepState.ints.(statenames{ss}));
%     fitISIs = ISIStats.allspikes.ISIs{cc}(fitISIs);
%     [ISIfits.(statenames{ss}).lambdas(:,cc),ISIfits.(statenames{ss}).ks(:,cc),...
%         ISIfits.(statenames{ss}).weights(:,cc),ISIfits.(statenames{ss}).fiterror(cc,:),...
%         ISIfits.(statenames{ss}).Nmodes(cc)] = ...
%         bz_FitISIGammaModes(fitISIs,...
%         'showfig',false,'returnNmodes',Nmodes,'maxNmodes',maxNmodes,...
%         'sequentialreduce',true,'Nestimatemethod','descending');
%     end
% 
%     ISIfits.(statenames{ss}).rates = 1./(ISIfits.(statenames{ss}).ks./ISIfits.(statenames{ss}).lambdas);
%     ISIfits.(statenames{ss}).CVs = 1./ISIfits.(statenames{ss}).ks;
%     ISIfits.(statenames{ss}).weights(ISIfits.(statenames{ss}).weights<0.01) = nan;
%     ISIfits.(statenames{ss}).rates(isnan(ISIfits.(statenames{ss}).weights)) = nan;
%     
% end
% 
% 
% %% Example cell
% 
% cc = randi(numcells);
% for ss =1:3
% %
% fitISIs = InIntervals(ISIStats.allspikes.times{cc},SleepState.ints.(statenames{ss}));
% fitISIs = ISIStats.allspikes.ISIs{cc}(fitISIs);
% [~] = ...
%     bz_FitISIGammaModes(fitISIs,...
%     'showfig',true,'sequentialreduce',true,...
%     'maxNmodes',maxNmodes,'returnNmodes',Nmodes,'autoNmodes','LargeInflection',...
%     'Nestimatemethod','descending');
% 
%     NiceSave(['ISImodefits_ExCell_',num2str(cc),'_',(statenames{ss})],figfolder,baseName)
% end
% %%
% for ss = 1:3
% figure
% subplot(2,2,1)
% hold on
% for cc = 1:2
%     for mm = 1:Nmodes
%     scatter(log10(ISIfits.(statenames{ss}).rates(mm,CellClass.(celltypes{cc}))),...
%         ISIfits.(statenames{ss}).CVs(mm,CellClass.(celltypes{cc})),...
%         10.*ISIfits.(statenames{ss}).weights(mm,CellClass.(celltypes{cc})),cellcolor{cc})
%     end
% %     scatter(log10(1./(ks(:,CellClass.(celltypes{cc}))./rates(:,CellClass.(celltypes{cc})))),...
% %         (1./ks(:,CellClass.(celltypes{cc}))),weights(:,CellClass.(celltypes{cc})),...
% %         repmat(ISIStats.summstats.NREMstate.meanrate(CellClass.(celltypes{cc})),4,1))
%     LogScale('x',10)
%     xlabel('Rate');ylabel('CV')
% end
% title((statenames{ss}))
% 
% subplot(2,2,2)
% %imagesc(log10(fiterror))
% hold on
% for cc = 1:2
%     plot(1:maxNmodes,nanmean(log10(ISIfits.(statenames{ss}).fiterror(CellClass.(celltypes{cc}),:)),1),...
%         '-o','linewidth',2,'color',cellcolor{cc})
%     errorshade(1:maxNmodes,nanmean(log10(ISIfits.(statenames{ss}).fiterror(CellClass.(celltypes{cc}),:)),1),...
%         nanstd(log10(ISIfits.(statenames{ss}).fiterror(CellClass.(celltypes{cc}),:)),[],1),...
%         nanstd(log10(ISIfits.(statenames{ss}).fiterror(CellClass.(celltypes{cc}),:)),[],1),cellcolor{cc},'scalar');
% end
% LogScale('y',10)
% xlabel('N Modes');ylabel('Total Squared Error')
% 
% subplot(4,2,5)
% hold on
% for cc = 1:2
%     plot(log10(ISIfits.(statenames{ss}).rates(:,CellClass.(celltypes{cc}))),...
%         log10(ISIfits.(statenames{ss}).weights(:,CellClass.(celltypes{cc}))),'.','color',cellcolor{cc})
% 
% end
%     LogScale('xy',10)
%     xlabel('Rate');ylabel('Weight')
% 
% subplot(4,2,7)
% hold on
% for cc = 1:2
%     for mm = 1:Nmodes
%     scatter(log10(ISIfits.(statenames{ss}).rates(mm,CellClass.(celltypes{cc}))),...
%         log10(ISIStats.summstats.(statenames{ss}).meanrate(CellClass.(celltypes{cc}))),...
%         10.*ISIfits.(statenames{ss}).weights(mm,CellClass.(celltypes{cc})),cellcolor{cc})
%     end
%     
% %    plot(log10(ISIfits.(statenames{ss}).rates(:,CellClass.(celltypes{cc}))),...
% %        repmat(log10(ISIStats.summstats.(statenames{ss}).meanrate(CellClass.(celltypes{cc}))),Nmodes,1),'.','color',cellcolor{cc})
% 
% end
%     LogScale('xy',10)
%     UnityLine
%     ylim(log10([min(ISIStats.summstats.(statenames{ss}).meanrate) max(ISIStats.summstats.(statenames{ss}).meanrate)]))
%     xlabel('Mode Rate');ylabel(' Cell Rate')
% 
%     for cc = 1:2
% subplot(4,2,6+(cc-1)*2)
%     hist(ISIfits.(statenames{ss}).Nmodes(CellClass.(celltypes{cc})))
%     end
%     
%     NiceSave(['ISImodefits',(statenames{ss})],figfolder,baseName)
% 
% end
%%
% figure
%
% subplot(2,2,1)
% hold on
% for cc = 1
%     plot(log10(rates(:,CellClass.(celltypes{cc}))),...
%         repmat(log10(ISIStats.summstats.WAKEstate.meanrate(CellClass.(celltypes{cc}))),Nmodes,1),'.','color',cellcolor{cc})
% %     scatter(log10(1./(ks(:,CellClass.(celltypes{cc}))./rates(:,CellClass.(celltypes{cc})))),...
% %         (1./ks(:,CellClass.(celltypes{cc}))),weights(:,CellClass.(celltypes{cc})),...
% %         repmat(ISIStats.summstats.NREMstate.meanrate(CellClass.(celltypes{cc})),4,1))
% end
%     LogScale('xy',10)
%     UnityLine
%     %ylim(log10([min(ISIStats.summstats.NREMstate.meanrate) max(ISIStats.summstats.NREMstate.meanrate)]))
%     xlabel('Mode Rate');ylabel(' Cell Rate')
%
% title('CA1 - NREM')
%
% subplot(2,2,2)
% hold on
% for cc = 1
%     for mm = 1:Nmodes
% %     plot(log10(weights(:,CellClass.(celltypes{cc}))),...
% %         repmat(log10(ISIStats.summstats.NREMstate.meanrate(CellClass.(celltypes{cc}))),3,1),'.','color',cellcolor{cc})
%     scatter(log10(rates(mm,CellClass.(celltypes{cc}))),...
%         log10(ISIStats.summstats.WAKEstate.meanrate(CellClass.(celltypes{cc}))),...
%         5,20*weights(mm,CellClass.(celltypes{cc})))
%     end
% end
%     LogScale('xy',10)
%     %UnityLine
%     xlabel('Weight');ylabel('mean rate')

end
