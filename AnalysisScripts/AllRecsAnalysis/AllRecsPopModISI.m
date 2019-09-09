reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/PopActivityModulationAnalysis'];

datasetPath.fCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/BW_CTX';
datasetPath.CA1 = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AG_HPC';
datasetPath.vCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/YS_CTX';
datasetPath.THAL = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AP_THAL';
regions = {'THAL','vCTX','fCTX','CA1'};


popthresh.pE = 25;
popthresh.pI = 5;
popthresh.ALL = 25;

for rr = 1:length(regions)
    disp(['Loading ',regions{rr}])
    %[ISIStats.(regions{rr}),baseNames] = bz_LoadCellinfo(datasetPath.(regions{rr}),'ISIStats','dataset',true,'catall',true);
    %CellClass.(regions{rr}) = bz_LoadCellinfo(datasetPath.(regions{rr}),'CellClass','dataset',true,'catall',true,'baseNames',baseNames);

    [PopModAll,baseNames] = bz_LoadAnalysisResults(datasetPath.(regions{rr}),'PopActivityModulationAnalysis','dataset',true);
    %PopActivityAll = GetMatResults(figfolder,'SpikeStatsbyPopActivityAnalysis','baseNames',baseNames);
    PopModAll = bz_CollapseStruct(PopModAll);
    
    recinfo.(regions{rr}).baseName = PopModAll.baseName;
    %recinfo.(regions{rr}).Ncells = PopModAll.Ncells;
    recinfo.(regions{rr}).cellinfofiles = baseNames;


    PopCellCorr.(regions{rr}) = bz_CollapseStruct(PopModAll.PopCellCorr,'match','justcat',true);
    PopMod.(regions{rr}) = bz_CollapseStruct(PopModAll.PopMod,'match','justcat',true);
    cellinfo.(regions{rr}) = bz_CollapseStruct(PopModAll.cellinfo,'match','justcat',true);
end
%%
statenames = {'WAKEstate','NREMstate','REMstate'};
statecolors = {[0 0 0],[0 0 1],[1 0 0]};
numstates = length(statenames);

%%
synchtypes = {'pE','pI','ALL'};
celltypes = {'pE','pI'};
%%
%% Population average modulation
for rr = 1:4
    PopMod.(regions{rr}).Ncells.ALL = PopMod.(regions{rr}).Ncells.pE + PopMod.(regions{rr}).Ncells.pI;
for ss = 1:3
    for st = 1:length(synchtypes)
        keepcells = PopMod.(regions{rr}).Ncells.(synchtypes{st})>popthresh.(synchtypes{st});
        for cc = 1:length(celltypes)
                PopMod.(regions{rr}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc}) =...
                    nanmean(PopMod.(regions{rr}).(synchtypes{st}).(statenames{ss}).allcells(:,:,cellinfo.(regions{rr}).CellClass.(celltypes{cc})&keepcells),3); 
%                PopMod_MTO.(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc}) =...
%                    nanmean(PopMod_MTO.(synchtypes{st}).(statenames{ss}).allcells(:,:,CellClass.(celltypes{cc})),3); 
%                 PopCellCorr.(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc}) =...
%                     nanmean(PopCellCorr.(synchtypes{st}).(statenames{ss}).allcells(:,CellClass.(celltypes{cc})),3); 
        end
    end
end
end
%%
for rr = 1:4
figure
for ss = 1:3
for st = 1:2
for cc = 1:length(celltypes)
    subplot(4,3,ss+(st-1)*3+(cc-1)*6)
        imagesc((PopMod.(regions{rr}).bins.ISIbins(1,:)),log10(PopMod.(regions{rr}).bins.BinSizeBins(1,:)),...
            log10(PopMod.(regions{rr}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc})))
        %colorbar
        
        %caxis([-0.1 0.1])
        %crameri('vik','pivot',1)
      
        %LogScale('y',2)
         LogScale('xy',10,'exp',true);
        axis xy
        
        if cc==1 &st==1
           title(statenames{ss}) 
        end
        if ss == 1
            ylabel({[(synchtypes{st}),' Modulation'],'Bin Size'})
        end
        if st == 2
            xlabel([(celltypes{cc}),' ISI (s)'])
        end
        
end
end 
end
NiceSave(['PopMod_',(regions{rr})],figfolder,[])
end


%% FIgure : ALLMUA modulation
maxR = [0.1 0.25 0.15 0.3];

figure
for rr = 1:4
for ss = 1:3
for st = 3
for cc = 1:length(celltypes)
    subplot(6,4,(ss-1)*4+(cc-1)*12+rr)
        imagesc((PopMod.(regions{rr}).bins.ISIbins(1,:)),log10(PopMod.(regions{rr}).bins.BinSizeBins(1,:)),...
            log10(PopMod.(regions{rr}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc})))
        colorbar
        
        %caxis([-0.1 0.1])
        caxis([-0.05 maxR(rr)])
        %crameri('vik','pivot',0)
      
        %LogScale('y',2)
         LogScale('xy',10,'exp',true);
        axis xy
        
        if cc==1 &st==1
           title(statenames{ss}) 
        end
        if ss == 1
            ylabel({[(synchtypes{st}),' Modulation'],'Bin Size'})
        end
        if st == 2
            xlabel([(celltypes{cc}),' ISI (s)'])
        end
        
end
end 
end
%NiceSave(['PopMod_',(regions{rr})],figfolder,[])
end