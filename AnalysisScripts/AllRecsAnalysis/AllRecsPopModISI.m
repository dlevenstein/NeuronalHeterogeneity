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
    PopMod_MTO.(regions{rr}) = bz_CollapseStruct(PopModAll.PopMod_MTO,'match','justcat',true);
    cellinfo.(regions{rr}) = bz_CollapseStruct(PopModAll.cellinfo,'match','justcat',true);
    MutInfo.(regions{rr}) = bz_CollapseStruct(PopModAll.MutInfo,'match','justcat',true);
end
%%
statenames = {'WAKEstate','NREMstate','REMstate'};
statecolors = {[0 0 0],[0 0 1],[1 0 0]};
numstates = length(statenames);

%%
synchtypes = {'pE','pI','ALL'};
celltypes = {'pE','pI'};
cellcolor = {[0 0 0],[1 0 0]};

%% Population average modulation
for rr = 1:4
    PopMod.(regions{rr}).Ncells.ALL = PopMod.(regions{rr}).Ncells.pE + PopMod.(regions{rr}).Ncells.pI;
for ss = 1:3
    for st = 1:length(synchtypes)
        keepcells = PopMod.(regions{rr}).Ncells.(synchtypes{st})>popthresh.(synchtypes{st});
        for cc = 1:length(celltypes)
            
            %Remove ISI columns that don't have enoughcells
            nspkthresh = 50;
            ncellthresh.pE = 125;
            ncellthresh.pI = 50;
            %sum(ISIbytheta.(regions{rr}).Xhist>nspkthresh,3)
            %sum(ISIbyPSS.(regions{rr}).Xhist>nspkthresh,3)
            ISIcolsunderthresh = sum(PopMod.(regions{rr}).(synchtypes{st}).(statenames{ss}).nISI(1,:,cellinfo.(regions{rr}).CellClass.(celltypes{cc})&keepcells)>nspkthresh,3)<ncellthresh.(celltypes{cc});
            PopMod.(regions{rr}).(synchtypes{st}).(statenames{ss}).allcells(:,ISIcolsunderthresh,cellinfo.(regions{rr}).CellClass.(celltypes{cc}))=nan;
            ISIcolsunderthresh = sum(PopMod_MTO.(regions{rr}).(synchtypes{st}).(statenames{ss}).nISI(1,:,cellinfo.(regions{rr}).CellClass.(celltypes{cc})&keepcells)>nspkthresh,3)<ncellthresh.(celltypes{cc});
            PopMod_MTO.(regions{rr}).(synchtypes{st}).(statenames{ss}).allcells(:,ISIcolsunderthresh,cellinfo.(regions{rr}).CellClass.(celltypes{cc}))=nan;

            
            
                PopMod.(regions{rr}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc}) =...
                    nanmean(PopMod.(regions{rr}).(synchtypes{st}).(statenames{ss}).allcells(:,:,cellinfo.(regions{rr}).CellClass.(celltypes{cc})&keepcells),3); 
                PopMod_MTO.(regions{rr}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc}) =...
                    nanmean(PopMod_MTO.(regions{rr}).(synchtypes{st}).(statenames{ss}).allcells(:,:,cellinfo.(regions{rr}).CellClass.(celltypes{cc})&keepcells),3); 
                PopCellCorr.(regions{rr}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc}) =...
                    nanmean(PopCellCorr.(regions{rr}).(synchtypes{st}).(statenames{ss}).allcells(:,cellinfo.(regions{rr}).CellClass.(celltypes{cc})&keepcells),2); 
                MutInfo.(regions{rr}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc}) = ...
                    nanmean(MutInfo.(regions{rr}).(synchtypes{st}).(statenames{ss}).allcells(:,cellinfo.(regions{rr}).CellClass.(celltypes{cc})&keepcells),2);
                
                ISIdist.(regions{rr}).(statenames{ss}).(celltypes{cc}) = nanmean(PopMod.(regions{rr}).(synchtypes{st}).(statenames{ss}).pISI(1,:,cellinfo.(regions{rr}).CellClass.(celltypes{cc})&keepcells),3);
                MTOISIdist.(regions{rr}).(statenames{ss}).(celltypes{cc}) = nanmean(PopMod_MTO.(regions{rr}).(synchtypes{st}).(statenames{ss}).pISI(1,:,cellinfo.(regions{rr}).CellClass.(celltypes{cc})&keepcells),3);

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
maxR = [0.125 0.2 0.1 0.3];

figure
for rr = 1:4
for ss = 1:3
for st = 3
for cc = 1:length(celltypes)
    subplot(6,4,(ss-1)*4+(cc-1)*12+rr)
        imagesc((PopMod.(regions{rr}).bins.ISIbins(1,:)),log10(PopMod.(regions{rr}).bins.BinSizeBins(1,:)),...
            (PopMod.(regions{rr}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc})))
        %colorbar
        hold on
        plot(PopMod.(regions{rr}).bins.ISIbins(1,:),...
            bz_NormToRange(ISIdist.(regions{rr}).(statenames{ss}).(celltypes{cc}),0.4),...
        'color',statecolors{ss})
        %caxis([-0.1 0.1])
        caxis(10.^[0 maxR(rr)])
        %crameri('vik','pivot',1)
      
        %LogScale('y',2)
        xlim([-2.7 1.7])
         LogScale('xy',10,'exp',true,'nohalf',true);
        axis xy
        
        
        if ss==1 && cc==1
           title(regions{rr})
        end
        if rr == 1
            ylabel({(statenames{ss}),'Bin Size (s)'})
        end
        if ss == 3 | ss==2
            xlabel([(celltypes{cc}),' ISI (s)'])
        end
        
end
end 
end

end
NiceSave('PopMod',figfolder,[])
%% Sort by rate and plot all cells
% for ss = 1:3
%     
% end
% %%
% exbin = 7; %~20ms
% exbin = 14; %~350ms
% rr=3
% figure
% for ss = 1:3
% for st = 1:3
% subplot(3,3,ss+(st-1)*3)
% h = imagesc(squeeze(log10(PopMod.(regions{rr}).(synchtypes{st}).(statenames{ss}).allcells(exbin,:,:)))');
%  set(h,'AlphaData',50*(squeeze(PopMod.(regions{rr}).(synchtypes{st}).(statenames{ss}).pISI(exbin,:,:))'));
% % PopMod.(synchtypes{st}).(statenames{ss}).pISI(bb,:,:)
% %axis xy
% colorbar
% %caxis([0 0.5])
% end
% end
%% FIgure : ALLMUA modulation MTO
maxR = [0.125 0.225 0.12 0.3];

figure
for rr = 1:4
for ss = 1:3
for st = 3
for cc = 1:length(celltypes)
    subplot(6,4,(ss-1)*4+(cc-1)*12+rr)
        imagesc((PopMod_MTO.(regions{rr}).bins.ISIbins(1,:)),log10(PopMod_MTO.(regions{rr}).bins.BinSizeBins(1,:)),...
            log10(PopMod_MTO.(regions{rr}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc})))
        colorbar
        hold on
        plot([0 0],get(gca,'ylim'),'w--')
        plot(PopMod_MTO.(regions{rr}).bins.ISIbins(1,:),...
            bz_NormToRange(MTOISIdist.(regions{rr}).(statenames{ss}).(celltypes{cc}),0.4),...
        'color',statecolors{ss})
        %caxis([-0.1 0.1])
        caxis([0 maxR(rr)])
        %crameri('vik','pivot',0)
      
        %LogScale('y',2)
         LogScale('xy',10,'exp',true);
        axis xy
        
        if ss==1 && cc==1
           title(regions{rr})
        end
        if rr == 1
            ylabel({(statenames{ss}),'Bin Size (s)'})
        end
        if ss == 3
            xlabel([(celltypes{cc}),' ISI (MTONorm)'])
        end
        
end
end 
end
%
end
NiceSave('PopModMTO',figfolder,[])

%%
for rr = 1:4
figure
for ss = 1:3
for st = 1:3

    subplot(3,3,ss+(st-1)*3)
    hold on
    for cc = 1:length(celltypes)
                if all(~(cellinfo.(regions{rr}).CellClass.(celltypes{cc})))
            continue
                end
        keepcells = PopMod.(regions{rr}).Ncells.(synchtypes{st})>popthresh.(synchtypes{st});
        plotcells = cellinfo.(regions{rr}).CellClass.(celltypes{cc})&keepcells;
% plot(log10(MutInfo.(regions{rr}).bins.BinSizeBins),MutInfo.(regions{rr}).(synchtypes{st}).(statenames{ss}).allcells(:,plotcells)',...
%     'color',min(1,cellcolor{cc}+0.5),'linewidth',0.1)

plot(log10(MutInfo.(regions{rr}).bins.BinSizeBins(1,:)),MutInfo.(regions{rr}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc}),...
    'color',cellcolor{cc},'linewidth',2)
    end
%caxis([-1 1])
%ylim([-1 1])
xlim(log10(MutInfo.(regions{rr}).bins.BinSizeBins([1 end])))
xlabel('TimeScale');
LogScale('x',10)

        if ss == 1
            ylabel([(synchtypes{st}),' MI'])
        end
        if cc==1 &st==1
           title(statenames{ss}) 
        end
%crameri('vik','pivot',0)
end
%end
end
end
%NiceSave('ISIMITimeScale',figfolder,baseName)


%%
st = 3;
figure
for rr = 1:4
for ss = 1:3


    subplot(4,4,rr+(ss-1)*4)
    hold on
    for cc = 1:length(celltypes)
                if all(~(cellinfo.(regions{rr}).CellClass.(celltypes{cc})))
            continue
                end
        keepcells = PopMod.(regions{rr}).Ncells.(synchtypes{st})>popthresh.(synchtypes{st});
        plotcells = cellinfo.(regions{rr}).CellClass.(celltypes{cc})&keepcells;
% plot(log10(MutInfo.(regions{rr}).bins.BinSizeBins),MutInfo.(regions{rr}).(synchtypes{st}).(statenames{ss}).allcells(:,plotcells)',...
%     'color',min(1,cellcolor{cc}+0.5),'linewidth',0.1)

plot(log10(MutInfo.(regions{rr}).bins.BinSizeBins(1,:)),MutInfo.(regions{rr}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc}),...
    'color',cellcolor{cc},'linewidth',2)
    end
%caxis([-1 1])
%ylim([-1 1])
xlim(log10(MutInfo.(regions{rr}).bins.BinSizeBins([1 end])))

LogScale('x',10)

        if rr == 1
            ylabel({(statenames{ss}),' MI'})
        end
        if ss ==3
            xlabel('PopRate Bin (s)')
        end
        if ss==1
           title(regions{rr})
        end
%crameri('vik','pivot',0)
end
%end
end
NiceSave('ISIMITimeScale',figfolder,[])

%%
%%
st = 3;
figure
for rr = 1:4
for ss = 1:3


    subplot(4,4,rr+(ss-1)*4)
    hold on
    for cc = 1:length(celltypes)
                if all(~(cellinfo.(regions{rr}).CellClass.(celltypes{cc})))
            continue
                end
        keepcells = PopMod.(regions{rr}).Ncells.(synchtypes{st})>popthresh.(synchtypes{st});
        plotcells = cellinfo.(regions{rr}).CellClass.(celltypes{cc})&keepcells;
% plot(log10(MutInfo.(regions{rr}).bins.BinSizeBins),MutInfo.(regions{rr}).(synchtypes{st}).(statenames{ss}).allcells(:,plotcells)',...
%     'color',min(1,cellcolor{cc}+0.5),'linewidth',0.1)

plot(log10(PopCellCorr.(regions{rr}).bins.BinSizeBins(1,:)),PopCellCorr.(regions{rr}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{cc}),...
    'color',cellcolor{cc},'linewidth',2)
    end
%caxis([-1 1])
%ylim([-1 1])
xlim(log10(PopCellCorr.(regions{rr}).bins.BinSizeBins([1 end])))
xlabel('TimeScale');
LogScale('x',10)

        if rr == 1
            ylabel({(statenames{ss}),' MI'})
        end
        if cc==1 &st==1
           title(statenames{ss}) 
        end
%crameri('vik','pivot',0)
end
%end
end