reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/SpikeStatsbyPopActivityAnalysis'];

datasetPath.fCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/BW_CTX';
datasetPath.CA1 = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AG_HPC';
datasetPath.vCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/YS_CTX';
datasetPath.THAL = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AP_THAL';
regions = {'THAL','vCTX','fCTX','CA1'};

for rr = 1:length(regions)
    [ISIStats.(regions{rr}),baseNames] = bz_LoadCellinfo(datasetPath.(regions{rr}),'ISIStats','dataset',true,'catall',true);
    CellClass.(regions{rr}) = bz_LoadCellinfo(datasetPath.(regions{rr}),'CellClass','dataset',true,'catall',true,'baseNames',baseNames);

    PopActivityAll = GetMatResults(figfolder,'SpikeStatsbyPopActivityAnalysis','baseNames',baseNames);
    PopActivityAll = bz_CollapseStruct(PopActivityAll);
    
    recinfo.(regions{rr}).baseName = PopActivityAll.name;
    recinfo.(regions{rr}).Ncells = PopActivityAll.Ncells;
    recinfo.(regions{rr}).cellinfofiles = baseNames;



    popratehist_joint.(regions{rr}) = bz_CollapseStruct(PopActivityAll.popratehist_joint,3,'justcat',true);
    popratehist.(regions{rr}) = bz_CollapseStruct(PopActivityAll.popratehist,'match','justcat',true);
    ISIbySynch.(regions{rr}) = bz_CollapseStruct(PopActivityAll.ISIbySynch,'match','justcat',true);
    SynchbyISI.(regions{rr}) = bz_CollapseStruct(PopActivityAll.SynchbyISI,'match','justcat',true);
    CV2popcorr.(regions{rr}) = bz_CollapseStruct(PopActivityAll.CV2popcorr,'match','justcat',true);
    ratepopcorr.(regions{rr}) = bz_CollapseStruct(PopActivityAll.ratepopcorr,'match','justcat',true);

    popratehist_joint_mean.(regions{rr}) = bz_CollapseStruct(PopActivityAll.popratehist_joint,3,'mean',true);

end
%%
statenames = fieldnames(ISIStats.(regions{1}).summstats);
statecolors = {[0 0 0],[0 0 1],[1 0 0]};
numstates = length(statenames);

popthresh.pE = 15;
popthresh.pI = 5;
popthresh.ALL = 20;

for rr = 1:length(regions)
    celltypes = fieldnames(ISIbySynch.(regions{rr}).pE.NREMstate.celltypeidx);
    synchtypes = fieldnames(ISIbySynch.(regions{rr}));
    for ss = 1:3
    for tt = 1:length(celltypes)
        inclass = ISIbySynch.(regions{rr}).pE.NREMstate.celltypeidx.(celltypes{tt});
        popratehist_joint.(regions{rr}).(statenames{ss}).(celltypes{tt}).cellCV2s = nanmean(popratehist_joint.(regions{rr}).(statenames{ss}).cellCV2s(:,:,inclass),3);
        popratehist_joint.(regions{rr}).(statenames{ss}).(celltypes{tt}).pSpk = nanmean(popratehist_joint.(regions{rr}).(statenames{ss}).pSpk(:,:,inclass),3);
       % popratehist_joint.(regions{rr}).(statenames{ss}).(celltypes{tt}).geomeanISIs = nanmean(popratehist_joint.(regions{rr}).(statenames{ss}).geomeanISIs(:,:,inclass),3);
        
        for st = 1:length(synchtypes)
            enoughpopcells = popratehist.(regions{rr}).(synchtypes{tt})>popthresh.(synchtypes{tt});
            ISIbySynch.(regions{rr}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{tt}) = nanmean(ISIbySynch.(regions{rr}).(synchtypes{st}).(statenames{ss}).pYX(:,:,inclass&enoughpopcells),3);
            SynchbyISI.(regions{rr}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{tt}) = nanmean(SynchbyISI.(regions{rr}).(synchtypes{st}).(statenames{ss}).pYX(:,:,inclass&enoughpopcells),3);
        end
    end
    end
end


%%
figure

for rr = 1:4
   hold on
        plot(cat(1,recinfo.(regions{rr}).Ncells.pE),cat(1,recinfo.(regions{rr}).Ncells.pI),'.','markersize',10)
        xlabel('# E Cells');ylabel('# I Cells')
        %title(regions{rr})
end
legend(regions,'Location','northwest')
NiceSave('CellCounts',figfolder,[])

%%
for rr = 1:length(regions)
figure
for ss = 1:3
    subplot(3,3,ss)
        h = imagesc(popratehist_joint.(regions{rr}).Ebins(1,:),popratehist_joint.(regions{rr}).Ibins(1,:),...
            popratehist_joint_mean.(regions{rr}).(statenames{ss}).alltime');
        axis xy
        set(h,'AlphaData',~(popratehist_joint_mean.(regions{rr}).(statenames{ss}).alltime'==0));

        title(statenames{ss})
    for tt = 1:length(celltypes)
    subplot(6,6,(ss-1)*2+12+tt)
        h = imagesc(popratehist_joint.(regions{rr}).Ebins(1,:),popratehist_joint.(regions{rr}).Ibins(1,:),log10(popratehist_joint.(regions{rr}).(statenames{ss}).(celltypes{tt}).pSpk)');
        axis xy
        set(h,'AlphaData',~isnan(popratehist_joint.(regions{rr}).(statenames{ss}).(celltypes{tt}).pSpk'));
        colorbar
        %crameri lajolla
        caxis([-0.5 1.75])
        
%     subplot(6,6,(ss-1)*2+18+tt)
%             h = imagesc(popratehist_joint.(regions{rr}).Ebins(1,:),popratehist_joint.(regions{rr}).Ibins(1,:),1./(popratehist_joint.(regions{rr}).(statenames{ss}).(celltypes{tt}).geomeanISIs)');
%         axis xy
%         set(h,'AlphaData',~isnan(popratehist_joint.(regions{rr}).(statenames{ss}).(celltypes{tt}).geomeanISIs'));
%         colorbar
%         %crameri lajolla
%         %caxis([-1 1])
        
    subplot(6,6,(ss-1)*2+24+tt)
        h = imagesc(popratehist_joint.(regions{rr}).Ebins(1,:),popratehist_joint.(regions{rr}).Ibins(1,:),popratehist_joint.(regions{rr}).(statenames{ss}).(celltypes{tt}).cellCV2s');
        axis xy
        set(h,'AlphaData',~isnan(popratehist_joint.(regions{rr}).(statenames{ss}).(celltypes{tt}).cellCV2s'));
        colorbar
        crameri berlin
        caxis([0.7 1.3])
    end
end
NiceSave(['popratehist_joints_',(regions{rr})],figfolder,[])
end
%
%%
for rr = 1:length(regions)
    
figure
for ss = 1:3
for tt = 1:length(celltypes)
for st = 1:2
subplot(4,3,(ss-1)+(tt-1)*3+(st-1)*6+1)
    imagesc(ISIbySynch.(regions{rr}).(synchtypes{st}).(statenames{ss}).Xbins(1,:,1),ISIbySynch.(regions{rr}).(synchtypes{st}).(statenames{ss}).Ybins(1,:,1), ISIbySynch.(regions{rr}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{tt})')
    %hold on
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    %LogScale('y',10)
    ylabel([(celltypes{tt}),' ISI (log(s))']);xlabel([(synchtypes{st}),' Synch'])
    if tt==1 & st == 1
        title(statenames{ss})
    end
    %title((celltypes{tt}))
   % colorbar
    if tt ==1 
        caxis([0 0.02])
    elseif tt==2
         caxis([0 0.03])
    end
    xlim(ISIbySynch.(regions{rr}).(synchtypes{st}).(statenames{ss}).Xbins(1,[1 end],1))
  
end
end
end
end
%NiceSave('ISIbySynch.(regions{rr})',figfolder,baseName)
%%
figure
for ss = 1:3
for tt = 1:length(celltypes)
for st = 1:2
subplot(4,3,(ss-1)+(tt-1)*3+(st-1)*6+1)
    
    imagesc(SynchbyISI.(synchtypes{st}).(statenames{ss}).Xbins(1,:,1),SynchbyISI.(synchtypes{st}).(statenames{ss}).Ybins(1,:,1), SynchbyISI.(synchtypes{st}).(statenames{ss}).pop.(celltypes{tt})')
    %hold on
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    %LogScale('y',10)
    xlabel([(celltypes{tt}),' ISI (log(s))']);ylabel([(synchtypes{st}),' Synch'])
    %title((celltypes{tt}))
    if tt==1 & st == 1
        title(statenames{ss})
    end
%     if tt ==1 
%         caxis([0 0.02])
%     elseif tt==2
%          caxis([0 0.03])
%     end
    xlim(SynchbyISI.(synchtypes{tt}).(statenames{ss}).Xbins(1,[1 end],1))
end 
end
end
%NiceSave('SynchbyISI',figfolder,baseName)

%%

    
figure
for rr = 1:length(regions)
for ss = 1:3
for tt = 1:length(celltypes)
for st = 3
subplot(6,4,(ss-1)*4+(tt-1)*12+rr)
    imagesc(ISIbySynch.(regions{rr}).(synchtypes{st}).(statenames{ss}).Xbins(1,:,1),ISIbySynch.(regions{rr}).(synchtypes{st}).(statenames{ss}).Ybins(1,:,1), ISIbySynch.(regions{rr}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{tt})')
    %hold on
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    %LogScale('y',10)
    if rr==1
    ylabel({(statenames{ss}),[(celltypes{tt}),' ISI (log(s))']});
    end
    if ss==3
        xlabel(['Pop Rate, ',(synchtypes{st}),' cells (Hz/cell)'])
    end
    if tt==1 & ss==1
        title(regions{rr})
    end
    %title((celltypes{tt}))
   % colorbar
    if tt ==1 
        caxis([0 0.02])
    elseif tt==2
         caxis([0 0.03])
    end
    xlim(ISIbySynch.(regions{rr}).(synchtypes{st}).(statenames{ss}).Xbins(1,[1 end],1))
  
end
end
end
end
NiceSave('SynchbyISI',figfolder,[])
%%
cellcolor = {'k','r'};
for ss = 1:3
figure
    subplot(3,3,1)
        for tt = 1:length(celltypes)
            plot(log10(ISIStats.(regions{rr}).summstats.(statenames{ss}).meanrate(CellClass.(regions{rr}).(celltypes{tt}))),...
                ratepopcorr.(regions{rr}).(statenames{ss}).pE(CellClass.(regions{rr}).(celltypes{tt})),...
                '.','color',cellcolor{tt})
            hold on
        end
        plot(get(gca,'xlim'),[0 0],'k')
        xlabel('Mean Rate (Hz)');ylabel('Cell Rate-pE Rate Corr.')
        LogScale('x',10)
        title((statenames{ss}))

    subplot(3,3,2)
        for tt = 1:length(celltypes)
            plot(log10(ISIStats.(regions{rr}).summstats.(statenames{ss}).meanrate(CellClass.(regions{rr}).(celltypes{tt}))),ratepopcorr.(regions{rr}).(statenames{ss}).pI(CellClass.(regions{rr}).(celltypes{tt})),...
                '.','color',cellcolor{tt})
            hold on
        end
        plot(get(gca,'xlim'),[0 0],'k')
        xlabel('Mean Rate (Hz)');ylabel('Cell Rate-pI Rate Corr.')
        LogScale('x',10)
        
        
    subplot(3,3,3)
        for tt = 1:length(celltypes)
            plot(ratepopcorr.(regions{rr}).(statenames{ss}).pE(CellClass.(regions{rr}).(celltypes{tt})),ratepopcorr.(regions{rr}).(statenames{ss}).pI(CellClass.(regions{rr}).(celltypes{tt})),...
                '.','color',cellcolor{tt})
            hold on
        end
        xlabel('Rate-pE Corr');ylabel('Rate-pI Corr')
        axis tight
        plot(get(gca,'xlim'),[0 0],'k')
        plot([0 0],get(gca,'ylim'),'k')
        
        
    subplot(3,3,4)
    for tt = 1:length(celltypes)
        plot(log10(ISIStats.(regions{rr}).summstats.NREMstate.meanrate(CellClass.(regions{rr}).(celltypes{tt}))),CV2popcorr.(regions{rr}).(statenames{ss}).pE(CellClass.(regions{rr}).(celltypes{tt})),...
            '.','color',cellcolor{tt})
        hold on
    end
    plot(get(gca,'xlim'),[0 0],'k')
    xlabel('Mean Rate (Hz)');ylabel('CV2-pE Rate Corr.')
    LogScale('x',10)
    %title(statenames{ss})
    
    subplot(3,3,5)
    for tt = 1:length(celltypes)
        plot(log10(ISIStats.(regions{rr}).summstats.NREMstate.meanrate(CellClass.(regions{rr}).(celltypes{tt}))),CV2popcorr.(regions{rr}).(statenames{ss}).pI(CellClass.(regions{rr}).(celltypes{tt})),...
            '.','color',cellcolor{tt})
        hold on
    end
    plot(get(gca,'xlim'),[0 0],'k')
    xlabel('Mean Rate (Hz)');ylabel('CV2-pI Rate Corr.')
    LogScale('x',10)
    
    
    subplot(3,3,6)
    for tt = 1:length(celltypes)
        plot(CV2popcorr.(regions{rr}).(statenames{ss}).pE(CellClass.(regions{rr}).(celltypes{tt})),CV2popcorr.(regions{rr}).(statenames{ss}).pI(CellClass.(regions{rr}).(celltypes{tt})),...
            '.','color',cellcolor{tt})
        hold on
    end
    plot(get(gca,'xlim'),[0 0],'k')
    plot([0 0],get(gca,'ylim'),'k')
    xlabel('CV2-pE Corr.');ylabel('CV2-pI Corr.')
    
    subplot(3,3,7)
    for tt = 1:length(celltypes)
        plot(ratepopcorr.(regions{rr}).(statenames{ss}).pE(CellClass.(regions{rr}).(celltypes{tt})),CV2popcorr.(regions{rr}).(statenames{ss}).pE(CellClass.(regions{rr}).(celltypes{tt})),...
            '.','color',cellcolor{tt})
        hold on
    end
    plot(get(gca,'xlim'),[0 0],'k')
    plot([0 0],get(gca,'ylim'),'k')
    xlabel('Rate-pE Corr.');ylabel('CV2-pE Corr.')
    
    subplot(3,3,8)
    for tt = 1:length(celltypes)
        plot(ratepopcorr.(regions{rr}).(statenames{ss}).pI(CellClass.(regions{rr}).(celltypes{tt})),CV2popcorr.(regions{rr}).(statenames{ss}).pI(CellClass.(regions{rr}).(celltypes{tt})),...
            '.','color',cellcolor{tt})
        hold on
    end
    plot(get(gca,'xlim'),[0 0],'k')
    plot([0 0],get(gca,'ylim'),'k')
    xlabel('Rate-pI Corr.');ylabel('CV2-pI Corr.')
    
    %NiceSave(['RateCV2PopCorr_',(statenames{ss})],figfolder,baseName)
end