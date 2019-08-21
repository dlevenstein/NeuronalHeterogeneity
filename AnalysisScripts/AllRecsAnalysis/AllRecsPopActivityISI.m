reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/SpikeStatsbyPopActivityAnalysis'];

datasetPath.fCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/BW_CTX';
datasetPath.CA1 = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AG_HPC';
datasetPath.vCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/YS_CTX';
datasetPath.THAL = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AP_THAL';
regions = {'THAL','vCTX','fCTX','CA1'};


popthresh.pE = 15;
popthresh.pI = 5;
popthresh.ALL = 15;

for rr = 1:length(regions)
    [ISIStats.(regions{rr}),baseNames] = bz_LoadCellinfo(datasetPath.(regions{rr}),'ISIStats','dataset',true,'catall',true);
    CellClass.(regions{rr}) = bz_LoadCellinfo(datasetPath.(regions{rr}),'CellClass','dataset',true,'catall',true,'baseNames',baseNames);

    PopActivityAll.(regions{rr}) = GetMatResults(figfolder,'SpikeStatsbyPopActivityAnalysis','baseNames',baseNames);
    PopActivityAll.(regions{rr}) = bz_CollapseStruct(PopActivityAll.(regions{rr}));
    
    recinfo.(regions{rr}).baseName = PopActivityAll.(regions{rr}).name;
    recinfo.(regions{rr}).Ncells = PopActivityAll.(regions{rr}).Ncells;
    recinfo.(regions{rr}).cellinfofiles = baseNames;



    popratehist_joint.(regions{rr}) = bz_CollapseStruct(PopActivityAll.(regions{rr}).popratehist_joint,3,'justcat',true);
    popratehist.(regions{rr}) = bz_CollapseStruct(PopActivityAll.(regions{rr}).popratehist,'match','justcat',true);
    ISIbySynch.(regions{rr}) = bz_CollapseStruct(PopActivityAll.(regions{rr}).ISIbySynch,'match','justcat',true);
    SynchbyISI.(regions{rr}) = bz_CollapseStruct(PopActivityAll.(regions{rr}).SynchbyISI,'match','justcat',true);
    CV2popcorr.(regions{rr}) = bz_CollapseStruct(PopActivityAll.(regions{rr}).CV2popcorr,'match','justcat',true);
    ratepopcorr.(regions{rr}) = bz_CollapseStruct(PopActivityAll.(regions{rr}).ratepopcorr,'match','justcat',true);

    keeprecs = [recinfo.(regions{rr}).Ncells.pE]>popthresh.pE & [recinfo.(regions{rr}).Ncells.pI]>popthresh.pI; 
    if all(~keeprecs)
        disp(['Region ',regions{rr},' has no recordings with enough cells'])
        continue
    else
        popratehist_joint_mean.(regions{rr}) = bz_CollapseStruct(PopActivityAll.(regions{rr}).popratehist_joint(keeprecs),3,'mean',true);
    end
end
%%
statenames = fieldnames(ISIStats.(regions{1}).summstats);
statecolors = {[0 0 0],[0 0 1],[1 0 0]};
numstates = length(statenames);



for rr = 1:length(regions)
    celltypes = fieldnames(ISIbySynch.(regions{rr}).pE.NREMstate.celltypeidx);
    synchtypes = fieldnames(ISIbySynch.(regions{rr}));
    
    popratehist.(regions{rr}).ALL = popratehist.(regions{rr}).pE + popratehist.(regions{rr}).pI;
    for n = 1:length(recinfo.(regions{rr}).Ncells)
    recinfo.(regions{rr}).Ncells(n).ALL = [recinfo.(regions{rr}).Ncells(n).pE]+ [recinfo.(regions{rr}).Ncells(n).pI];
    end
    for ss = 1:3
    for tt = 1:length(celltypes)
        inclass = ISIbySynch.(regions{rr}).pE.NREMstate.celltypeidx.(celltypes{tt});
        popratehist_joint.(regions{rr}).(statenames{ss}).(celltypes{tt}).cellCV2s = nanmean(popratehist_joint.(regions{rr}).(statenames{ss}).cellCV2s(:,:,inclass),3);
        popratehist_joint.(regions{rr}).(statenames{ss}).(celltypes{tt}).pSpk = nanmean(popratehist_joint.(regions{rr}).(statenames{ss}).pSpk(:,:,inclass),3);
       % popratehist_joint.(regions{rr}).(statenames{ss}).(celltypes{tt}).geomeanISIs = nanmean(popratehist_joint.(regions{rr}).(statenames{ss}).geomeanISIs(:,:,inclass),3);
        
        for st = 1:length(synchtypes)
            popratehist.(regions{rr}).enoughpopcells.(synchtypes{st}) = popratehist.(regions{rr}).(synchtypes{st})>popthresh.(synchtypes{st});
            ISIbySynch.(regions{rr}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{tt}) = nanmean(ISIbySynch.(regions{rr}).(synchtypes{st}).(statenames{ss}).pYX(:,:,inclass&popratehist.(regions{rr}).enoughpopcells.(synchtypes{st})),3);
            SynchbyISI.(regions{rr}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{tt}) = nanmean(SynchbyISI.(regions{rr}).(synchtypes{st}).(statenames{ss}).pYX(:,:,inclass&popratehist.(regions{rr}).enoughpopcells.(synchtypes{st})),3);
            
            enoughpopcellsrec = [recinfo.(regions{rr}).Ncells.(synchtypes{st})]>popthresh.(synchtypes{st});
            popratehist_mean.(regions{rr}).(statenames{ss}).(synchtypes{st}) = nanmean(popratehist.(regions{rr}).(statenames{ss}).(synchtypes{st})(enoughpopcellsrec,:),1);
        end
    end
    end
end

%%
clear Noverthresh
clear sortncells
for rr = 1:length(regions)
    for st = 1:length(synchtypes)
        [~,sortncells.(regions{rr}).(synchtypes{st})]=sort([recinfo.(regions{rr}).Ncells(:).(synchtypes{st})]);
        Noverthresh.(regions{rr}).(synchtypes{st}) = sum([recinfo.(regions{rr}).Ncells(:).(synchtypes{st})]>popthresh.(synchtypes{st}));
        frac.(regions{rr}).(synchtypes{st}) = 1-((Noverthresh.(regions{rr}).(synchtypes{st})+0.5)./length([recinfo.(regions{rr}).Ncells(:).(synchtypes{st})]));
    end
end
%%
figure
subplot(4,2,1)
for rr = 1:4
   hold on
        plot(cat(1,recinfo.(regions{rr}).Ncells.pE),cat(1,recinfo.(regions{rr}).Ncells.pI),'.','markersize',10)
        xlabel('# E Cells');ylabel('# I Cells')
        %title(regions{rr})
end
axis tight
plot(popthresh.pE.*[1 1],get(gca,'ylim'),'k')
plot(get(gca,'xlim'),popthresh.pI.*[1 1],'k')
plot(linspace(0,popthresh.ALL,100),popthresh.ALL-linspace(0,popthresh.ALL,100),'k')

legend(regions,'Location','eastoutside')
for st = 1:length(synchtypes)
for rr = 1:length(regions)
    subplot(4,4,rr+4*st)
    imagesc(popratehist.(regions{rr}).bins.(synchtypes{st})(1,:),[0 1],...
        popratehist.(regions{rr}).(statenames{ss}).(synchtypes{st})(sortncells.(regions{rr}).(synchtypes{st}),:))  
    axis xy
    hold on
    %for st = 1:length(synchtypes)
        plot(popratehist.(regions{rr}).bins.(synchtypes{st})(1,:),...
            bz_NormToRange(popratehist_mean.(regions{rr}).(statenames{ss}).(synchtypes{st}),[0 1]))
        plot(popratehist.(regions{rr}).bins.(synchtypes{st})(1,[1 end]),frac.(regions{rr}).(synchtypes{st}).*[1 1],'k')
    %end
    if rr==1
    ylabel((synchtypes{st}))
    end
    caxis([0 2*max(popratehist_mean.(regions{rr}).(statenames{ss}).ALL)])
    crameri bilbao
    if st ==1
        title(regions{rr})
    end
end
end
NiceSave('CellCounts',figfolder,[])

%%
ss = 1
figure
for rr = 1:length(regions)
    subplot(4,4,rr)
    imagesc(popratehist.(regions{rr}).bins.ALL(1,:),[0 0.1],popratehist.(regions{rr}).(statenames{ss}).ALL(sortncells.(regions{rr}),:))  
    axis xy
    hold on
    for st = 1:length(synchtypes)
        plot(popratehist.(regions{rr}).bins.(synchtypes{st})(1,:),popratehist_mean.(regions{rr}).(statenames{ss}).(synchtypes{st}))
    end
    
    subplot(4,4,rr+4)
        imagesc(ISIbySynch.(regions{rr}).ALL.(statenames{ss}).Xbins(1,:,1),[0 1],...
            squeeze(ISIbySynch.(regions{rr}).ALL.(statenames{ss}).pX(1,:,popratehist.(regions{rr}).enoughpopcells.ALL))')
end

%%
for rr = 2:length(regions)
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