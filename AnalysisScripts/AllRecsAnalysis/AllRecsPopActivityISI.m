reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/SpikeStatsbyPopActivityAnalysis'];

datasetPath.fCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/BW_CTX';
datasetPath.CA1 = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AG_HPC';
datasetPath.vCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/YS_CTX';
datasetPath.THAL = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AP_THAL';
regions = {'THAL','vCTX','fCTX','CA1'};

synchtypes = {'pE','pI','ALL'};
normtypes = {'lin','log','lognorm','norm'};
popthresh.pE = 25;
popthresh.pI = 5;
popthresh.ALL = 25;

for rr = 1:length(regions)
    disp(['Loading ',regions{rr}])
    %[ISIStats.(regions{rr}),baseNames] = bz_LoadCellinfo(datasetPath.(regions{rr}),'ISIStats','dataset',true,'catall',true);
    %CellClass.(regions{rr}) = bz_LoadCellinfo(datasetPath.(regions{rr}),'CellClass','dataset',true,'catall',true,'baseNames',baseNames);

    [PopActivityAll,baseNames] = bz_LoadAnalysisResults(datasetPath.(regions{rr}),'SpikeStatsbyPopActivityAnalysis','dataset',true);
    CellClass.(regions{rr}) = bz_LoadCellinfo(datasetPath.(regions{rr}),'CellClass','dataset',true,'catall',true,'baseNames',baseNames);
    %PopActivityAll = GetMatResults(figfolder,'SpikeStatsbyPopActivityAnalysis','baseNames',baseNames);
    PopActivityAll = bz_CollapseStruct(PopActivityAll);
    
    recinfo.(regions{rr}).baseName = PopActivityAll.baseName;
    recinfo.(regions{rr}).cellinfofiles = baseNames;
    recinfo.(regions{rr}).Ncells.pE = [PopActivityAll.Ncells.pE]';
    recinfo.(regions{rr}).Ncells.pI = [PopActivityAll.Ncells.pI]';
    recinfo.(regions{rr}).Ncells.ALL = recinfo.(regions{rr}).Ncells.pE+recinfo.(regions{rr}).Ncells.pI;

    popratehist_joint.(regions{rr}) = bz_CollapseStruct(PopActivityAll.popratehist_joint,3,'justcat',true);
    popratehist.(regions{rr}) = bz_CollapseStruct(PopActivityAll.popratehist,'match','justcat',true);
    ISIbySynch.(regions{rr}) = bz_CollapseStruct(PopActivityAll.ISIbySynch,'match','justcat',true);
    SynchbyISI.(regions{rr}) = bz_CollapseStruct(PopActivityAll.SynchbyISI,'match','justcat',true);
    normISIbySynch.(regions{rr}) = bz_CollapseStruct(PopActivityAll.normISIbySynch,'match','justcat',true);
    CV2popcorr.(regions{rr}) = bz_CollapseStruct(PopActivityAll.CV2popcorr,'match','justcat',true);
    ratepopcorr.(regions{rr}) = bz_CollapseStruct(PopActivityAll.ratepopcorr,'match','justcat',true);
    
    for st = 3
        keeprecs = recinfo.(regions{rr}).Ncells.(synchtypes{st})>popthresh.(synchtypes{st});
        PopRatebyPSS.(regions{rr}).(synchtypes{st}) = bz_CollapseStruct(PopActivityAll.PopRatebyPSS(keeprecs),3,'mean',true);
        PopRatebyTheta.(regions{rr}).(synchtypes{st}) = bz_CollapseStruct(PopActivityAll.PopRatebyTheta(keeprecs),3,'mean',true);
    end
    
    keeprecs = [recinfo.(regions{rr}).Ncells.pE]>popthresh.pE & [recinfo.(regions{rr}).Ncells.pI]>popthresh.pI; 
    if all(~keeprecs)
        disp(['Region ',regions{rr},' has no recordings with enough cells'])
        continue
    else
        popratehist_joint_mean.(regions{rr}) = bz_CollapseStruct(PopActivityAll.popratehist_joint(keeprecs),3,'mean',true);
    end

    clear PopActivityAll
end
%%
statenames = {'WAKEstate','NREMstate','REMstate'};
statecolors = {[0 0 0],[0 0 1],[1 0 0]};
numstates = length(statenames);



for nn = 1:length(normtypes)

for rr = 1:length(regions)
    celltypes = fieldnames(ISIbySynch.(regions{rr}).norm.pE.NREMstate.celltypeidx);
    %synchtypes = fieldnames(ISIbySynch.(regions{rr}).norm);
    
    popratehist.(regions{rr}).ALL = popratehist.(regions{rr}).pE + popratehist.(regions{rr}).pI;
    for n = 1:length(recinfo.(regions{rr}).Ncells)
    recinfo.(regions{rr}).Ncells(n).ALL = [recinfo.(regions{rr}).Ncells(n).pE]+ [recinfo.(regions{rr}).Ncells(n).pI];
    end
    
    
    
    
    for ss = 1:3
    for tt = 1:length(celltypes)
        inclass = ISIbySynch.(regions{rr}).norm.pE.NREMstate.celltypeidx.(celltypes{tt});
        popratehist_joint.(regions{rr}).(normtypes{nn}).(statenames{ss}).(celltypes{tt}).cellCV2s = nanmean(popratehist_joint.(regions{rr}).(normtypes{nn}).(statenames{ss}).cellCV2s(:,:,inclass),3);
        popratehist_joint.(regions{rr}).(normtypes{nn}).(statenames{ss}).(celltypes{tt}).pSpk = nanmean(popratehist_joint.(regions{rr}).(normtypes{nn}).(statenames{ss}).pSpk(:,:,inclass),3);
       % popratehist_joint.(regions{rr}).(normtypes{nn}).(statenames{ss}).(celltypes{tt}).geomeanISIs = nanmean(popratehist_joint.(regions{rr}).(normtypes{nn}).(statenames{ss}).geomeanISIs(:,:,inclass),3);
        
        for st = 1:length(synchtypes)
            
           

            
            popratehist.(regions{rr}).enoughpopcells.(synchtypes{st}) = popratehist.(regions{rr}).(synchtypes{st})>popthresh.(synchtypes{st});
            if nn>2
            %How many cells are contributing?
            nspkthresh = 100;
            ncellthresh = 250;
            %sum(ISIbytheta.(regions{rr}).Xhist>nspkthresh,3)
            %sum(ISIbyPSS.(regions{rr}).Xhist>nspkthresh,3)
            ISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pYX(sum(ISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).Xhist>nspkthresh,3)<ncellthresh,:,:)=nan;
            normISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pYX(sum(normISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).Xhist>nspkthresh,3)<ncellthresh,:,:)=nan;
            
            
            ISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{tt}) = nanmean(ISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pYX(:,:,inclass&popratehist.(regions{rr}).enoughpopcells.(synchtypes{st})),3);
            ISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pop_pX.(celltypes{tt}) = nanmean(ISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pX(:,:,inclass&popratehist.(regions{rr}).enoughpopcells.(synchtypes{st})),3);
            SynchbyISI.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{tt}) = nanmean(SynchbyISI.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pYX(:,:,inclass&popratehist.(regions{rr}).enoughpopcells.(synchtypes{st})),3);
            
            normISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{tt}) = ...
                nanmean(normISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pYX(:,:,inclass&popratehist.(regions{rr}).enoughpopcells.(synchtypes{st})),3);
            normISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pop_pX.(celltypes{tt}) = ...
                nanmean(normISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pX(:,:,inclass&popratehist.(regions{rr}).enoughpopcells.(synchtypes{st})),3);
            end     
            enoughpopcellsrec = [recinfo.(regions{rr}).Ncells.(synchtypes{st})]>popthresh.(synchtypes{st});
            popratehist_mean.(regions{rr}).(normtypes{nn}).(statenames{ss}).(synchtypes{st}) = nanmean(popratehist.(regions{rr}).(normtypes{nn}).(statenames{ss}).(synchtypes{st})(enoughpopcellsrec,:),1);
            popratehist_std.(regions{rr}).(normtypes{nn}).(statenames{ss}).(synchtypes{st}) = nanstd(popratehist.(regions{rr}).(normtypes{nn}).(statenames{ss}).(synchtypes{st})(enoughpopcellsrec,:),[],1);
            
        end
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
%% Pop Rate between regions
display('Making figure - Regional Pop Rate')
for nn = 1:4
figure
for rr = 1:4
for st = 1:length(synchtypes)
    subplot(4,3,st+(rr-1)*3)
        hold on
    for ss = 3:-1:1
         plot(popratehist.(regions{rr}).(normtypes{nn}).bins.(synchtypes{st})(1,:),...
                popratehist_mean.(regions{rr}).(normtypes{nn}).(statenames{ss}).(synchtypes{st}),...
                'color',statecolors{ss},'linewidth',2)

    end
    %legend(regions)
    axis tight 
    xlabel(['Rate (',(normtypes{nn}),')'])
    set(gca,'ytick',[])
    if st==1
       ylabel({(regions{rr}),'P[time]'}) 
    end
    if rr ==1
       title((synchtypes{st}))
    end
end 

end

NiceSave(['PopRateDist_',(normtypes{nn})],figfolder,[])


end
%% EI poprate
display('Making figure - Regional E/I Pop Rate')

for nn = 1:4
figure
    for rr = 2:length(regions)
        for ss = 1:3
    subplot(3,3,ss+(rr-2)*3)
        h = imagesc(popratehist_joint.(regions{rr}).(normtypes{nn}).bins.pE(1,:),...
            popratehist_joint.(regions{rr}).(normtypes{nn}).bins.pI(1,:),...
            popratehist_joint_mean.(regions{rr}).(normtypes{nn}).(statenames{ss}).alltime');
        axis xy
        set(h,'AlphaData',~(popratehist_joint_mean.(regions{rr}).(normtypes{nn}).(statenames{ss}).alltime'==0));
        if ss==1
            ylabel({(regions{rr}),'I Rate'})
        end
        if rr==2
            title(statenames{ss})
        end
        if rr==4
           xlabel('E Rate') 
        end
        
        non0 = popratehist_joint_mean.(regions{rr}).(normtypes{nn}).(statenames{ss}).alltime([2:end],[2:end]);
        caxis([0 max(non0(:))])
        
        end
    end
    
NiceSave(['PopRateDistEI_',(normtypes{nn})],figfolder,[])

    
end
%% Variability between recordings
display('Making figure - Recording Variation')

for nn = 1:4
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
    imagesc(popratehist.(regions{rr}).(normtypes{nn}).bins.(synchtypes{st})(1,:),[0 1],...
        popratehist.(regions{rr}).(normtypes{nn}).(statenames{ss}).(synchtypes{st})(sortncells.(regions{rr}).(synchtypes{st}),:))  
    axis xy
    hold on
    %for st = 1:length(synchtypes)
        plot(popratehist.(regions{rr}).(normtypes{nn}).bins.(synchtypes{st})(1,:),...
            bz_NormToRange(popratehist_mean.(regions{rr}).(normtypes{nn}).(statenames{ss}).(synchtypes{st}),[frac.(regions{rr}).(synchtypes{st}) 1]))
        plot(popratehist.(regions{rr}).(normtypes{nn}).bins.(synchtypes{st})(1,[1 end]),frac.(regions{rr}).(synchtypes{st}).*[1 1],'k')
    %end
    if rr==1
    ylabel({[(synchtypes{st}),' Cells'],'Recording'})
    end
    caxis([0 2*max(popratehist_mean.(regions{rr}).(normtypes{nn}).(statenames{ss}).ALL)])
    crameri bilbao
    set(gca,'ytick',[])
    if st ==1
        title(regions{rr})
    end
    if st==3
        xlabel(['Rate (',normtypes{nn},')'])
    end
end
end
NiceSave(['CellCounts_',(normtypes{nn})],figfolder,[])
end

%%
figure
for rr = 1:4
subplot(4,4,rr)
imagesc(PopRatebyPSS.(regions{rr}).(synchtypes{st}).(normtypes{nn}).(synchtypes{st}).Xbins,...
    PopRatebyPSS.(regions{rr}).(synchtypes{st}).(normtypes{nn}).(synchtypes{st}).Ybins,...
    PopRatebyPSS.(regions{rr}).(synchtypes{st}).(normtypes{nn}).(synchtypes{st}).pYX')
hold on
plot(PopRatebyPSS.(regions{rr}).(synchtypes{st}).(normtypes{nn}).(synchtypes{st}).Xbins,...
    bz_NormToRange(PopRatebyPSS.(regions{rr}).(synchtypes{st}).(normtypes{nn}).(synchtypes{st}).pX),'k')
axis xy
xlabel('PSS');ylabel('Pop Rate (norm)')
title(regions{rr})

subplot(4,4,rr+4)
imagesc(PopRatebyTheta.(regions{rr}).(synchtypes{st}).notNREM.(normtypes{nn}).(synchtypes{st}).Xbins,...
    PopRatebyTheta.(regions{rr}).(synchtypes{st}).notNREM.(normtypes{nn}).(synchtypes{st}).Ybins,...
    PopRatebyTheta.(regions{rr}).(synchtypes{st}).notNREM.(normtypes{nn}).(synchtypes{st}).pYX')
hold on
plot(PopRatebyTheta.(regions{rr}).(synchtypes{st}).notNREM.(normtypes{nn}).(synchtypes{st}).Xbins,...
    bz_NormToRange(PopRatebyTheta.(regions{rr}).(synchtypes{st}).notNREM.(normtypes{nn}).(synchtypes{st}).pX),'k')
axis xy
xlabel('Theta');ylabel('Pop Rate (norm)')
end
NiceSave(['PopDistbyState',(normtypes{nn})],figfolder,[])

%%
figure
subplot(4,2,6)
imagesc(PopRatebyTheta.(regions{rr}).(synchtypes{st}).WAKE.(normtypes{nn}).(synchtypes{st}).Xbins,...
    PopRatebyTheta.(regions{rr}).(synchtypes{st}).WAKE.(normtypes{nn}).(synchtypes{st}).Ybins,...
    PopRatebyTheta.(regions{rr}).(synchtypes{st}).WAKE.(normtypes{nn}).(synchtypes{st}).pYX')
hold on
plot(PopRatebyTheta.(regions{rr}).(synchtypes{st}).WAKE.(normtypes{nn}).(synchtypes{st}).Xbins,...
    bz_NormToRange(PopRatebyTheta.(regions{rr}).(synchtypes{st}).WAKE.(normtypes{nn}).(synchtypes{st}).pX),'k')
axis xy
ylabel({'WAKE','Pop Rate'})

subplot(4,2,8)
imagesc(PopRatebyTheta.(regions{rr}).(synchtypes{st}).REM.(normtypes{nn}).(synchtypes{st}).Xbins,...
    PopRatebyTheta.(regions{rr}).(synchtypes{st}).REM.(normtypes{nn}).(synchtypes{st}).Ybins,...
    PopRatebyTheta.(regions{rr}).(synchtypes{st}).REM.(normtypes{nn}).(synchtypes{st}).pYX')
hold on
plot(PopRatebyTheta.(regions{rr}).(synchtypes{st}).REM.(normtypes{nn}).(synchtypes{st}).Xbins,...
    bz_NormToRange(PopRatebyTheta.(regions{rr}).(synchtypes{st}).REM.(normtypes{nn}).(synchtypes{st}).pX),'k')
axis xy
xlabel('Theta');ylabel({'REM','Pop Rate'})




%%
%Pop Rate Distirbutioon 
%Other cell pop rate distirbution
%Other cell pop rade distribution | spike
for nn = 3:4
figure
for rr = 1:4
for st = 1:length(synchtypes)
    subplot(4,3,st+(rr-1)*3)
        hold on
    for ss = 3:-1:1
        plot(ISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).Xbins(1,:,1),...
        ISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pop_pX.(celltypes{tt}),...
                'color',statecolors{ss},'linewidth',2)

    end
    %legend(regions)
    axis tight 
    xlabel(['Rate (',(normtypes{nn}),')'])
    set(gca,'ytick',[])
    if st==1
       ylabel({(regions{rr}),'P[time]'}) 
    end
    if rr ==1
       title((synchtypes{st}))
    end
    
end 

end

%NiceSave(['PopRateDist_',(normtypes{nn})],figfolder,[])


end
%hold on
%          plot(popratehist.(regions{rr}).(normtypes{nn}).bins.(synchtypes{st})(1,:),...
%                 popratehist_mean.(regions{rr}).(normtypes{nn}).(statenames{ss}).(synchtypes{st}),...
%                 'color',statecolors{ss},'linewidth',2)
%plot(ISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).Xbins(1,:,1),ISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pop_pX.(celltypes{tt}))


%%
display('Making figure - Gross ISI stats byRegional Pop Rate')
for nn = 3:4

figure
for ss = 1:3
for rr = 2:length(regions)
    for tt = 1:length(celltypes)
    subplot(6,6,(ss-1)*6+(tt-1)*18+rr-1)
        h = imagesc(popratehist_joint.(regions{rr}).(normtypes{nn}).bins.pE(1,:),...
            popratehist_joint.(regions{rr}).(normtypes{nn}).bins.pI(1,:),...
            log10(popratehist_joint.(regions{rr}).(normtypes{nn}).(statenames{ss}).(celltypes{tt}).pSpk)');
        axis xy
        set(h,'AlphaData',~isnan(popratehist_joint.(regions{rr}).(normtypes{nn}).(statenames{ss}).(celltypes{tt}).pSpk'));
        %colorbar
        %crameri lajolla
        if tt == 1
            caxis([-0.5 1])
        elseif tt == 2
            caxis([0 1.75])
        end
        if ss ==1 & tt ==1
            title(regions{rr})
        end
        if ss==3&tt==2
           xlabel('pE Rate') 
        end
        if rr == 2
           ylabel('pI Rate') 
        end
%     subplot(6,6,(ss-1)*2+18+tt)
%             h = imagesc(popratehist_joint.(regions{rr}).(normtypes{nn}).bins.pE(1,:),popratehist_joint.(regions{rr}).(normtypes{nn}).bins.pI(1,:),1./(popratehist_joint.(regions{rr}).(normtypes{nn}).(statenames{ss}).(celltypes{tt}).geomeanISIs)');
%         axis xy
%         set(h,'AlphaData',~isnan(popratehist_joint.(regions{rr}).(normtypes{nn}).(statenames{ss}).(celltypes{tt}).geomeanISIs'));
%         colorbar
%         %crameri lajolla
%         %caxis([-1 1])
        
    subplot(6,6,(ss-1)*6+(tt-1)*18++rr+2)
        h = imagesc(popratehist_joint.(regions{rr}).(normtypes{nn}).bins.pE(1,:),popratehist_joint.(regions{rr}).(normtypes{nn}).bins.pI(1,:),popratehist_joint.(regions{rr}).(normtypes{nn}).(statenames{ss}).(celltypes{tt}).cellCV2s');
        axis xy
        set(h,'AlphaData',~isnan(popratehist_joint.(regions{rr}).(normtypes{nn}).(statenames{ss}).(celltypes{tt}).cellCV2s'));
        %colorbar
        crameri berlin
        caxis([0.7 1.3])
        if ss ==1 & tt ==1
            title(regions{rr})
        end
        if ss==3&tt==2
           xlabel('pE Rate') 
        end
    end
end
end
NiceSave(['popratehist_jointstats_',(normtypes{nn}),'_',(regions{rr})],figfolder,[])

end
%
%%
display('Making figure - ISI by Pop Rate')
for nn = 3:4
%nn=3;
for rr = 1:length(regions)
    
figure
for ss = 1:3
for tt = 1:length(celltypes)
for st = 1:2
subplot(4,3,(ss-1)+(tt-1)*6+(st-1)*3+1)
    imagesc(ISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).Xbins(1,:,1),ISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).Ybins(1,:,1), ISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{tt})')
    hold on
	plot(ISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).Xbins(1,:,1),...
        bz_NormToRange(ISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pop_pX.(celltypes{tt})),...
                'color',statecolors{ss},'linewidth',1)
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    %LogScale('y',10)
    ylabel([(celltypes{tt}),' ISI (log(s))']);xlabel([(synchtypes{st}),' Rate (norm)'])
    if tt==1 & st == 1
        title(statenames{ss})
    end
    
    %title((celltypes{tt}))
   % colorbar
%     if tt ==1 
%         caxis([0 0.02])
%     elseif tt==2
%          caxis([0 0.03])
%     end
    %xlim(ISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).Xbins(1,[1 end],1))
  
end
end
end
NiceSave(['ISIbySynch_',(normtypes{nn}),'_',(regions{rr})],figfolder,[])
end
end

%%
display('Making figure - normISI by Pop Rate')

for nn = 3:4
%nn=3;
for rr = 1:length(regions)
    
figure
for ss = 1:3
for tt = 1:length(celltypes)
for st = 1:2
subplot(4,3,(ss-1)+(tt-1)*3+(st-1)*6+1)
    imagesc(normISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).Xbins(1,:,1),...
        normISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).Ybins(1,:,1),...
        normISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{tt})')
    hold on
	plot(normISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).Xbins(1,:,1),...
        bz_NormToRange(normISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pop_pX.(celltypes{tt})),...
                'color',statecolors{ss},'linewidth',1)
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    %LogScale('y',10)
    ylabel([(celltypes{tt}),' ISI (MTOnorm)']);xlabel([(synchtypes{st}),' Rate (norm)'])
    if tt==1 & st == 1
        title(statenames{ss})
    end
    
    %title((celltypes{tt}))
   % colorbar
%     if tt ==1 
%         caxis([0 0.03])
%     elseif tt==2
%          caxis([0 0.04])
%     end
    %xlim(ISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).Xbins(1,[1 end],1))
  
end
end
end
NiceSave(['normISIbySynch_',(normtypes{nn}),'_',(regions{rr})],figfolder,[])
end
end
%%
display('Making figure - Pop Rate by ISI')

figure
for ss = 1:3
for tt = 1:length(celltypes)
for st = 1:2
subplot(4,3,(ss-1)+(tt-1)*3+(st-1)*6+1)
    
    imagesc(SynchbyISI.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).Xbins(1,:,1),...
        SynchbyISI.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).Ybins(1,:,1),...
        SynchbyISI.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{tt})')
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
    xlim(SynchbyISI.(regions{rr}).(normtypes{nn}).(synchtypes{tt}).(statenames{ss}).Xbins(1,[1 end],1))
end 
end
end
%NiceSave('SynchbyISI',figfolder,baseName)

%%
display('Making figure - ISI by Pop Rate (ALL)')

for nn = 3:4
    
figure
for rr = 1:length(regions)
for ss = 1:3
for tt = 1:length(celltypes)
for st = 3
subplot(6,4,(ss-1)*4+(tt-1)*12+rr)
    imagesc(ISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).Xbins(1,:,1),ISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).Ybins(1,:,1), ISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{tt})')
    hold on
	plot(ISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).Xbins(1,:,1),...
        bz_NormToRange(ISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pop_pX.(celltypes{tt})),...
                'color',statecolors{ss},'linewidth',1)
    axis xy
    %LogScale('y',10)
    if rr==1
    ylabel({(statenames{ss}),[(celltypes{tt}),' ISI (log(s))']});
    end
    if ss==3
        xlabel(['Pop Rate, ',(synchtypes{st}),' cells (Norm)'])
    end
    if tt==1 & ss==1
        title(regions{rr})
    end
    %title((celltypes{tt}))
   % colorbar
%     if tt ==1 
%         caxis([0 0.025])
%     elseif tt==2
%          caxis([0 0.03])
%     end
   % xlim(ISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).Xbins(1,[1 end],1))
  
end
end
end
end
NiceSave(['ISIbySynch_',(normtypes{nn})],figfolder,[])
end

%%
display('Making figure - normISI by Pop Rate (ALL)')

for nn = 3:4
    
figure
for rr = 1:length(regions)
for ss = 1:3
for tt = 1:length(celltypes)
for st = 3
subplot(6,4,(ss-1)*4+(tt-1)*12+rr)
    imagesc(normISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).Xbins(1,:,1),...
        normISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).Ybins(1,:,1),...
        normISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pop.(celltypes{tt})')
    hold on
	plot(normISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).Xbins(1,:,1),...
        bz_NormToRange(normISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).pop_pX.(celltypes{tt})),...
                'color',statecolors{ss},'linewidth',1)
    axis xy
    %LogScale('y',10)
    if rr==1
    ylabel({(statenames{ss}),[(celltypes{tt}),' ISI (MTOnorm)']});
    end
    if ss==3
        xlabel(['Pop Rate, ',(synchtypes{st}),' cells (Norm)'])
    end
    if tt==1 & ss==1
        title(regions{rr})
    end
    %title((celltypes{tt}))
   % colorbar
%     if tt ==1 
%         caxis([0 0.03])
%     elseif tt==2
%          caxis([0 0.04])
%     end
    xlim(normISIbySynch.(regions{rr}).(normtypes{nn}).(synchtypes{st}).(statenames{ss}).Xbins(1,[1 end],1))
  
end
end
end
end
NiceSave(['normISIbySynch_',(normtypes{nn})],figfolder,[])
end
%%
display('Making figure - PopRate corr stats')

cellcolor = {'k','r'};
for rr = 1:4
for ss = 1:3
figure
    subplot(3,3,1)
        for tt = 1:length(celltypes)
            plotcells = CellClass.(regions{rr}).(celltypes{tt})&popratehist.(regions{rr}).enoughpopcells.(synchtypes{st});
            plot(log10(ratepopcorr.(regions{rr}).(statenames{ss}).cellrate(plotcells)),...
                ratepopcorr.(regions{rr}).(statenames{ss}).pE(plotcells),...
                '.','color',cellcolor{tt})
            hold on
        end
        axis tight;box off
        plot(get(gca,'xlim'),[0 0],'k')
        xlabel('Mean Rate (Hz)');ylabel('Cell Rate-pE Rate Corr.')
        LogScale('x',10)
        title((regions{rr}))

    subplot(3,3,2)
        for tt = 1:length(celltypes)
            plotcells = CellClass.(regions{rr}).(celltypes{tt})&popratehist.(regions{rr}).enoughpopcells.(synchtypes{st});
            plot(log10(ratepopcorr.(regions{rr}).(statenames{ss}).cellrate(plotcells)),...
                ratepopcorr.(regions{rr}).(statenames{ss}).pI(plotcells),...
                '.','color',cellcolor{tt})
            hold on
        end
        axis tight;box off
        plot(get(gca,'xlim'),[0 0],'k')
        xlabel('Mean Rate (Hz)');ylabel('Cell Rate-pI Rate Corr.')
        LogScale('x',10)
        title((statenames{ss}))
        
    subplot(3,3,3)
        for tt = 1:length(celltypes)
            plotcells = CellClass.(regions{rr}).(celltypes{tt})&popratehist.(regions{rr}).enoughpopcells.(synchtypes{st});
            plot(ratepopcorr.(regions{rr}).(statenames{ss}).pE(plotcells),...
                ratepopcorr.(regions{rr}).(statenames{ss}).pI(plotcells),...
                '.','color',cellcolor{tt})
            hold on
        end
        xlabel('Rate-pE Corr');ylabel('Rate-pI Corr')
        axis tight;box off
        plot(get(gca,'xlim'),[0 0],'k')
        plot([0 0],get(gca,'ylim'),'k')
        
        
    subplot(3,3,4)
    for tt = 1:length(celltypes)
        plotcells = CellClass.(regions{rr}).(celltypes{tt})&popratehist.(regions{rr}).enoughpopcells.(synchtypes{st});
        plot(log10(ratepopcorr.(regions{rr}).(statenames{ss}).cellrate(plotcells)),...
            CV2popcorr.(regions{rr}).(statenames{ss}).pE(plotcells),...
            '.','color',cellcolor{tt})
        hold on
    end
    axis tight;box off
    plot(get(gca,'xlim'),[0 0],'k')
    xlabel('Mean Rate (Hz)');ylabel('CV2-pE Rate Corr.')
    LogScale('x',10)
    %title(statenames{ss})
    
    subplot(3,3,5)
    for tt = 1:length(celltypes)
        plotcells = CellClass.(regions{rr}).(celltypes{tt})&popratehist.(regions{rr}).enoughpopcells.(synchtypes{st});
        plot(log10(ratepopcorr.(regions{rr}).(statenames{ss}).cellrate(plotcells)),...
            CV2popcorr.(regions{rr}).(statenames{ss}).pI(plotcells),...
            '.','color',cellcolor{tt})
        hold on
    end
    axis tight;box off
    plot(get(gca,'xlim'),[0 0],'k')
    xlabel('Mean Rate (Hz)');ylabel('CV2-pI Rate Corr.')
    LogScale('x',10)
    
    
    subplot(3,3,6)
    for tt = 1:length(celltypes)
        plotcells = CellClass.(regions{rr}).(celltypes{tt})&popratehist.(regions{rr}).enoughpopcells.(synchtypes{st});
        plot(CV2popcorr.(regions{rr}).(statenames{ss}).pE(plotcells),...
            CV2popcorr.(regions{rr}).(statenames{ss}).pI(plotcells),...
            '.','color',cellcolor{tt})
        hold on
    end
    plot(get(gca,'xlim'),[0 0],'k')
    plot([0 0],get(gca,'ylim'),'k')
    xlabel('CV2-pE Corr.');ylabel('CV2-pI Corr.')
    
    subplot(3,3,7)
    for tt = 1:length(celltypes)
        plotcells = CellClass.(regions{rr}).(celltypes{tt})&popratehist.(regions{rr}).enoughpopcells.(synchtypes{st});
        plot(ratepopcorr.(regions{rr}).(statenames{ss}).pE(plotcells),...
            CV2popcorr.(regions{rr}).(statenames{ss}).pE(plotcells),...
            '.','color',cellcolor{tt})
        hold on
    end
    plot(get(gca,'xlim'),[0 0],'k')
    plot([0 0],get(gca,'ylim'),'k')
    xlabel('Rate-pE Corr.');ylabel('CV2-pE Corr.')
    
    subplot(3,3,8)
    for tt = 1:length(celltypes)
        plotcells = CellClass.(regions{rr}).(celltypes{tt})&popratehist.(regions{rr}).enoughpopcells.(synchtypes{st});
        plot(ratepopcorr.(regions{rr}).(statenames{ss}).pI(plotcells),...
            CV2popcorr.(regions{rr}).(statenames{ss}).pI(plotcells),...
            '.','color',cellcolor{tt})
        hold on
    end
    plot(get(gca,'xlim'),[0 0],'k')
    plot([0 0],get(gca,'ylim'),'k')
    xlabel('Rate-pI Corr.');ylabel('CV2-pI Corr.')
    
    NiceSave(['RateCV2PopCorr_',(regions{rr}),'_',(statenames{ss})],figfolder,[])
end
end
%%
st =3; %all MUA
figure
for rr = 1:4
    for ss = 1:3
    subplot(3,4,rr+(ss-1)*4)
        for tt = 1:length(celltypes)
            plotcells = CellClass.(regions{rr}).(celltypes{tt})&popratehist.(regions{rr}).enoughpopcells.(synchtypes{st});
            plot(log10(ratepopcorr.(regions{rr}).(statenames{ss}).cellrate(plotcells)),...
                ratepopcorr.(regions{rr}).(statenames{ss}).(synchtypes{st})(plotcells),...
                '.','color',cellcolor{tt})
            hold on
        end
        axis tight;box off
        xlim([-2 2]);ylim([-0.6 0.75])
        plot(get(gca,'xlim'),[0 0],'k')
        if ss ==3
        xlabel('Mean Rate (Hz)');
        end
        LogScale('x',10)
        if ss == 1
        title((regions{rr}))
        end
        if rr == 1
            ylabel({statenames{ss},'Pop-Rate Corr'})
        end
    end
end
NiceSave('RatePopCorr',figfolder,[])
%%
st =3; %all MUA
figure
for rr = 1:4
    for ss = 1:3
    subplot(3,4,rr+(ss-1)*4)
        for tt = 1:length(celltypes)
            plotcells = CellClass.(regions{rr}).(celltypes{tt})&popratehist.(regions{rr}).enoughpopcells.(synchtypes{st});
            plot(log10(ratepopcorr.(regions{rr}).(statenames{ss}).cellrate(plotcells)),...
                CV2popcorr.(regions{rr}).(statenames{ss}).(synchtypes{st})(plotcells),...
                '.','color',cellcolor{tt})
            hold on
        end
        axis tight;box off
        xlim([-2 2]);ylim([-0.25 0.25])
        plot(get(gca,'xlim'),[0 0],'k')
        if ss ==3
        xlabel('Mean Rate (Hz)');
        end
        LogScale('x',10)
        if ss == 1
        title((regions{rr}))
        end
        if rr == 1
            ylabel({statenames{ss},'Pop-CV2 Corr'})
        end
    end
end
NiceSave('CV2PopCorr',figfolder,[])
%% 
for rr = 1:4
    for ss = 1:3
        for st = 1:3
            for tt = 1:length(celltypes)
[~,choristersort.(regions{rr}).(celltypes{tt}).(statenames{ss}).(synchtypes{st})] = sort(ratepopcorr.(regions{rr}).(statenames{ss}).(synchtypes{st}));

usecells = find(CellClass.(regions{rr}).(celltypes{tt})&popratehist.(regions{rr}).enoughpopcells.(synchtypes{st}));
choristersort.(regions{rr}).(celltypes{tt}).(statenames{ss}).(synchtypes{st}) = ...
    intersect(choristersort.(regions{rr}).(celltypes{tt}).(statenames{ss}).(synchtypes{st}),usecells,...
    'stable');
            end
        end
        
        
    end
end

%%
nn= 3
ss = 1
figure
for rr = 1:length(regions)
    for ss = 1:3
    for tt = 1:2
    subplot(6,4,rr+(ss-1)*4+(tt-1)*12)
        imagesc(ISIbySynch.(regions{rr}).(normtypes{nn}).ALL.(statenames{ss}).Xbins(1,:,1),[0 1],...
            squeeze(ISIbySynch.(regions{rr}).(normtypes{nn}).ALL.(statenames{ss}).pX(1,:,choristersort.(regions{rr}).(celltypes{tt}).(statenames{ss}).(synchtypes{st})))')
        if ss ==3
            xlabel('All MUA (lognorm)')
        end
        if rr == 1
            ylabel({statenames{ss},[(celltypes{tt}), ' Cells, Chor. Sort']})
        end
    end
    end
end
NiceSave('SpikeTriggeredPopDist',figfolder,[])