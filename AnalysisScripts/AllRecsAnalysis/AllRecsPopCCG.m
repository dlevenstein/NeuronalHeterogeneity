%reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/PopCCGAnalysis'];

% datasetPath.fCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/BW_CTX';
% datasetPath.CA1 = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AG_HPC';
% datasetPath.vCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/YS_CTX';
% datasetPath.THAL = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AP_THAL';
regions = {'CA1'};

%%

for rr = 1:length(regions)
    
    [PopCCGAll,baseNames] = GetMatResults(figfolder,'PopCCGAnalysis');
    PopCCGAll = bz_CollapseStruct(PopCCGAll);
    %thisregion = 'fCTX';
    
    ISICCG.(regions{rr}) = bz_CollapseStruct(PopCCGAll.ISICCG,'match','justcat',true);
    popCCG.(regions{rr}) = bz_CollapseStruct(PopCCGAll.popCCG,3,'mean',true);
    popCCG_all.(regions{rr}) = bz_CollapseStruct(PopCCGAll.popCCG,'match','justcat',true);

    CellClass.(regions{rr}) = bz_CollapseStruct(PopCCGAll.CellClass,'match','justcat',true);
   % ISIStats.(regions{rr}) = bz_CollapseStruct(PopCCGAll.ISIStats,'match','justcat',true);
   clear PopCCGAll
end

%%
states = {'WAKEstate','NREMstate','REMstate'};
celltypes = {'pE','pI'};
for ss = 1:3
for tt = 1:2
    for tt2 = 1:2
        ISICCG.(regions{rr}).(states{ss}).popmean.(celltypes{tt}).(celltypes{tt2}) = ...
            nanmean(ISICCG.(regions{rr}).(states{ss}).(celltypes{tt})(:,:,CellClass.(regions{rr}).(celltypes{tt2})),3);
    end
end
end

%%
numISIbins = 60;
logISIbounds = [0.001 200];
ISICCG.(regions{rr}).logISIbins = linspace(log10(logISIbounds(1)),log10(logISIbounds(2)),numISIbins);
%%

figure
for ss = 1:3
for tt = 1:2
    for tt2 = 1:2
        subplot(4,3,(tt-1)*3+(tt2-1)*6+ss)
            imagesc(ISICCG.(regions{rr}).t_ccg,...
                ISICCG.(regions{rr}).logISIbins,ISICCG.(regions{rr}).(states{ss}).popmean.(celltypes{tt}).(celltypes{tt2})')
            hold on
            %plot(ISIStats.ISIhist.(state).log(cc,:),ISIStats.ISIhist.logbins,'k')
            LogScale('y',10)
            
            if tt == 1
                ColorbarWithAxis([0 2.5],'E Pop Rate')
            else 
                ColorbarWithAxis([0 30],'I Pop Rate')
            end
            if tt==2
            xlabel(['t lag (s) - ',(celltypes{tt2})]);
            end
            if tt ==1 & tt2==1
                title(states{ss})
            end
            if ss == 1
                ylabel('ISI (s)')
            end
            
    end
    
end
end
NiceSave(['CCGbyISI_',(regions{rr})],figfolder,[])

%%
excell = 50;
figure
for ss = 1:3
for tt = 1:2

        subplot(4,3,(tt-1)*3+ss)
            imagesc(ISICCG.(regions{rr}).t_ccg,...
                ISICCG.(regions{rr}).logISIbins,ISICCG.(regions{rr}).(states{ss}).(celltypes{tt})(:,:,excell)')
            hold on
            %plot(ISIStats.ISIhist.(state).log(cc,:),ISIStats.ISIhist.logbins,'k')
            LogScale('y',10)
            
            if tt == 1
                ColorbarWithAxis([0 2.5],'E Pop Rate')
            else 
                ColorbarWithAxis([0 30],'I Pop Rate')
            end
            
            xlabel('t lag (s)');
            %end
            if tt ==1
                title(states{ss})
            end
            if ss == 1
                ylabel('ISI (s)')
            end
            
    end
    
end
NiceSave(['CCGbyISI_example',(regions{rr})],figfolder,[])
%%

figure
for ss = 1:3
    for tt = 1:2
        subplot(3,3,ss+(tt-1)*3)
        plot(popCCG.CA1.(states{ss}).t_ccg,popCCG.CA1.(states{ss}).pop.(celltypes{tt}),'linewidth',2)

        if tt == 1
            title(states{ss})
        end
        
        if tt == 1
            xlabel('t lag (s)')
        elseif tt == 2
            xlabel('t lag (s)')
        end
        if ss == 1
            if tt == 1
                ylabel('pE Rate')
            elseif tt == 2
                ylabel('pI Rate')
            end
        end


    end
end
NiceSave(['PopCCG',(regions{rr})],figfolder,[])


%%
figure
imagesc(popCCG_all.CA1.NREMstate.cells.pE)
colorbar
caxis([0 2.5])