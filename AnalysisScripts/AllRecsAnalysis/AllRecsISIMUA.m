reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/ISIModesbyPopActivityAnalysis']; 

[baseNames] = getDatasetBasenames();
regions = {'THAL','vCTX','fCTX','BLA','PIR','CA1'};
for rr = 1:6
    %Get baseNames
    if rr == 4 || rr == 5
        ISIPop_ALL = GetMatResults([figfolder,'_',regions{rr}],['ISIModesbyPopActivityAnalysis','_',regions{rr}]);
    else
        ISIPop_ALL = GetMatResults(figfolder,'ISIModesbyPopActivityAnalysis',...
            'baseNames',baseNames.(regions{rr}));
    end
    ISIPop_ALL = bz_CollapseStruct(ISIPop_ALL);
    
    PopCorr.(regions{rr}) = bz_CollapseStruct(ISIPop_ALL.PopCorr,'match','justcat',true);
    MUAConditionalISIDist.(regions{rr}) = bz_CollapseStruct(ISIPop_ALL.MUAConditionalISIDist,3,'justcat',true);
    %MUAConditionalISIDist_all.(regions{rr}) = bz_CollapseStruct(ISIPop_ALL.MUAConditionalISIDist_all,'match','justcat');

end
%%
synchrate = {'rate','synch'};
statenames = {'WAKEstate','NREMstate','REMstate'};
celltypes = {'pE','pI'};
cellcolor = {'k','r'};
statecolors = {[0 0 0],[0 0 1],[1 0 0]};

cellthresh.pE = 20;
cellthresh.pI = 5;
%%
for sr = 1:2
for rr = 2:6
for ss = 1:3
    for tt = 1:length(celltypes)
        
        if ~isfield(MUAConditionalISIDist.(regions{rr}).(statenames{ss}).(synchrate{sr}),(celltypes{tt}))
            continue
        end
        for tt2 = 1:length(celltypes)
            usecells = PopCorr.(regions{rr}).CellClass.(celltypes{tt2}) & ...
                PopCorr.(regions{rr}).cellcount.(celltypes{tt})>cellthresh.(celltypes{tt});
            
            
            MeanCondISI.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}).Dist.pYX = ...
                nanmean(MUAConditionalISIDist.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt}).Dist.pYX(:,:,...
                (usecells)),3);
            
            MeanCondISI.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}).Dist.SpikeRate = ...
                nanmean(MUAConditionalISIDist.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt}).Dist.SpikeRate(:,:,...
                (usecells)),3);
            
            MeanCondISI.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}).Dist.Xbins = ...
                nanmean(MUAConditionalISIDist.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt}).Dist.Xbins(:,:,...
                (usecells)),3);
            MeanCondISI.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}).Dist.Ybins = ...
                nanmean(MUAConditionalISIDist.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt}).Dist.Ybins(:,:,...
                (usecells)),3);
            
        end
        
    end
end
end

end
%%
for sr = 1:2
figure
for rr = 1:6
for tt = 1:2
    for ss = 1:3
        if ~isfield(PopCorr.(regions{rr}).(statenames{ss}).(synchrate{sr}),(celltypes{tt}))
            continue
        end
        
        usecells = PopCorr.(regions{rr}).cellcount.(celltypes{tt})>cellthresh.(celltypes{tt});
        
        subplot(6,6,(tt-1)*18+(ss-1)*6+rr)
        plot(PopCorr.(regions{rr}).(statenames{ss}).GSrate(usecells),...
            PopCorr.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt})(usecells),'k.','markersize',1)
        hold on
        plot(log10(PopCorr.(regions{rr}).(statenames{ss}).meanRate(PopCorr.(regions{rr}).CellClass.pI&usecells)),...
            PopCorr.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt})(PopCorr.(regions{rr}).CellClass.pI&usecells),...
            'r.','markersize',1)
        axis tight
        box off
        plot(xlim(gca),[0 0],'k--')
        xlabel('GS/Mean Rate (Hz)');%ylabel([(celltypes{tt}),' Corr'])
        LogScale('x',10,'exp',true,'nohalf',true)
        if ss == 1
            title(regions{rr})
        end
        if rr ==1
            ylabel([(celltypes{tt}),' ',(synchrate{sr}),' Corr'])
        end
%         subplot(4,4,(tt-1)*4+ss+8)
%         plot(PopCorr.(regions{rr}).(statenames{ss}).GSCV,PopCorr.(regions{rr}).(statenames{ss}).(celltypes{tt}),'.')
%         hold on
%         axis tight
%         box off
%         plot(xlim(gca),[0 0],'k--')
%         xlabel('GS CV');ylabel([(celltypes{tt}),' Corr'])
%         %LogScale('x',10,'exp',false,'nohalf',true)
        
    end 
    
end
end
NiceSave(['MUACorrandGSRate_',(synchrate{sr})],figfolder,[])
end
%%
%% 
for sr = 1:2
for tt2 = 1:length(celltypes)
figure
for rr = 1:6
for ss = 1:3
    for tt = 1:length(celltypes) %Pop
        if ~isfield( MeanCondISI.(regions{rr}).(statenames{ss}).(synchrate{sr}),(celltypes{tt})) | ...
            ~isfield(MeanCondISI.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt}),(celltypes{tt2})) | ...
            all(isnan(MeanCondISI.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}).Dist.pYX))
            continue
        end
     %Ref

    subplot(6,6,(tt-1)*18+(ss-1)*6+rr)
        imagesc(MeanCondISI.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}).Dist.Xbins,...
            MeanCondISI.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}).Dist.Ybins,...
            MeanCondISI.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}).Dist.pYX')
        hold on
        plot(MeanCondISI.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}).Dist.Xbins,...
            -log10(MeanCondISI.(regions{rr}).(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}).Dist.SpikeRate),...
            'r','LineWidth',2)
           
        
        LogScale('y',10,'exp',true,'nohalf',true)
        if rr == 1
            ylabel({(statenames{ss}),[(celltypes{tt2}),' ISI (s)']})
        elseif rr == 6
            set (gca,'yticklabels',[])
            bz_AddRightRateAxis
        else
            set (gca,'yticklabels',[])
        end
        
        if ss == 3
            xlabel([(celltypes{tt}),' ',(synchrate{sr})]);
        else
        set (gca,'xticklabels',[])
        end
        if ss == 1
            title(regions{rr})  
        end
        
    end
    end 
end
NiceSave(['ISIbyMUA_',(synchrate{sr})],figfolder,celltypes{tt2})
end
end



%%
