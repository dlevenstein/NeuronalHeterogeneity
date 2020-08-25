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
statenames = {'WAKEstate','NREMstate','REMstate'};
celltypes = {'pE','pI'};
cellcolor = {'k','r'};
statecolors = {[0 0 0],[0 0 1],[1 0 0]};
%%
for rr = 1:6
for ss = 1:3
    for tt = 1:length(celltypes)
        if ~isfield(MUAConditionalISIDist.(regions{rr}).(statenames{ss}),(celltypes{tt}))
            continue
        end
        for tt2 = 1:length(celltypes)
            MeanCondISI.(regions{rr}).(statenames{ss}).(celltypes{tt}).(celltypes{tt2}).Dist.pYX = ...
                nanmean(MUAConditionalISIDist.(regions{rr}).(statenames{ss}).(celltypes{tt}).Dist.pYX(:,:,...
                (PopCorr.(regions{rr}).CellClass.(celltypes{tt2}))),3);
            
            MeanCondISI.(regions{rr}).(statenames{ss}).(celltypes{tt}).(celltypes{tt2}).Dist.SpikeRate = ...
                nanmean(MUAConditionalISIDist.(regions{rr}).(statenames{ss}).(celltypes{tt}).Dist.SpikeRate(:,:,...
                (PopCorr.(regions{rr}).CellClass.(celltypes{tt2}))),3);
            
            MeanCondISI.(regions{rr}).(statenames{ss}).(celltypes{tt}).(celltypes{tt2}).Dist.Xbins = ...
                nanmean(MUAConditionalISIDist.(regions{rr}).(statenames{ss}).(celltypes{tt}).Dist.Xbins(:,:,...
                (PopCorr.(regions{rr}).CellClass.(celltypes{tt2}))),3);
            MeanCondISI.(regions{rr}).(statenames{ss}).(celltypes{tt}).(celltypes{tt2}).Dist.Ybins = ...
                nanmean(MUAConditionalISIDist.(regions{rr}).(statenames{ss}).(celltypes{tt}).Dist.Ybins(:,:,...
                (PopCorr.(regions{rr}).CellClass.(celltypes{tt2}))),3);
            
        end
        
    end
end
end


%%

figure
for rr = 1:6
for tt = 1:2
    for ss = 1:3
        if ~isfield(PopCorr.(regions{rr}).(statenames{ss}),(celltypes{tt}))
            continue
        end
        subplot(6,6,(tt-1)*18+(ss-1)*6+rr)
        plot(PopCorr.(regions{rr}).(statenames{ss}).GSrate,PopCorr.(regions{rr}).(statenames{ss}).(celltypes{tt}),'.')
        hold on
        axis tight
        box off
        plot(xlim(gca),[0 0],'k--')
        xlabel('GS Rate (Hz)');ylabel([(celltypes{tt}),' Corr'])
        LogScale('x',10,'exp',true,'nohalf',true)
        if ss == 1
            title(regions{rr})
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
NiceSave('MUACorrandGSRate',figfolder,[])

%%
%% 
for tt2 = 1:length(celltypes)
figure
for rr = 1:6
for ss = 1:3
    for tt = 1:length(celltypes) %Pop
        if ~isfield( MeanCondISI.(regions{rr}).(statenames{ss}),(celltypes{tt})) | ...
            ~isfield(MeanCondISI.(regions{rr}).(statenames{ss}).(celltypes{tt}),(celltypes{tt2})) | ...
            all(isnan(MeanCondISI.(regions{rr}).(statenames{ss}).(celltypes{tt}).(celltypes{tt2}).Dist.pYX))
            continue
        end
     %Ref

    subplot(6,6,(tt-1)*18+(ss-1)*6+rr)
        imagesc(MeanCondISI.(regions{rr}).(statenames{ss}).(celltypes{tt}).(celltypes{tt2}).Dist.Xbins,...
            MeanCondISI.(regions{rr}).(statenames{ss}).(celltypes{tt}).(celltypes{tt2}).Dist.Ybins,...
            MeanCondISI.(regions{rr}).(statenames{ss}).(celltypes{tt}).(celltypes{tt2}).Dist.pYX')
        hold on
        plot(MeanCondISI.(regions{rr}).(statenames{ss}).(celltypes{tt}).(celltypes{tt2}).Dist.Xbins,...
            -log10(MeanCondISI.(regions{rr}).(statenames{ss}).(celltypes{tt}).(celltypes{tt2}).Dist.SpikeRate),...
            'r','LineWidth',2)
           
        xlabel([(celltypes{tt}),' Rate']);ylabel([(celltypes{tt2}),' ISI (s)'])
        LogScale('y',10,'exp',true,'nohalf',true)
        bz_AddRightRateAxis
        if tt ==1 & tt2 == 1
            title(statenames{ss})
        end
    end
    end 
end
NiceSave('ISIbyMUA',figfolder,celltypes{tt2})
end




%%
