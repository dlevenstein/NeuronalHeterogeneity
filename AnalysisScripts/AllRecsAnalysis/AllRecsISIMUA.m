reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/ISIModesbyPopActivityAnalysis']; 

[baseNames] = getDatasetBasenames();
regions = {'THAL','vCTX','fCTX','CA1','BLA','PIR'};
for rr = 2:4
    %Get baseNames
    if rr == 5 || rr == 6
        ISIPop_ALL = GetMatResults([figfolder,'_',regions{rr}],['ISIModesbyPopActivityAnalysis','_',regions{rr}]);
    else
        ISIPop_ALL = GetMatResults(figfolder,'ISIModesbyPopActivityAnalysis',...
            'baseNames',baseNames.(regions{rr}));
    end
    ISIPop_ALL = bz_CollapseStruct(ISIPop_ALL);
    
    PopCorr.(regions{rr}) = bz_CollapseStruct(ISIPop_ALL.PopCorr,'match','justcat',true);
    MUAConditionalISIDist.(regions{rr}) = bz_CollapseStruct(ISIPop_ALL.MUAConditionalISIDist,3,'justcat',true);
    MUAConditionalISIDist_all.(regions{rr}) = bz_CollapseStruct(ISIPop_ALL.MUAConditionalISIDist_all,'match','justcat');

end

%%
statenames = {'WAKEstate','NREMstate','REMstate'};
celltypes = {'pE','pI'};
cellcolor = {'k','r'};
statecolors = {[0 0 0],[0 0 1],[1 0 0]};

%%
for rr = 2:4
figure
for tt = 1:2
    for ss = 1:3
        subplot(4,4,(tt-1)*4+ss)
        plot(PopCorr.(regions{rr}).(statenames{ss}).GSrate,PopCorr.(regions{rr}).(statenames{ss}).(celltypes{tt}),'.')
        hold on
        axis tight
        box off
        plot(xlim(gca),[0 0],'k--')
        xlabel('GS Rate (Hz)');ylabel([(celltypes{tt}),' Corr'])
        LogScale('x',10,'exp',true,'nohalf',true)
        if tt == 1
            title(statenames{ss})
        end
        
        subplot(4,3,(tt-1)*3+ss+6)
        plot(PopCorr.(regions{rr}).(statenames{ss}).GSCV,PopCorr.(regions{rr}).(statenames{ss}).(celltypes{tt}),'.')
        hold on
        axis tight
        box off
        plot(xlim(gca),[0 0],'k--')
        xlabel('GS CV');ylabel([(celltypes{tt}),' Corr'])
        %LogScale('x',10,'exp',false,'nohalf',true)
        
    end 
    
end
end
%NiceSave('MUACorrandGSRate',figfolder,[])
