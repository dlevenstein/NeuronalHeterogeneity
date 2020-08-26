reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/ISIPopTimescaleAnalysis']; 

[baseNames] = getDatasetBasenames();
regions = {'THAL','vCTX','fCTX','BLA','PIR','CA1'};
for rr = 1:6
    %Get baseNames
    if rr == 4 || rr == 5
        ISIPopTimescale_ALL = GetMatResults([figfolder,'_',regions{rr}],['ISIPopTimescaleAnalysis','_',regions{rr}]);
    else
        ISIPopTimescale_ALL = GetMatResults(figfolder,'ISIPopTimescaleAnalysis',...
            'baseNames',baseNames.(regions{rr}));
    end
    ISIPopTimescale_ALL = bz_CollapseStruct(ISIPopTimescale_ALL);
    
    MutInf.(regions{rr}) = bz_CollapseStruct(ISIPopTimescale_ALL.MutInf,'match','justcat',true);
    %MUAConditionalISIDist.(regions{rr}) = bz_CollapseStruct(ISIPopTimescale_ALL.MUAConditionalISIDist,3,'justcat',true);
    %MUAConditionalISIDist_all.(regions{rr}) = bz_CollapseStruct(ISIPop_ALL.MUAConditionalISIDist_all,'match','justcat');

end
%%
synchrate = {'rate','synch'};
spikerate = {'cellrate','cellspike','cellISI'};
statenames = {'WAKEstate','NREMstate','REMstate'};
celltypes = {'pE','pI'};
cellcolor = {'k','r'};
statecolors = {[0 0 0],[0 0 1],[1 0 0]};

cellthresh.pE = 20;
cellthresh.pI = 5;

%%
for rr = 1:6
for cellsr = 1:3
    for ss = 1:3
        
        for popsr = 1:2
            for poptt = 1:2
                for celltt = 1:2
                    usecells = MutInf.(regions{rr}).CellClass.(celltypes{celltt}) & 
meanMI.(regions{rr}).(spikerate{cellsr}).(statenames{ss}).(synchrate{popsr}).(celltypes{poptt}).(celltypes{celltt}) = ...
    nanmedian(MutInf.(regions{rr}).(spikerate{cellsr}).(statenames{ss}).(synchrate{popsr}).(celltypes{poptt})(usecells,:),1);
                end
            end
        end
    end
end
end

%%
figure
% for sr = 1:2
% subplot(2,2,sr)
% plot(MutInf.(spikerate{sr}).(statenames{2}).synch.(celltypes{1}),...
%     MutInf.(spikerate{sr}).(statenames{2}).rate.(celltypes{1}),'.')
% axis tight
% box off
% hold on
% legend
% UnityLine
% xlabel('Pop Synch');ylabel('Pop Rate');title(spikerate{sr})
% end

for popsr = 1:2
for cellsr = 1:3
    for celltt = 1:2
subplot(4,3,(popsr-1)*6+cellsr+(celltt-1)*3)
hold on
%plot(log10(timescales),(MutInf.(spikerate{sr}).(statenames{1}).synch.(celltypes{1})(MutInf.CellClass.(celltypes{tt}),:)),'k')
%plot(log10(timescales),(MutInf.(spikerate{sr}).(statenames{1}).synch.(celltypes{2})(MutInf.CellClass.(celltypes{tt}),:)),'r')

for poptt = 1:2
plot(log10(timescales),(meanMI.(regions{rr}).(spikerate{cellsr}).(statenames{ss}).(synchrate{popsr}).(celltypes{poptt}).(celltypes{celltt})),...
    'color',cellcolor{poptt})
end

LogScale('x',10)
xlabel('Time Window (s)');
if cellsr == 1
    ylabel({celltypes{celltt},'MI'})
end
if celltt==1
title([ spikerate{cellsr}, ' by pop ',(synchrate{popsr})])
end
    end
end
end
NiceSave('CellModulationByPop',figfolder,[])