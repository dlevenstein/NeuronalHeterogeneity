reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/SpikeStatsbyBrainState'];

datasetPath.fCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onDesktop/BW_CTX';
datasetPath.CA1 = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onDesktop/AG_HPC';
%regions = {'fCTX','CA1'};
%datasetPath = figfolder;

[BrainStateandSpikingAll,baseNames] = GetMatResults(figfolder,'SpikeStatsbyBrainState','select',true);
BrainStateandSpikingAll = bz_CollapseStruct(BrainStateandSpikingAll);
thisregion = 'CA1';

%%
ISIbyPSS = bz_CollapseStruct(BrainStateandSpikingAll.ISIbyPSS,'match','justcat',true);
ISIbytheta = bz_CollapseStruct(BrainStateandSpikingAll.ISIbytheta,'match','justcat',true);
BShist = bz_CollapseStruct(BrainStateandSpikingAll.BShist,1,'mean',true);
BShist.std = bz_CollapseStruct(BrainStateandSpikingAll.BShist,1,'std',true);
BShist_all = bz_CollapseStruct(BrainStateandSpikingAll.BShist,1,'justcat',true);

celltypes = fieldnames(ISIbyPSS.celltypeidx);

for tt = 1:length(celltypes)
    ISIbyPSS.pop.(celltypes{tt}) = nanmean(ISIbyPSS.pYX(:,:,ISIbyPSS.celltypeidx.(celltypes{tt})),3);
    ISIbytheta.pop.(celltypes{tt}) = nanmean(ISIbytheta.pYX(:,:,ISIbytheta.celltypeidx.(celltypes{tt})),3);
end

%%
figure
subplot(2,2,1)
imagesc(BShist_all.bins(1,:),[0 1],BShist_all.ALL.PSS)
subplot(2,2,2)
imagesc(BShist_all.bins(1,:),[0 1],BShist_all.ALL.thratio)
%%
states = {'WAKEstate','NREMstate','REMstate'};
statecolors = {'k','b','r'};
%%
figure
   
for tt = 1:length(celltypes)
subplot(4,3,tt*3-2)
    imagesc(ISIbyPSS.Xbins(1,:,1),ISIbyPSS.Ybins(1,:,1), ISIbyPSS.pop.(celltypes{tt})')
    %hold on
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    %LogScale('y',10)
    ylabel({(celltypes{tt}),'ISI (log(s))'});
    %title((celltypes{tt}))
   % colorbar
    if tt ==1 
        caxis([0 0.023])
    elseif tt==2
         caxis([0 0.032])
    end
    xlim(ISIbyPSS.Xbins(1,[1 end],1))
    if strcmp(thisregion,'CA1')
        xlim([ISIbyPSS.Xbins(1,1,1) 1.9])
    end
end 


for tt = 1:length(celltypes)
subplot(4,3,tt*3-1)
    imagesc(ISIbytheta.Xbins(1,:,1),ISIbytheta.Ybins(1,:,1), ISIbytheta.pop.(celltypes{tt})')
    %hold on
    %plot(CONDXY.Xbins(1,:,1),meanthetabyPOP.(celltypes{tt}),'w')
    axis xy
    %LogScale('y',10)
    %ylabel({(celltypes{tt}),'ISI (s)'});xlabel('Theta Ratio')
    %title((celltypes{tt}))
   % colorbar
    if tt ==1 
        caxis([0 0.025])
    elseif tt==2
         caxis([0 0.035])
    end
    xlim(ISIbytheta.Xbins(1,[1 end],1))
end 

subplot(6,3,10)
    hold on
    for ss = 1:3
        errorshade(BShist.bins,BShist.(states{ss}).PSS,...
            BShist.std.(states{ss}).PSS,BShist.std.(states{ss}).PSS,statecolors{ss},'scalar')
    plot(BShist.bins,BShist.(states{ss}).PSS,'color',statecolors{ss})
    end
    xlabel('PSS')
    
    box off
    axis tight
    xlim(ISIbyPSS.Xbins(1,[1 end],1))
        if strcmp(thisregion,'CA1')
        xlim([ISIbyPSS.Xbins(1,1,1) 1.9])
    end
    
subplot(6,3,11)
    hold on
    for ss = [1 3]
        errorshade(BShist.bins,BShist.(states{ss}).thratio,...
            BShist.std.(states{ss}).thratio,BShist.std.(states{ss}).thratio,statecolors{ss},'scalar')
    plot(BShist.bins,BShist.(states{ss}).thratio,'color',statecolors{ss})
    end
    xlabel('Theta Ratio')
        box off
    axis tight
    xlim(ISIbytheta.Xbins(1,[1 end],1))


 NiceSave('ISIbyStateVars',figfolder,thisregion)
