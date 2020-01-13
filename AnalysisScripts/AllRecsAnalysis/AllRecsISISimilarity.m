reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/ISIDistSimilarityAnalysis'];

datasetPath.fCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/BW_CTX';
datasetPath.CA1 = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AG_HPC';
datasetPath.vCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/YS_CTX';
datasetPath.THAL = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AP_THAL';
regions = {'THAL','vCTX','fCTX','CA1'};

celltypes = {'pE','pI'};

for rr = 1:length(regions)
    disp(['Loading ',regions{rr}])
    %[ISIStats.(regions{rr}),baseNames] = bz_LoadCellinfo(datasetPath.(regions{rr}),'ISIStats','dataset',true,'catall',true);
    %CellClass.(regions{rr}) = bz_LoadCellinfo(datasetPath.(regions{rr}),'CellClass','dataset',true,'catall',true,'baseNames',baseNames);

    [ISISimilarityAll,baseNames] = bz_LoadAnalysisResults(datasetPath.(regions{rr}),'ISIDistSimilarityAnalysis','dataset',true);
    
    
    %PopActivityAll = GetMatResults(figfolder,'SpikeStatsbyPopActivityAnalysis','baseNames',baseNames);
    ISISimilarityAll = bz_CollapseStruct(ISISimilarityAll);
    
    allpairs.(regions{rr}) = bz_CollapseStruct(ISISimilarityAll.allpairs,3,'justcat',true);
    allISIs.(regions{rr}) =bz_CollapseStruct(ISISimilarityAll.allISIs,2,'justcat',true);
    lowestpairISI.(regions{rr}) =bz_CollapseStruct(ISISimilarityAll.lowestpairISI,3,'justcat',true);
    
    clear ISISimilarityAll
end

%%
ISIthreshold = 800;

for rr = 1:length(regions)
    for ii = 1:3
        for jj = 1:3
            for cc = 1:2
                newpairs.(regions{rr}).(celltypes{cc}){ii,jj} = ...
                    cat(1,allpairs.(regions{rr}).(celltypes{cc}){ii,jj,:});

                numISIs.(regions{rr}).(celltypes{cc}){ii,jj} = ...
                    cat(1,lowestpairISI.(regions{rr}).(celltypes{cc}){ii,jj,:});
                abovethresh = numISIs.(regions{rr}).(celltypes{cc}){ii,jj}>ISIthreshold;
                
                simmatrices.(regions{rr}).(celltypes{cc})(ii,jj) = ...
                    nanmedian(newpairs.(regions{rr}).(celltypes{cc}){ii,jj}(abovethresh));
            
                if rr ==1
                    newpairs.ALL.(celltypes{cc}){ii,jj} = ...
                        newpairs.(regions{rr}).(celltypes{cc}){ii,jj}(abovethresh);
                else
                    newpairs.ALL.(celltypes{cc}){ii,jj} = ...
                        [newpairs.ALL.(celltypes{cc}){ii,jj};newpairs.(regions{rr}).(celltypes{cc}){ii,jj}(abovethresh)];

                    simmatrices.ALL.(celltypes{cc})(ii,jj) = ...
                        nanmedian(newpairs.ALL.(celltypes{cc}){ii,jj});
                end
            end
        end
    end
end
%
for rr = 1:length(regions)
    for ii = 1:6
        for jj = 1:6
            newpairs.(regions{rr}).difft{ii,jj} = cat(1,allpairs.(regions{rr}).difft{ii,jj,:});
            
            numISIs.(regions{rr}).difft{ii,jj} = ...
                cat(1,lowestpairISI.(regions{rr}).difft{ii,jj,:});
            abovethresh = numISIs.(regions{rr}).difft{ii,jj}>ISIthreshold;
            
            simmatrices.(regions{rr}).difft(ii,jj) = nanmedian(newpairs.(regions{rr}).difft{ii,jj}(abovethresh));
            
            if rr ==1
                newpairs.ALL.difft{ii,jj} = newpairs.(regions{rr}).difft{ii,jj}(abovethresh);
            else
                newpairs.ALL.difft{ii,jj} = [newpairs.ALL.difft{ii,jj};newpairs.(regions{rr}).difft{ii,jj}(abovethresh)];
                
                simmatrices.ALL.difft(ii,jj) = nanmedian(newpairs.ALL.difft{ii,jj});
            end
        end
    end
end


%%
figure
for rr = 1:4

for cc = 1:length(celltypes)
subplot(4,4,(cc-1).*4+rr)
imagesc(simmatrices.(regions{rr}).(celltypes{cc}))
alpha(gca,single(~isnan(simmatrices.(regions{rr}).(celltypes{cc}))))
caxis([0.05 0.4])
box off
set(gca,'ytick',[1:3]);set(gca,'xtick',[1:3]);
set(gca,'yticklabel',{'N','W','R'})
set(gca,'xticklabel',{'N','W','R'})
ylabel(celltypes{cc});
crameri tokyo

    title((regions{rr}))

    
subplot(2,2,2+cc)
imagesc(simmatrices.ALL.(celltypes{cc}))
alpha(gca,single(~isnan(simmatrices.ALL.(celltypes{cc}))))
caxis([0.05 0.4])
colorbar
box off
set(gca,'ytick',[1:3]);set(gca,'xtick',[1:3]);
set(gca,'yticklabel',{'N','W','R'})
set(gca,'xticklabel',{'N','W','R'})
crameri tokyo

    title(celltypes{cc})
end
end

NiceSave('ISISimilarity_WithinCell',figfolder,'')
%%
figure
for rr = 1:4
subplot(2,2,rr)
imagesc(simmatrices.(regions{rr}).difft)
alpha(gca,single(~isnan(simmatrices.(regions{rr}).difft)))
colorbar
caxis([0.25 2.5])
title('Diff. Cells')
box off
crameri tokyo
set(gca,'ytick',[1:6]);set(gca,'xtick',[1:6]);
set(gca,'yticklabel',{'N','W','R','N','W','R'})
set(gca,'xticklabel',{'N','W','R','N','W','R'})
title((regions{rr}))
end
NiceSave('ISISimilarity_BetweenCell',figfolder,'')
%%
figure

imagesc(simmatrices.ALL.difft)
alpha(gca,single(~isnan(simmatrices.ALL.difft)))
colorbar
caxis([0.25 2.5])
title('Diff. Cells')
box off
crameri tokyo
set(gca,'ytick',[1:6]);set(gca,'xtick',[1:6]);
set(gca,'yticklabel',{'N','W','R','N','W','R'})
set(gca,'xticklabel',{'N','W','R','N','W','R'})
title('All Pairs')

NiceSave('ISISimilarity_AllCells',figfolder,'')


%% Build the All ISI matrix
%Subsample ISIs
keepISIs.numPerregion = 1800;
keepISIs.numISIthresh = 800;
for rr = 1:length(regions)
    OKISIS = (allISIs.(regions{rr}).celltype.pE | allISIs.(regions{rr}).celltype.pI) & ...
        ~any(isnan(allISIs.(regions{rr}).hists),1) &  allISIs.(regions{rr}).numISIs>keepISIs.numISIthresh;
   keepISIs.(regions{rr}) = randsample(find(OKISIS), keepISIs.numPerregion);
end
% allISIs = [ISIs.THAL(keepISIs.THAL),ISIs.vCTX(keepISIs.vCTX),...
%     ISIs.fCTX(keepISIs.fCTX),ISIs.CA1(keepISIs.CA1)];
%numISIs = cellfun(@length,allISIs);
allISIhists.hists = [allISIs.THAL.hists(:,keepISIs.THAL),allISIs.vCTX.hists(:,keepISIs.vCTX),...
     allISIs.fCTX.hists(:,keepISIs.fCTX),allISIs.CA1.hists(:,keepISIs.CA1)];
 
 %%
ALLregions = [1.*ones(1,keepISIs.numPerregion),2.*ones(1,keepISIs.numPerregion),...
    3.*ones(1,keepISIs.numPerregion),4.*ones(1,keepISIs.numPerregion)];
ALLcelltypes.pI = [allISIs.THAL.celltype.pI(keepISIs.THAL),allISIs.vCTX.celltype.pI(keepISIs.vCTX),...
    allISIs.fCTX.celltype.pI(keepISIs.fCTX),allISIs.CA1.celltype.pI(keepISIs.CA1)];
ALLcelltypes.pE = [allISIs.THAL.celltype.pE(keepISIs.THAL),allISIs.vCTX.celltype.pE(keepISIs.vCTX),...
    allISIs.fCTX.celltype.pE(keepISIs.fCTX),allISIs.CA1.celltype.pE(keepISIs.CA1)];
states.ALL = [allISIs.THAL.state(keepISIs.THAL),allISIs.vCTX.state(keepISIs.vCTX),...
    allISIs.fCTX.state(keepISIs.fCTX),allISIs.CA1.state(keepISIs.CA1)];

[jregions,iregions] = meshgrid(ALLregions,ALLregions);
[jiscelltype.pE,iiscelltype.pE] = meshgrid(ALLcelltypes.pE,ALLcelltypes.pE);
[jiscelltype.pI,iiscelltype.pI] = meshgrid(ALLcelltypes.pI,ALLcelltypes.pI);
[jstates,istates] = meshgrid(states.ALL,states.ALL);
%%
numcells = length(ALLregions);
KLDIST = nan(numcells);
%parpool
for ii = 1:numcells
    bz_Counter(ii,numcells,'Cell');
            KLDIST(:,ii) = KLDiv(allISIhists.hists',allISIhists.hists(:,ii)','symmetric',true);

end

%% For KS Statistic

% numcells = length(allISIs);
% KSSTAT = nan(numcells);
% %parpool
% for ii = 1:numcells
%     bz_Counter(ii,numcells,'Cell');
%     parfor jj = ii:numcells
%         if ii==jj %For Same Cell, calculate first/last half spikes
%             [~,~,KSSTAT(ii,jj)] = kstest2(allISIs{ii}(1:round(end/2)),allISIs{jj}(round(end/2):end));
%         else
%             [~,~,KSSTAT(ii,jj)] = kstest2(allISIs{ii},allISIs{jj});
%         end
%     end
%     for jj = ii:numcells
%         KSSTAT(jj,ii) = KSSTAT(ii,jj);
%     end
% end

%%
statenames = {'WAKEstate','NREMstate','REMstate'};
statecolors = {[0 0 0],[0 0 1],[1 0 0]};


%%
%simmatrices.(statenames{ss}).(celltypes{cc})

for ss = 1:3
    for cc = 1:2
        simmatrices.(statenames{ss}).(celltypes{cc}) = nan(4);
        simmatrices.ALLregions.(celltypes{cc}) = nan(4);
    end
end
simmatrices.ALLregions.(celltypes{cc}) = nan(4);
for rr1 = 1:length(regions)
    for rr2 = rr1:-1:1
        for cc = 1:2
            for ss = 1:3
            
                allpairs.(statenames{ss}).(celltypes{cc}){rr1,rr2} = ...
                    KLDIST(istates==ss & jstates==ss & ...
                    iiscelltype.(celltypes{cc}) & jiscelltype.(celltypes{cc}) & ...
                    iregions == rr1 & jregions == rr2);

                simmatrices.(statenames{ss}).(celltypes{cc})(rr1,rr2) = ...
                    nanmedian(allpairs.(statenames{ss}).(celltypes{cc}){rr1,rr2});
            end
                allpairs.ALLregions.(celltypes{cc}){rr1,rr2} = ...
                    KLDIST(iiscelltype.(celltypes{cc}) & jiscelltype.(celltypes{cc}) & ...
                    iregions == rr1 & jregions == rr2);

                simmatrices.ALLregions.(celltypes{cc})(rr1,rr2) = ...
                    nanmedian(allpairs.(statenames{ss}).(celltypes{cc}){rr1,rr2});
        end
    end
end

%%
figure
for cc = 1:2
    for ss = 1:3

        subplot(3,3,(cc-1)*3+ss)
            imagesc(simmatrices.(statenames{ss}).(celltypes{cc}))
            %colorbar
            alpha(gca,single(~isnan(simmatrices.(statenames{ss}).(celltypes{cc}))))
            colorbar
            caxis([0.25 1.75])
            if cc == 1
            title(statenames{ss})
            end
            box off
            crameri tokyo
            set(gca,'ytick',[1:4]);set(gca,'xtick',[1:4])
            set(gca,'yticklabels',regions);set(gca,'xticklabels',regions)
    end
        subplot(3,3,6+cc)
            imagesc(simmatrices.ALLregions.(celltypes{cc}))
            %colorbar
            alpha(gca,single(~isnan(simmatrices.ALLregions.(celltypes{cc}))))
            colorbar
            caxis([0.25 1.75])
            if cc == 1
            title('ALL')
            end
            box off
            crameri tokyo
            set(gca,'ytick',[1:4]);set(gca,'xtick',[1:4])
            set(gca,'yticklabels',regions);set(gca,'xticklabels',regions)
end
NiceSave('ISISimilarity_BetweenRegions',figfolder,'')

%% tSNE Map
valid = KLDIST;
for ii =1:numcells
        valid(ii,ii)=0;
end
%%
% %MDS Map
% Y = cmdscale(valid);
%Good: perp=20, valid not squared
%Try different perp parameters
clear Y
clear Ysq
perps = 10:10:40;
for pp = 1:length(perps)
close all

P = d2p(valid, perps(pp), 1e-6); 
Y(:,:,pp) = tsne_p(P, ALLcelltypes.pE, 2, 2000);

close all
P = d2p(valid.^2, perps(pp), 1e-6); 
Ysq(:,:,pp) = tsne_p(P, ALLcelltypes.pE, 2, 2000);
end
%% 2D
figure
for pp = 1:length(perps)
%legend('NREM','WAKE','REM')
subplot(4,4,pp)
hold on
for rr = 1:length(regions)
    plot(Y(ALLregions==rr,1,pp),Y(ALLregions==rr,2,pp),'.','markersize',3)
end
%legend(regions,'location','northoutside')
axis tight
box on
xticks(gca,[]);yticks(gca,[])
title(perps(pp))
if pp==1
    ylabel('KS^1')
end
subplot(4,4,pp+4)
hold on
axis tight
box on
xticks(gca,[]);yticks(gca,[])
plot(Y(ALLcelltypes.pE==1,1,pp),Y(ALLcelltypes.pE==1,2,pp),'.k','markersize',3)
plot(Y(ALLcelltypes.pI==1,1,pp),Y(ALLcelltypes.pI==1,2,pp),'.r','markersize',3)
%legend(celltypes,'location','northoutside')


subplot(4,4,pp+8)
hold on
for rr = 1:length(regions)
    plot(Ysq(ALLregions==rr,1,pp),Ysq(ALLregions==rr,2,pp),'.','markersize',3)
end
%legend(regions,'location','northoutside')
axis tight
box on
if pp==1
    ylabel('KS^2')
end
xticks(gca,[]);yticks(gca,[])
subplot(4,4,pp+12)
hold on
axis tight
box on

xticks(gca,[]);yticks(gca,[])
plot(Ysq(ALLcelltypes.pE==1,1,pp),Ysq(ALLcelltypes.pE==1,2,pp),'.k','markersize',3)
plot(Ysq(ALLcelltypes.pI==1,1,pp),Ysq(ALLcelltypes.pI==1,2,pp),'.r','markersize',3)
end
NiceSave('tSNE_ComparePerplexity',figfolder,'')
%%
close all
clear tSNEmap
perplexity = 20;
P = d2p(valid.^2, perplexity, 1e-6); 
tSNEmap(:,:) = tsne_p(P, ALLcelltypes.pE, 2, 2000);
%%

figure
subplot(3,3,[5 8])
hold on
plot(tSNEmap(states.ALL==2,1),tSNEmap(states.ALL==2,2),'k.','markersize',3)
plot(tSNEmap(states.ALL==1,1),tSNEmap(states.ALL==1,2),'b.','markersize',3)
plot(tSNEmap(states.ALL==3,1),tSNEmap(states.ALL==3,2),'r.','markersize',3)
axis tight
box on
xticks(gca,[]);yticks(gca,[])
legend(statenames,'location','northoutside')

%legend('NREM','WAKE','REM')
subplot(3,3,[6 9])
hold on
for rr = 1:length(regions)
    plot(tSNEmap(ALLregions==rr,1),tSNEmap(ALLregions==rr,2),'.','markersize',3)
end
legend(regions,'location','northoutside')
axis tight
box on
xticks(gca,[]);yticks(gca,[])
subplot(3,3,[4 7])
hold on
axis tight
box on
xticks(gca,[]);yticks(gca,[])
plot(tSNEmap(ALLcelltypes.pE==1,1),tSNEmap(ALLcelltypes.pE==1,2),'.k','markersize',3)
plot(tSNEmap(ALLcelltypes.pI==1,1),tSNEmap(ALLcelltypes.pI==1,2),'.r','markersize',3)
legend(celltypes,'location','northoutside')
NiceSave('tSNE_Map',figfolder,'')


%% Subsets...

% Each State Only
for ss = 1:3
    close all
    P = d2p(valid(states.ALL==ss,states.ALL ==ss).^2, perplexity, 1e-6); 
    statetSNEmap.(statenames{ss}) = tsne_p(P, ALLcelltypes.pE(states.ALL==ss), 2, 1000);
end

%%
figure

for ss = 1:3
subplot(2,3,ss)
hold on
for rr = 1:length(regions)
    plot(statetSNEmap.(statenames{ss})(ALLregions(states.ALL==ss)==rr,1),...
        statetSNEmap.(statenames{ss})(ALLregions(states.ALL==ss)==rr,2),'.','markersize',3)
end
legend(regions,'location','northoutside')
axis tight
box on
xticks(gca,[]);yticks(gca,[])
subplot(2,3,ss+3)
hold on
axis tight

xticks(gca,[]);yticks(gca,[])
plot(statetSNEmap.(statenames{ss})(ALLcelltypes.pE(states.ALL==ss)==1,1),...
    statetSNEmap.(statenames{ss})(ALLcelltypes.pE(states.ALL==ss)==1,2),'.k','markersize',3)
plot(statetSNEmap.(statenames{ss})(ALLcelltypes.pI(states.ALL==ss)==1,1),...
    statetSNEmap.(statenames{ss})(ALLcelltypes.pI(states.ALL==ss)==1,2),'.r','markersize',3)
legend(celltypes,'location','northoutside')
box on
end
NiceSave('tSNE_Map_SeparateStates',figfolder,'')
%%

% E only
    close all
    P = d2p(valid(ALLcelltypes.pE==1,ALLcelltypes.pE==1).^2, perplexity, 1e-6); 
    statetSNEmap.pE.all = tsne_p(P, states.ALL(ALLcelltypes.pE==1), 2, 1000);

for ss = 1:3
    close all
    P = d2p(valid(ALLcelltypes.pE==1 & states.ALL==ss,ALLcelltypes.pE==1 & states.ALL ==ss).^2, perplexity, 1e-6); 
    statetSNEmap.pE.(statenames{ss}) = tsne_p(P, ALLregions(ALLcelltypes.pE==1 & states.ALL==ss), 2, 1000);
end
%%
figure
for ss = 1:3
    subplot(3,3,ss)
    hold on
    for rr = 1:length(regions)
        plot(statetSNEmap.pE.(statenames{ss})(ALLregions(ALLcelltypes.pE==1 & states.ALL==ss)==rr,1),...
            statetSNEmap.pE.(statenames{ss})(ALLregions(ALLcelltypes.pE==1 & states.ALL==ss)==rr,2),'.','markersize',3)
    end
    axis tight
    box on
    xticks(gca,[]);yticks(gca,[])
    title(statenames{ss})
end

subplot(3,3,[5 8])
hold on
plot(statetSNEmap.pE.all(states.ALL(ALLcelltypes.pE==1)==2,1),statetSNEmap.pE.all(states.ALL(ALLcelltypes.pE==1)==2,2),'k.','markersize',3)
plot(statetSNEmap.pE.all(states.ALL(ALLcelltypes.pE==1)==1,1),statetSNEmap.pE.all(states.ALL(ALLcelltypes.pE==1)==1,2),'b.','markersize',3)
plot(statetSNEmap.pE.all(states.ALL(ALLcelltypes.pE==1)==3,1),statetSNEmap.pE.all(states.ALL(ALLcelltypes.pE==1)==3,2),'r.','markersize',3)
axis tight
box on
xticks(gca,[]);yticks(gca,[])
legend(statenames,'location','northoutside')

%legend('NREM','WAKE','REM')
subplot(3,3,[6 9])
hold on
for rr = 1:length(regions)
    plot(statetSNEmap.pE.all(ALLregions(ALLcelltypes.pE==1)==rr,1),statetSNEmap.pE.all(ALLregions(ALLcelltypes.pE==1)==rr,2),'.','markersize',3)
end
legend(regions,'location','northoutside')
axis tight
box on
xticks(gca,[]);yticks(gca,[])

NiceSave('tSNE_Map_pEOnly',figfolder,'')

