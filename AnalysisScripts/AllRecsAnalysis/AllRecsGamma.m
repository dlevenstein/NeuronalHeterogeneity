reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/GammaModeFitAnalysis'];

datasetPath.fCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/BW_CTX';
datasetPath.CA1 = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AG_HPC';
datasetPath.vCTX = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/YS_CTX';
datasetPath.THAL = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/AP_THAL';
datasetPath.BLA = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/GG_BLA';
datasetPath.PIR = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/Datasets/onProbox/GG_BLA';
regions = {'THAL','vCTX','fCTX','BLA','PIR','CA1'};
rnames =  {''    ,''    ,''    ,'bla','pir',''   };
regioncolors = crameri('batlow',length(regions));
celltypes = {'pE','pI'};
cellcolor = {'k','r'};
statenames = {'NREMstate','WAKEstate','REMstate'};

for rr = 1:length(regions)
    disp(['Loading ',regions{rr}])
    %[ISIStats.(regions{rr}),baseNames] = bz_LoadCellinfo(datasetPath.(regions{rr}),'ISIStats','dataset',true,'catall',true);
    %CellClass.(regions{rr}) = bz_LoadCellinfo(datasetPath.(regions{rr}),'CellClass','dataset',true,'catall',true,'baseNames',baseNames);

    [GammaFitAll,baseNames] = bz_LoadAnalysisResults(datasetPath.(regions{rr}),'GammaModeFitAnalysis','dataset',true);
    CellClass.(regions{rr}) = bz_LoadCellinfo(datasetPath.(regions{rr}),'CellClass','dataset',true,'catall',true,'baseNames',baseNames);
    ISIStats.(regions{rr}) = bz_LoadCellinfo(datasetPath.(regions{rr}),'ISIStats','dataset',true,'catall',true,'baseNames',baseNames);
    ISIStats.(regions{rr}) = rmfield( ISIStats.(regions{rr}),'allspikes');
    
    %PopActivityAll = GetMatResults(figfolder,'SpikeStatsbyPopActivityAnalysis','baseNames',baseNames);
    GammaFitAll = bz_CollapseStruct(GammaFitAll);
    
    ISIfits.(regions{rr}) = bz_CollapseStruct(GammaFitAll.ISIfits,'match','justcat',true);
    %Remove cells not in the proper region by removing their cell class!
    if ismember(rr,[4 5])
        inregion = cellfun(@(X) strcmp(X,rnames{rr}),ISIStats.(regions{rr}).cellinfo.regions);
        CellClass.(regions{rr}).label(~inregion)={[]};
        CellClass.(regions{rr}).pE(~inregion)=false;
        CellClass.(regions{rr}).pI(~inregion)=false;
        %clear ISIStats
    end
    
    clear GammaFitAll
end

%%

maxNmodes = 12;

%%
figure
for ss = 1:3
for cc = 1:2
subplot(3,2,(cc)+(ss-1)*2)
%imagesc(log10(fiterror))
hold on
for rr = 1:length(regions)
    plot(1:maxNmodes,nanmean((ISIfits.(regions{rr}).(statenames{ss}).fiterror(CellClass.(regions{rr}).(celltypes{cc}),:)),1),...
        '-o','linewidth',2,'color',regioncolors(rr,:))
   % errorshade(1:maxNmodes,nanmean(log10(ISIfits.(regions{rr}).(statenames{ss}).fiterror(CellClass.(regions{rr}).(celltypes{cc}),:)),1),...
    %    nanstd(log10(ISIfits.(regions{rr}).(statenames{ss}).fiterror(CellClass.(regions{rr}).(celltypes{cc}),:)),[],1),...
    %    nanstd(log10(ISIfits.(regions{rr}).(statenames{ss}).fiterror(CellClass.(regions{rr}).(celltypes{cc}),:)),[],1),regioncolors(rr,:),'scalar');
end
%LogScale('y',10)
xlabel('N Modes');ylabel({statenames{ss},'MSE'})
if ss == 1
    title(celltypes{cc})
    if cc == 1
        legend(regions)
    end
end
end
end





%%
for cc = 1:2
figure
for ss = 1:3
    for rr = 1:length(regions)
        ISIfits.(regions{rr}).(statenames{ss}).weights(ISIfits.(regions{rr}).(statenames{ss}).weights<0.05) = nan;
    ISIfits.(regions{rr}).(statenames{ss}).rates(isnan(ISIfits.(regions{rr}).(statenames{ss}).weights))=nan;

subplot(length(regions),3,(rr-1)*3+ss)
hold on

    plot(log10(ISIfits.(regions{rr}).(statenames{ss}).rates(:,CellClass.(regions{rr}).(celltypes{cc}))),...
        log10(ISIfits.(regions{rr}).(statenames{ss}).CVs(:,CellClass.(regions{rr}).(celltypes{cc}))),...
        '.','color',cellcolor{cc},'markersize',1)
    plot(xlim(gca),[0 0],'k--')
    LogScale('xy',10)
    ylim([-2 1])
    %ylim([0 4])
    xlim([-2 2.5])
    

if rr == 1
title((statenames{ss}))
elseif rr ==length(regions)
    xlabel('Rate (Hz)');
end
if ss == 1
    ylabel({(regions{rr}),'CV'})
end

    end 

end
NiceSave(['AllISImodes',(celltypes{cc})],figfolder,[])
end

%%

figure
for ss = 1:3
    for rr = 1:length(regions)
        subplot(length(regions),3,(rr-1)*3+ss)
            hist(ISIfits.(regions{rr}).(statenames{ss}).Nmodes(CellClass.(regions{rr}).(celltypes{cc})))
    end
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
               
            %Keeping the right region within each recording (i.e. bla/pir separation)    
            if ismember(rr,[4 5]) & ~isempty(abovethresh)
                pairregions.(regions{rr}).(celltypes{cc}){ii,jj} = ...
                    cat(1,whichregion.(regions{rr}).(celltypes{cc}){ii,jj,:});
                inregion = cellfun(@(X) strcmp(X,rnames{rr}),pairregions.(regions{rr}).(celltypes{cc}){ii,jj});
            else
                inregion = true(size(abovethresh));
            end
                
                simmatrices.(regions{rr}).(celltypes{cc})(ii,jj) = ...
                    nanmedian(newpairs.(regions{rr}).(celltypes{cc}){ii,jj}(abovethresh & inregion));
            
                if rr ==1
                    newpairs.ALL.(celltypes{cc}){ii,jj} = ...
                        newpairs.(regions{rr}).(celltypes{cc}){ii,jj}(abovethresh & inregion);
                else
                    newpairs.ALL.(celltypes{cc}){ii,jj} = ...
                        [newpairs.ALL.(celltypes{cc}){ii,jj};newpairs.(regions{rr}).(celltypes{cc}){ii,jj}(abovethresh & inregion)];

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
            
            if ismember(rr,[4 5]) & ~isempty(abovethresh)
                pairregions.(regions{rr}).difft{ii,jj} = ...
                    cat(1,whichregion.(regions{rr}).difft{ii,jj,:});
                inregion = cellfun(@(X) strcmp(X,rnames{rr}),pairregions.(regions{rr}).difft{ii,jj});
            else
                inregion = true(size(abovethresh));
            end
            
            simmatrices.(regions{rr}).difft(ii,jj) = nanmedian(newpairs.(regions{rr}).difft{ii,jj}(abovethresh & inregion));
            
            if rr ==1
                newpairs.ALL.difft{ii,jj} = newpairs.(regions{rr}).difft{ii,jj}(abovethresh & inregion);
            else
                newpairs.ALL.difft{ii,jj} = [newpairs.ALL.difft{ii,jj};newpairs.(regions{rr}).difft{ii,jj}(abovethresh & inregion)];
                
                simmatrices.ALL.difft(ii,jj) = nanmedian(newpairs.ALL.difft{ii,jj});
            end
        end
    end
end


%%
figure
for rr = 1:length(regions)

for cc = 1:length(celltypes)
subplot(4,6,(cc-1).*6+rr)
imagesc(simmatrices.(regions{rr}).(celltypes{cc}))
alpha(gca,single(~isnan(simmatrices.(regions{rr}).(celltypes{cc}))))
caxis([0.1 0.35])
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
caxis([0.1 0.35])
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
for rr = 1:length(regions)
subplot(2,3,rr)
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
subplot(2,2,1)
imagesc(simmatrices.ALL.difft)
alpha(gca,single(~isnan(simmatrices.ALL.difft)))
colorbar
caxis([0.5 2.5])
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
keepISIs.numPerregion = 2000;
keepISIs.numISIthresh = 800;
keepISIs.numISIthresh = 300; %Try with smoothing
for rr = 1:length(regions)
    OKISIS = (allISIs.(regions{rr}).celltype.pE | allISIs.(regions{rr}).celltype.pI) & ...
        ~any(isnan(allISIs.(regions{rr}).hists),1) &  allISIs.(regions{rr}).numISIs>keepISIs.numISIthresh;
    %Here: inregion
    if ismember(rr,[4 5])
        INREGION = cellfun(@(X) strcmp(X,rnames{rr}),allISIs.(regions{rr}).region);
    else
        INREGION = true(size(OKISIS));
    end
    if length(find(OKISIS & INREGION))>keepISIs.numPerregion
        keepISIs.(regions{rr}) = randsample(find(OKISIS & INREGION), keepISIs.numPerregion);
    else
        keepISIs.(regions{rr}) = find(OKISIS & INREGION);
    end
end
% allISIs = [ISIs.THAL(keepISIs.THAL),ISIs.vCTX(keepISIs.vCTX),...
%     ISIs.fCTX(keepISIs.fCTX),ISIs.CA1(keepISIs.CA1)];
%numISIs = cellfun(@length,allISIs);
allISIhists.hists = [allISIs.THAL.hists(:,keepISIs.THAL),allISIs.vCTX.hists(:,keepISIs.vCTX),...
     allISIs.fCTX.hists(:,keepISIs.fCTX),allISIs.CA1.hists(:,keepISIs.CA1),...
     allISIs.BLA.hists(:,keepISIs.BLA),allISIs.PIR.hists(:,keepISIs.PIR)];
 
 %%
ALLregions = [1.*ones(size(keepISIs.THAL)),2.*ones(size(keepISIs.vCTX)),...
    3.*ones(size(keepISIs.fCTX)),6.*ones(size(keepISIs.CA1)),...
    4.*ones(size(keepISIs.BLA)),5.*ones(size(keepISIs.PIR))];
ALLcelltypes.pI = [allISIs.THAL.celltype.pI(keepISIs.THAL),allISIs.vCTX.celltype.pI(keepISIs.vCTX),...
    allISIs.fCTX.celltype.pI(keepISIs.fCTX),allISIs.CA1.celltype.pI(keepISIs.CA1),...
    allISIs.BLA.celltype.pI(keepISIs.BLA),allISIs.PIR.celltype.pI(keepISIs.PIR)];
ALLcelltypes.pE = [allISIs.THAL.celltype.pE(keepISIs.THAL),allISIs.vCTX.celltype.pE(keepISIs.vCTX),...
    allISIs.fCTX.celltype.pE(keepISIs.fCTX),allISIs.CA1.celltype.pE(keepISIs.CA1),...
    allISIs.BLA.celltype.pE(keepISIs.BLA),allISIs.PIR.celltype.pE(keepISIs.PIR)];
states.ALL = [allISIs.THAL.state(keepISIs.THAL),allISIs.vCTX.state(keepISIs.vCTX),...
    allISIs.fCTX.state(keepISIs.fCTX),allISIs.CA1.state(keepISIs.CA1),...
    allISIs.BLA.state(keepISIs.BLA),allISIs.PIR.state(keepISIs.PIR)];

[jregions,iregions] = meshgrid(ALLregions,ALLregions);
[jiscelltype.pE,iiscelltype.pE] = meshgrid(ALLcelltypes.pE,ALLcelltypes.pE);
[jiscelltype.pI,iiscelltype.pI] = meshgrid(ALLcelltypes.pI,ALLcelltypes.pI);
[jstates,istates] = meshgrid(states.ALL,states.ALL);

%% Try denoising ISIdists....

nvisbins = 20;
[Xdn, Sigma, npars] = bz_PCAdenoise(allISIhists.hists', nvisbins,'percvar',0.99);

%%
ex = randsample(size(allISIhists.hists,2),15);
figure
for ss = 1:15
subplot(3,5,ss)
plot(allISIhists.hists(:,ex(ss)))
hold on
plot(Xdn(ex(ss),:))
end
%%
figure
imagesc(Xdn)
%%
allISIhists.hists = Xdn';
allISIhists.hists(allISIhists.hists<0) = 0;

for cc = 1:length(allISIhists.hists)
    allISIhists.hists(:,cc) = allISIhists.hists(:,cc)./sum(allISIhists.hists(:,cc));
end

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
        simmatrices.(statenames{ss}).(celltypes{cc}) = nan(length(regions));
        simmatrices.ALLregions.(celltypes{cc}) = nan(length(regions));
    end
end
simmatrices.ALLregions.(celltypes{cc}) = nan(length(regions));
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
                    nanmedian(allpairs.ALLregions.(celltypes{cc}){rr1,rr2});
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
            caxis([0.75 1.75])
            if cc == 1
            title(statenames{ss})
            end
            box off
            crameri tokyo
            set(gca,'ytick',[1:length(regions)]);set(gca,'xtick',[1:length(regions)])
            set(gca,'yticklabels',regions);set(gca,'xticklabels',regions)
    end
        subplot(3,3,6+cc)
            imagesc(simmatrices.ALLregions.(celltypes{cc}))
            %colorbar
            alpha(gca,single(~isnan(simmatrices.ALLregions.(celltypes{cc}))))
            colorbar
            caxis([0.75 1.75])
            if cc == 1
            title('ALL')
            end
            box off
            crameri tokyo
            set(gca,'ytick',[1:length(regions)]);set(gca,'xtick',[1:length(regions)])
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
perps = 10:10:50;
for pp = 1:length(perps)
close all

P = d2p(valid, perps(pp), 1e-6); 
Y(:,:,pp) = tsne_p(P, ALLcelltypes.pE, 2, 1000);

close all
P = d2p(valid.^2, perps(pp), 1e-6); 
Ysq(:,:,pp) = tsne_p(P, ALLcelltypes.pE, 2, 1000);
end
%% 2D
figure
for pp = 1:length(perps)
%legend('NREM','WAKE','REM')
subplot(4,5,pp)
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
subplot(4,5,pp+5)
hold on
axis tight
box on
xticks(gca,[]);yticks(gca,[])
plot(Y(ALLcelltypes.pE==1,1,pp),Y(ALLcelltypes.pE==1,2,pp),'.k','markersize',3)
plot(Y(ALLcelltypes.pI==1,1,pp),Y(ALLcelltypes.pI==1,2,pp),'.r','markersize',3)
%legend(celltypes,'location','northoutside')


subplot(4,5,pp+10)
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
subplot(4,5,pp+15)
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
perplexity = 30;
P = d2p(valid.^2, perplexity, 1e-8); 
tSNEmap(:,:) = tsne_p(P, ALLcelltypes.pE, 2, 2000);
%%

figure
subplot(2,2,1)
hold on
plot(tSNEmap(states.ALL==3,1),tSNEmap(states.ALL==3,2),'r.','markersize',3)
plot(tSNEmap(states.ALL==1,1),tSNEmap(states.ALL==1,2),'b.','markersize',3)
plot(tSNEmap(states.ALL==2,1),tSNEmap(states.ALL==2,2),'k.','markersize',3)
axis tight
box on
xticks(gca,[]);yticks(gca,[])
legend(statenames,'location','northwest')

%legend('NREM','WAKE','REM')
subplot(2,2,2)
hold on
for rr = 1:length(regions)
    plot(tSNEmap(ALLregions==rr,1),tSNEmap(ALLregions==rr,2),...
        '.','markersize',3,'color',regioncolors(rr,:))
end
plot(tSNEmap(ALLcelltypes.pI==1,1),tSNEmap(ALLcelltypes.pI==1,2),'.r','markersize',2)
legend([regions,'pI'],'location','northwest')
axis tight
box on
xticks(gca,[]);yticks(gca,[])

% subplot(3,3,1)
% hold on
% axis tight
% box on
% xticks(gca,[]);yticks(gca,[])
% plot(tSNEmap(ALLcelltypes.pE==1,1),tSNEmap(ALLcelltypes.pE==1,2),'.k','markersize',3)
% plot(tSNEmap(ALLcelltypes.pI==1,1),tSNEmap(ALLcelltypes.pI==1,2),'.r','markersize',3)
% legend(celltypes,'location','northwest')
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
        statetSNEmap.(statenames{ss})(ALLregions(states.ALL==ss)==rr,2),...
        '.','markersize',4,'color',regioncolors(rr,:))
end
plot(statetSNEmap.(statenames{ss})(ALLcelltypes.pI(states.ALL==ss)==1,1),...
    statetSNEmap.(statenames{ss})(ALLcelltypes.pI(states.ALL==ss)==1,2),'.r','markersize',4)
legend([regions,'pI'],'location','northwest')
axis tight
box on
xticks(gca,[]);yticks(gca,[])

% subplot(2,3,ss+3)
% hold on
% axis tight
% 
% xticks(gca,[]);yticks(gca,[])
% plot(statetSNEmap.(statenames{ss})(ALLcelltypes.pE(states.ALL==ss)==1,1),...
%     statetSNEmap.(statenames{ss})(ALLcelltypes.pE(states.ALL==ss)==1,2),'.k','markersize',4)
% plot(statetSNEmap.(statenames{ss})(ALLcelltypes.pI(states.ALL==ss)==1,1),...
%     statetSNEmap.(statenames{ss})(ALLcelltypes.pI(states.ALL==ss)==1,2),'.r','markersize',4)
% legend(celltypes,'location','northeast')
% box on
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
            statetSNEmap.pE.(statenames{ss})(ALLregions(ALLcelltypes.pE==1 & states.ALL==ss)==rr,2),...
            '.','markersize',4,'color',regioncolors(rr,:))
    end
    axis tight
    box on
    xticks(gca,[]);yticks(gca,[])
    title(statenames{ss})
end

subplot(2,2,3)
hold on
plot(statetSNEmap.pE.all(states.ALL(ALLcelltypes.pE==1)==2,1),statetSNEmap.pE.all(states.ALL(ALLcelltypes.pE==1)==2,2),...
    'k.','markersize',4)
plot(statetSNEmap.pE.all(states.ALL(ALLcelltypes.pE==1)==1,1),statetSNEmap.pE.all(states.ALL(ALLcelltypes.pE==1)==1,2),...
    'b.','markersize',4)
plot(statetSNEmap.pE.all(states.ALL(ALLcelltypes.pE==1)==3,1),statetSNEmap.pE.all(states.ALL(ALLcelltypes.pE==1)==3,2),...
    'r.','markersize',4)
axis tight
box on
xticks(gca,[]);yticks(gca,[])
legend(statenames,'location','northwest')

%legend('NREM','WAKE','REM')
subplot(2,2,4)
hold on
for rr = 1:length(regions)
    plot(statetSNEmap.pE.all(ALLregions(ALLcelltypes.pE==1)==rr,1),statetSNEmap.pE.all(ALLregions(ALLcelltypes.pE==1)==rr,2),...
        '.','markersize',4,'color',regioncolors(rr,:))
end
legend(regions,'location','northwest')
axis tight
box on
xticks(gca,[]);yticks(gca,[])

NiceSave('tSNE_Map_pEOnly',figfolder,'')

