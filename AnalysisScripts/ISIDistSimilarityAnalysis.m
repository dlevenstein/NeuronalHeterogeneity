function [allISIs,allpairs,lowestpairISI ] = ISIDistSimilarityAnalysis(basePath,figfolder)
% Date XX/XX/20XX
%
%Question: 
%
%Plots
%-
%-
%
%% Load Header
%Initiate Paths
%reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
%reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/Cicero_09102014';
%basePath = pwd;
%basePath = fullfile(reporoot,'Datasets/onProbox/AG_HPC/Achilles_11012013');
%figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%Load Stuff
sessionInfo = bz_getSessionInfo(basePath);
%spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
states = fieldnames(SleepState.ints);
states{4} = 'ALL';
SleepState.ints.ALL = [0 Inf];
statecolors = {'k','b','r',[0.6 0.6 0.6]};

try
    celltypes = CellClass.celltypes;
catch
    celltypes = unique(CellClass.label);
end
cellcolor = {'k','r'};



%% Restrict ISIs to state

for ss =1:3
ISIStats.allspikes.instate = cellfun(@(X) InIntervals(X,double(SleepState.ints.(states{ss}))),...
    ISIStats.allspikes.times,'UniformOutput',false);

AllStateISIs.(states{ss}) = cellfun(@(X,Y) log10(X(Y)),...
    ISIStats.allspikes.ISIs,ISIStats.allspikes.instate,...
    'UniformOutput',false);
end


%% Compile all states together
allISIs.ISIs = [AllStateISIs.WAKEstate,AllStateISIs.NREMstate,AllStateISIs.REMstate];
allISIs.numISIs = cellfun(@length,allISIs.ISIs);

allISIs.state = [2.*ones(size(CellClass.pE)),1.*ones(size(CellClass.pE)),3.*ones(size(CellClass.pE))];
allISIs.celltype.pE = [CellClass.pE,CellClass.pE,CellClass.pE];
allISIs.celltype.pI = [CellClass.pI,CellClass.pI,CellClass.pI];
allISIs.UIDs = [ISIStats.UID,ISIStats.UID,ISIStats.UID];

allISIs.hists = [ISIStats.ISIhist.WAKEstate.log',ISIStats.ISIhist.NREMstate.log',...
    ISIStats.ISIhist.REMstate.log'];
allISIs.histbins = ISIStats.ISIhist.logbins;
%% Calculate KS and KL distance between ISIs, all cells

numcells = length(allISIs.ISIs);
KSSTAT = nan(numcells);
KLDIST = nan(numcells);
%parpool
for ii = 1:numcells
    bz_Counter(ii,numcells,'Cell');
    KLDIST(:,ii) = KLDiv(allISIs.hists',allISIs.hists(:,ii)','symmetric',true);
%     parfor jj = ii:numcells
%         if ii==jj %For Same Cell, calculate first/last half spikes
%             [~,~,KSSTAT(ii,jj)] = kstest2(allISIs.ISIs{ii}(1:round(end/2)),allISIs.ISIs{jj}(round(end/2):end));
            for jj=ii
                if isempty(allISIs.ISIs{ii})
                    KLDIST(ii,jj) = nan;
                    continue
                end
            firsthalf = hist(allISIs.ISIs{ii}(1:round(end/2)),allISIs.histbins);
            secondhalf = hist(allISIs.ISIs{jj}(round(end/2):end),allISIs.histbins);
            
            KLDIST(ii,jj) = KLDiv(firsthalf,secondhalf);
            end
%         else
%             [~,~,KSSTAT(ii,jj)] = kstest2(allISIs.ISIs{ii},allISIs.ISIs{jj});
%         end
%     end
%     for jj = ii:numcells
%         KSSTAT(jj,ii) = KSSTAT(ii,jj);
%     end
end
%clear ISIs

%%
figure
subplot(2,2,1)
imagesc(KLDIST)
subplot(2,2,2)
%imagesc(KSSTAT)
%%
[jUID,iUID] = meshgrid(allISIs.UIDs,allISIs.UIDs);
[jispE,iispE] = meshgrid(allISIs.celltype.pE,allISIs.celltype.pE);
[jispI,iispI] = meshgrid(allISIs.celltype.pI,allISIs.celltype.pI);
[jstates,istates] = meshgrid(allISIs.state,allISIs.state);
[jnumISI,inumISI] = meshgrid(allISIs.numISIs,allISIs.numISIs);
samecell = iUID==jUID;

simmatrices.pE = nan(3);
simmatrices.pI = nan(3);
simmatrices.difft = nan(6);

%%
% figure
% BoxAndScatterPlot(allpairs.pE(:))
% 
% %%
% figure
% imagesc((samecell & istates==3 & jstates==3 & iispE))
% 
% %%
% ss1 = 3;
% ss2 = 1;
% figure
% plot(lowestpairISI.pI{ss1,ss2},allpairs.pI{ss1,ss2},'.')
%%
numISIthresh = 1000;
for ss1 = 1:3
    for ss2 = ss1:-1:1
        allpairs.pE{ss1,ss2} = KLDIST(samecell & iispE & istates==ss1 & jstates==ss2);
        allpairs.pI{ss1,ss2} = KLDIST(samecell & iispI & istates==ss1 & jstates==ss2);
        
        lowestpairISI.pE{ss1,ss2} = min([inumISI(samecell & istates==ss1 & jstates==ss2 & iispE),...
            jnumISI(samecell & istates==ss1 & jstates==ss2 & iispE)],[],2);
        lowestpairISI.pI{ss1,ss2} = min([inumISI(samecell & istates==ss1 & jstates==ss2 & iispI),...
            jnumISI(samecell & istates==ss1 & jstates==ss2 & iispI)],[],2);
        
        simmatrices.pE(ss1,ss2) = median(allpairs.pE{ss1,ss2}(lowestpairISI.pE{ss1,ss2}>numISIthresh));
        simmatrices.pI(ss1,ss2) = median(allpairs.pI{ss1,ss2}(lowestpairISI.pI{ss1,ss2}>numISIthresh));
    
        
        allpairs.difft{ss1,ss2} = KLDIST(~samecell & iispE & jispE);
        allpairs.difft{ss1+3,ss2} = KLDIST(~samecell & istates==ss1 & jstates==ss2 & iispI & jispE);
        allpairs.difft{ss1+3,ss2+3} = KLDIST(~samecell & istates==ss1 & jstates==ss2 & iispI & jispI);
        
        lowestpairISI.difft{ss1,ss2} = min([inumISI(~samecell & istates==ss1 & jstates==ss2 & jispE & iispE),...
            jnumISI(~samecell & istates==ss1 & jstates==ss2 & jispE & iispE)],[],2);
        lowestpairISI.difft{ss1+3,ss2} = min([inumISI(~samecell & istates==ss1 & jstates==ss2 & iispI & jispE),...
            jnumISI(~samecell & istates==ss1 & jstates==ss2 & iispI & jispE)],[],2);
        lowestpairISI.difft{ss1+3,ss2+3} = min([inumISI(~samecell & istates==ss1 & jstates==ss2 & iispI & jispI),...
            jnumISI(~samecell & istates==ss1 & jstates==ss2 & iispI & jispI)],[],2);

        
        simmatrices.difft(ss1,ss2) = median(allpairs.difft{ss1,ss2}(lowestpairISI.difft{ss1,ss2}>numISIthresh));
        simmatrices.difft(ss1+3,ss2) = median(allpairs.difft{ss1+3,ss2}(lowestpairISI.difft{ss1+3,ss2}>numISIthresh));
        simmatrices.difft(ss1+3,ss2+3) = median(allpairs.difft{ss1+3,ss2+3}(lowestpairISI.difft{ss1+3,ss2+3}>numISIthresh));
    end
end

%%
%simmatrices.difft(3,1) = 10;

%% MDS Map
valid = KLDIST;
for ii =1:numcells
        valid(ii,ii)=0;
end
valid(isnan(valid))=0;
% Y = cmdscale(valid);

%% tSNE
close all
perplexity = 20;
P = d2p(valid, perplexity, 1e-8); 
Y = tsne_p(P, [], 2);
%%

allISIs.rates = [ISIStats.summstats.WAKEstate.meanrate,...
    ISIStats.summstats.NREMstate.meanrate,ISIStats.summstats.REMstate.meanrate];
%%
figure
for cc = 1:length(celltypes)
subplot(3,3,(cc-1).*3+1)
imagesc(simmatrices.(celltypes{cc}))
alpha(gca,single(~isnan(simmatrices.(celltypes{cc}))))
caxis([0 3])
%colorbar
box off
set(gca,'ytick',[1:3]);set(gca,'xtick',[1:3]);
set(gca,'yticklabel',{'N','W','R'})
set(gca,'xticklabel',{'N','W','R'})
ylabel(celltypes{cc});
crameri tokyo
if cc == 1
    title('Same Cell')
end
end

subplot(3,3,[2 3 5 6])
imagesc(simmatrices.difft)
alpha(gca,single(~isnan(simmatrices.difft)))
colorbar
caxis([0 3])
%caxis([0.075 0.35])
title('Diff. Cells')
box off
crameri tokyo
set(gca,'ytick',[1:6]);set(gca,'xtick',[1:6]);
set(gca,'yticklabel',{'N','W','R','N','W','R'})
set(gca,'xticklabel',{'N','W','R','N','W','R'})

% subplot(2,2,1)
% imagesc(KSSTAT)
subplot(3,3,8)
hold on
plot(Y(allISIs.state==1,1),Y(allISIs.state==1,2),'b.')
plot(Y(allISIs.state==2,1),Y(allISIs.state==2,2),'k.')
plot(Y(allISIs.state==3,1),Y(allISIs.state==3,2),'r.')
axis tight
xticks(gca,[]);yticks(gca,[])
%legend('NREM','WAKE','REM')
subplot(3,3,9)
scatter(Y(:,1),Y(:,2),5,log10(allISIs.rates))
axis tight
xticks(gca,[]);yticks(gca,[])
subplot(3,3,7)
hold on
axis tight
xticks(gca,[]);yticks(gca,[])
plot(Y(allISIs.celltype.pE,1),Y(allISIs.celltype.pE,2),'.k')
plot(Y(allISIs.celltype.pI,1),Y(allISIs.celltype.pI,2),'.r')
NiceSave('ISISimilarity',figfolder,baseName)
