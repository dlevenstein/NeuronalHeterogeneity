function [ISIs,ispE,ispI,UIDs,states,rates,numISIs,allpairs ] = ISIDistSimilarityAnalysis(basePath,figfolder)
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



%% Restrict to state

for ss =1:3
ISIStats.allspikes.instate = cellfun(@(X) InIntervals(X,double(SleepState.ints.(states{ss}))),...
    ISIStats.allspikes.times,'UniformOutput',false);

AllStateISIs.(states{ss}) = cellfun(@(X,Y) log10(X(Y)),...
    ISIStats.allspikes.ISIs,ISIStats.allspikes.instate,...
    'UniformOutput',false);
end


%%
ISIs = [AllStateISIs.WAKEstate,AllStateISIs.NREMstate,AllStateISIs.REMstate];
numISIs = cellfun(@length,ISIs);
%%

numcells = length(ISIs);
KSSTAT = nan(numcells);
%parpool
for ii = 1:numcells
    bz_Counter(ii,numcells,'Cell');
    for jj = ii:numcells
        if ii==jj %For Same Cell, calculate first/last half spikes
            [~,~,KSSTAT(ii,jj)] = kstest2(ISIs{ii}(1:round(end/2)),ISIs{jj}(round(end/2):end));
        else
            [~,~,KSSTAT(ii,jj)] = kstest2(ISIs{ii},ISIs{jj});
        end
    end
    for jj = ii:numcells
        KSSTAT(jj,ii) = KSSTAT(ii,jj);
    end
end
%clear ISIs

%%
states = [2.*ones(size(CellClass.pE)),1.*ones(size(CellClass.pE)),3.*ones(size(CellClass.pE))];
ispE = [CellClass.pE,CellClass.pE,CellClass.pE];
ispI = [CellClass.pI,CellClass.pI,CellClass.pI];
UIDs = [ISIStats.UID,ISIStats.UID,ISIStats.UID];
[iUID,jUID] = meshgrid(UIDs,UIDs);
[iispE,jispE] = meshgrid(ispE,ispE);
[iispI,jispI] = meshgrid(ispI,ispI);
[istates,jstates] = meshgrid(states,states);
samecell = iUID==jUID;

simmatrices.pE = nan(3);
simmatrices.pI = nan(3);
simmatrices.difft = nan(6);


for ss1 = 1:3
    for ss2 = ss1:-1:1
        allpairs.pE{ss1,ss2} = KSSTAT(samecell & istates==ss1 & jstates==ss2 & iispE);
        allpairs.pI{ss1,ss2} = KSSTAT(samecell & istates==ss1 & jstates==ss2 & iispI);
        
        simmatrices.pE(ss1,ss2) = mean(allpairs.pE{ss1,ss2});
        simmatrices.pI(ss1,ss2) = mean(allpairs.pI{ss1,ss2});
    
        
        allpairs.difft{ss1,ss2} = KSSTAT(~samecell & istates==ss1 & jstates==ss2 & jispE & iispE);
        allpairs.difft{ss1+3,ss2} = KSSTAT(~samecell & istates==ss1 & jstates==ss2 & iispI & jispE);
        allpairs.difft{ss1+3,ss2+3} = KSSTAT(~samecell & istates==ss1 & jstates==ss2 & iispI & jispI);
        
        simmatrices.difft(ss1,ss2) = mean(allpairs.difft{ss1,ss2});
        simmatrices.difft(ss1+3,ss2) = mean(allpairs.difft{ss1+3,ss2});
        simmatrices.difft(ss1+3,ss2+3) = mean(allpairs.difft{ss1+3,ss2+3});
    end
end


%% MDS Map
valid = KSSTAT;
for ii =1:numcells
        valid(ii,ii)=0;
end
% Y = cmdscale(valid);

%% tSNE
perplexity = 30;
P = d2p(valid, perplexity, 1e-5); 
Y = tsne_p(P, [], 2);
%%

rates = [ISIStats.summstats.WAKEstate.meanrate,...
    ISIStats.summstats.NREMstate.meanrate,ISIStats.summstats.REMstate.meanrate];
%%
figure
for cc = 1:length(celltypes)
subplot(3,3,(cc-1).*3+1)
imagesc(simmatrices.(celltypes{cc}))
alpha(gca,single(~isnan(simmatrices.(celltypes{cc}))))
caxis([0.075 0.35])
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
caxis([0.075 0.35])
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
plot(Y(states==1,1),Y(states==1,2),'b.')
plot(Y(states==2,1),Y(states==2,2),'k.')
plot(Y(states==3,1),Y(states==3,2),'r.')
axis tight
xticks(gca,[]);yticks(gca,[])
%legend('NREM','WAKE','REM')
subplot(3,3,9)
scatter(Y(:,1),Y(:,2),5,log10(rates))
axis tight
xticks(gca,[]);yticks(gca,[])
subplot(3,3,7)
hold on
axis tight
xticks(gca,[]);yticks(gca,[])
plot(Y(ispE,1),Y(ispE,2),'.k')
plot(Y(ispI,1),Y(ispI,2),'.r')
NiceSave('ISISimilarity',figfolder,baseName)
