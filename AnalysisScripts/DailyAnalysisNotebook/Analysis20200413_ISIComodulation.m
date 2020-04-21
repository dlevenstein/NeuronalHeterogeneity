function [ ] = Analysis20200413_ISIComodulation(basePath,figfolder)
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
reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/';
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
basePath = '/Users/dl2820/Dropbox/Research/Datasets/Cicero_09102014';
%basePath = pwd;
%basePath = fullfile(reporoot,'Datasets/onProbox/AG_HPC/Achilles_11012013');
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%Load Stuff
sessionInfo = bz_getSessionInfo(basePath);
spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
states = fieldnames(SleepState.ints);
%states{4} = 'ALL';
%SleepState.ints.ALL = [0 Inf];
statecolors = {'k','b','r'};

try
    celltypes = CellClass.celltypes;
catch
    celltypes = unique(CellClass.label);
end
cellcolor = {'k','r'};

%%
ISIStats.allspikes.ISInp1 = cellfun(@(X) [X(2:end);nan],ISIStats.allspikes.ISIs,...
    'UniformOutput',false);

%%

twin = 0.05; %(s)
ISItol = 0.2; %log(s)
numspikethresh = 30; %Need this many spikes to calculate hist (i or j), otherwise nan.

numISIbins = 100;
logISIbins = linspace(-2.5,1.5,numISIbins);
smoothwin = ISItol./mode(diff(logISIbins));
%% Restrict to state
clear ISIcomod
for ss = 1:2
state = states{ss};
ISIStats.allspikes.instate = cellfun(@(X) InIntervals(X,double(SleepState.ints.(state))),...
    ISIStats.allspikes.times,'UniformOutput',false);

%%
usespikes.times = cellfun(@(X,Y) X(Y), ISIStats.allspikes.times,ISIStats.allspikes.instate,'UniformOutput',false);
usespikes.ISIs = cellfun(@(X,Y) X(Y), ISIStats.allspikes.ISIs,ISIStats.allspikes.instate,'UniformOutput',false);
usespikes.ISInp1 = cellfun(@(X,Y) X(Y), ISIStats.allspikes.ISInp1,ISIStats.allspikes.instate,'UniformOutput',false);




i_ISIhist_marj = nan(length(logISIbins),spikes.numcells);
i_ISIhist = nan(length(logISIbins),length(logISIbins),spikes.numcells,spikes.numcells);
i_ISIhist_0bin = nan(length(logISIbins),spikes.numcells,spikes.numcells); 
i_ISIhist_any = nan(length(logISIbins),spikes.numcells,spikes.numcells); 

normS_i_ISIhist = nan(length(logISIbins),length(logISIbins),spikes.numcells,spikes.numcells);
norm_i_ISIhist = nan(length(logISIbins),length(logISIbins),spikes.numcells,spikes.numcells);
norm_i_ISIhist_0bin = nan(length(logISIbins),spikes.numcells,spikes.numcells); 
norm_i_ISIhist_any = nan(length(logISIbins),spikes.numcells,spikes.numcells); 

%DEV
%%
for ii = 1
    bz_Counter(ii,spikes.numcells,'Cell i')
    %
    
    %Get the marginal ISI distribution for cell i
    [i_ISIhist_marj(:,ii)] = hist(log10([usespikes.ISIs{ii};usespikes.ISInp1{ii}]),...
        logISIbins);
    %notenough = i_ISIhist_marj<numspikethresh ;
    i_ISIhist_marj(:,ii) = movmean(i_ISIhist_marj(:,ii),smoothwin);
    i_ISIhist_marj(:,ii) = i_ISIhist_marj(:,ii)./nansum(i_ISIhist_marj(:,ii)); %Normalize
    %i_ISIhist_marj(
    %i_ISIhist_marj(notenough) = nan;
    
    
    %Get the ISI dist conditioned on spikes/no spikes
    % Here: find i spikes that aren't within tolerance of any j spikes
    cellispikes_0bin = cellfun(@(times_j) ~ismembertol(usespikes.times{ii},...
        times_j,twin,'DataScale', 1),...
        usespikes.times,'UniformOutput',false);

    temp = cellfun(@(countspikes) hist(log10([usespikes.ISIs{ii}(countspikes);...
        usespikes.ISInp1{ii}(countspikes)]),...
        logISIbins),...
        cellispikes_0bin,'UniformOutput',false);
	temp = cellfun(@(histo) movmean(histo,smoothwin),temp,'UniformOutput',false); %Normalize
    temp = cellfun(@(histo) histo./sum(histo),temp,'UniformOutput',false); %Normalize
    temp = cat(1,temp{:});
    i_ISIhist_0bin(:,:,ii) = temp';%./i_ISIhist_marj';
    norm_i_ISIhist_0bin(:,:,ii) = temp'-i_ISIhist_marj(:,ii);

    temp = cellfun(@(countspikes) hist(log10([usespikes.ISIs{ii}(~countspikes);...
        usespikes.ISInp1{ii}(~countspikes)]),...
        logISIbins),...
        cellispikes_0bin,'UniformOutput',false);
	temp = cellfun(@(histo) movmean(histo,smoothwin),temp,'UniformOutput',false); %Normalize
    temp = cellfun(@(histo) histo./sum(histo),temp,'UniformOutput',false); %Normalize
    temp = cat(1,temp{:});
    i_ISIhist_any(:,:,ii) = temp';%./i_ISIhist_marj';
    norm_i_ISIhist_any(:,:,ii) = temp'-i_ISIhist_marj(:,ii);
    
    
    %tic
    parfor bb = 1:numISIbins %For each ISI bin,
        %bb
        %find cell j spikes in that bin.   
        jbinspikes = cellfun(@(ISIn,ISInp1) ismembertol(log10(ISIn),logISIbins(bb),0.5.*ISItol,'DataScale', 1)...
             |ismembertol(log10(ISInp1),logISIbins(bb),0.5.*ISItol,'DataScale', 1),...
             usespikes.ISIs,usespikes.ISInp1,'UniformOutput',false);
        
        numjspikes = cellfun(@sum,jbinspikes);
        
        %Find all cell i spikes within twin of those spikes
        ibinspikes = cellfun(@(times_j,inbin) ismembertol(usespikes.times{ii},...
            times_j(inbin),twin,'DataScale', 1),...
            usespikes.times,jbinspikes,'UniformOutput',false);
        
        numispikes = cellfun(@sum,ibinspikes);

        %and histogram of their ISIs (prev/next)
        temp = ...
            cellfun(@(nearspikes) hist(log10([usespikes.ISIs{ii}(nearspikes);usespikes.ISInp1{ii}(nearspikes)]),...
            logISIbins),...
            ibinspikes,'UniformOutput',false);
        temp = cellfun(@(histo) movmean(histo,smoothwin),temp,'UniformOutput',false); %Normalize
        temp = cellfun(@(histo) histo./sum(histo),temp,'UniformOutput',false); %Normalize
        temp = cat(1,temp{:});
        temp = temp';
        %temp = temp'./i_ISIhist_marj'; %Compare to marginal
        %This is the matrix for all ISIs of cell i with all cells j in a
        %given ISI_j bin 
        
        %remove rows with not enough spikes (For j, not enough
        %spikes in bin), for i, not enough spikes close to a j spike
        temp(:,numjspikes<numspikethresh | numispikes<numspikethresh) = nan;
        
        i_ISIhist(bb,:,:,ii) = temp;
        %norm_i_ISIhist(bb,:,:,ii) = temp-i_ISIhist_marj(:,ii);
        normS_i_ISIhist(bb,:,:,ii) = temp-i_ISIhist_any(:,:,ii);
        norm_i_ISIhist(bb,:,:,ii) = temp-i_ISIhist_marj(:,ii);

    end


end
    %toc
    %%


%%
ISIcomod.(state).ISIpars = i_ISIhist;
ISIcomod.(state).anyspike = i_ISIhist_any;
ISIcomod.(state).nospike = i_ISIhist_0bin;
ISIcomod.(state).i_ISIhist_marj = i_ISIhist_marj;

ISIcomod.(state).norm.ISIpars = norm_i_ISIhist;
ISIcomod.(state).norm.ISIpars_spknorm = normS_i_ISIhist;
ISIcomod.(state).norm.anyspike = norm_i_ISIhist_any;
ISIcomod.(state).norm.nospike = norm_i_ISIhist_0bin;
end
%%

for ss = 1:2
state = states{ss};
%clear ISIcomod
for tt = 1:2
    for cc = 1:spikes.numcells
        othercells = CellClass.(celltypes{tt});
        othercells(cc) = false;
        ISIcomod.(state).(celltypes{tt})(:,:,cc) = nanmean(ISIcomod.(state).ISIpars(:,:,othercells,cc),3);
        ISIcomod.(state).zbin.(celltypes{tt})(:,cc) = nanmean(ISIcomod.(state).nospike(:,othercells,cc),2);
        ISIcomod.(state).abin.(celltypes{tt})(:,cc) = nanmean(ISIcomod.(state).anyspike(:,othercells,cc),2);
        
        ISIcomod.(state).norm.(celltypes{tt})(:,:,cc) = nanmean(ISIcomod.(state).norm.ISIpars(:,:,othercells,cc),3);
        ISIcomod.(state).norm.zbin.(celltypes{tt})(:,cc) = nanmean(ISIcomod.(state).norm.nospike(:,othercells,cc),2);
        ISIcomod.(state).norm.abin.(celltypes{tt})(:,cc) = nanmean(ISIcomod.(state).norm.anyspike(:,othercells,cc),2);
    end
        meanISIcomod.(state).norm.(celltypes{tt}) = nanmean(ISIcomod.(state).norm.(celltypes{tt})(:,:,CellClass.(celltypes{tt})),3);
        meanISIcomod.(state).norm.zbin.(celltypes{tt}) = nanmean(ISIcomod.(state).norm.zbin.(celltypes{tt})(:,CellClass.(celltypes{tt})),2);
        meanISIcomod.(state).norm.abin.(celltypes{tt}) = nanmean(ISIcomod.(state).norm.abin.(celltypes{tt})(:,CellClass.(celltypes{tt})),2);
        
        meanISIcomod.(state).(celltypes{tt}) = nanmean(ISIcomod.(state).(celltypes{tt})(:,:,CellClass.(celltypes{tt})),3);
        meanISIcomod.(state).zbin.(celltypes{tt}) = nanmean(ISIcomod.(state).zbin.(celltypes{tt})(:,CellClass.(celltypes{tt})),2);
        meanISIcomod.(state).abin.(celltypes{tt}) = nanmean(ISIcomod.(state).abin.(celltypes{tt})(:,CellClass.(celltypes{tt})),2);
end
end

%%
for ss = 1:2
    state = states{ss};
tt = 1;
figure
subplot(2,2,1)
imagesc(logISIbins,logISIbins,(meanISIcomod.(state).(celltypes{tt})))
alpha(double(~isnan(meanISIcomod.(state).(celltypes{tt})')))
hold on
UnityLine('linecolor','w')
axis xy
xlim([-2.5 1.5])
ylim([-2.5 1.5])
LogScale('xy',10,'nohalf',true)
colorbar
%caxis([-0.15 0.15])
%caxis([-0.002 0.002])
%crameri('berlin','pivot',0)
ylabel('ISI - cell j');xlabel('ISI - cell i')

subplot(2,2,2)
plot(logISIbins,(meanISIcomod.(state).zbin.(celltypes{tt})),'linewidth',2)
hold on
plot(logISIbins,(meanISIcomod.(state).abin.(celltypes{tt})),'linewidth',2)
axis tight
xlim([-2.5 1.5])
xlabel('ISI');ylabel('relative P[ISI]')
legend('No Spikes','Any Spikes','location','northwest')
LogScale('x',10,'nohalf',true)
hold on
plot(xlim(gca),[0 0],'k')

subplot(2,2,3)
imagesc(logISIbins,logISIbins,(meanISIcomod.(state).norm.(celltypes{tt})))
alpha(double(~isnan(meanISIcomod.(state).(celltypes{tt})')))
hold on
UnityLine('linecolor','w')
axis xy
xlim([-2.5 1.5])
ylim([-2.5 1.5])
LogScale('xy',10,'nohalf',true)
colorbar
%caxis([-0.15 0.15])
caxis([-0.0025 0.0025])
crameri('berlin','pivot',0)
ylabel('ISI - cell j');xlabel('ISI - cell i')

subplot(2,2,4)
plot(logISIbins,(meanISIcomod.(state).norm.zbin.(celltypes{tt})),'linewidth',2)
hold on
plot(logISIbins,(meanISIcomod.(state).norm.abin.(celltypes{tt})),'linewidth',2)
axis tight
xlim([-2.5 1.5])
xlabel('ISI');ylabel('relative P[ISI]')
legend('No Spikes','Any Spikes','location','northwest')
LogScale('x',10,'nohalf',true)
hold on
plot(xlim(gca),[0 0],'k')


NiceSave(['ISImod_pE_',state],figfolder,baseName)


end
%%
tt = 1;
excell = 1;
pair = 20;
figure
subplot(3,3,4)
imagesc(logISIbins,logISIbins,(ISIcomod.(state).(celltypes{tt})(:,:,excell)))
alpha(double(~isnan(ISIcomod.(state).(celltypes{tt})(:,:,excell))))
hold on
UnityLine
axis xy
LogScale('xy',10,'nohalf',true)
%colorbar
%caxis([-0.5 0.5])
%crameri('berlin','pivot',0)
ylabel('ISI - cell j (All Pairs)');xlabel('ISI - cell i')
%title(['Cell i: ',num2str(excell)])

subplot(3,3,7)
imagesc(logISIbins,logISIbins,(ISIcomod.(state).norm.(celltypes{tt})(:,:,excell)))
alpha(double(~isnan(ISIcomod.(state).norm.(celltypes{tt})(:,:,excell))))
hold on
UnityLine
axis xy
LogScale('xy',10,'nohalf',true)
%colorbar
caxis([-0.0025 0.0025])
crameri('berlin','pivot',0)
ylabel('ISI - cell j (All Pairs)');xlabel('ISI - cell i')



subplot(6,3,4)
plot(logISIbins,ISIcomod.(state).i_ISIhist_marj(:,excell),'k','linewidth',1)
hold on
plot(logISIbins,(ISIcomod.(state).abin.(celltypes{tt})(:,excell)),'b')
plot(logISIbins,(ISIcomod.(state).zbin.(celltypes{tt})(:,excell)),'r')
%legend('i ISI dist','j Spike','No Spike','location','northwest')
box off
axis tight
xlabel('ISI');ylabel('P[ISI]')
LogScale('x',10,'nohalf',true)
hold on



subplot(3,3,5)
imagesc(logISIbins,logISIbins,(ISIcomod.(state).ISIpars(:,:,pair(1),excell)))
alpha(double(~isnan(ISIcomod.(state).ISIpars(:,:,pair(1),excell))))
hold on
UnityLine
axis xy
LogScale('xy',10,'nohalf',true)
%colorbar
%caxis([-1 1])
%crameri('berlin','pivot',0)
ylabel('ISI - cell j');xlabel('ISI - cell i')


subplot(3,3,8)
imagesc(logISIbins,logISIbins,(ISIcomod.(state).norm.ISIpars_spknorm(:,:,pair(1),excell)))
alpha(double(~isnan(ISIcomod.(state).ISIpars(:,:,pair(1),excell))))
hold on
UnityLine
axis xy
LogScale('xy',10,'nohalf',true)
%colorbar
%caxis([-1 1])
crameri('berlin','pivot',0)
ylabel('ISI - cell j');xlabel('ISI - cell i')



subplot(3,6,11)
plot(ISIStats.ISIhist.(state).log(pair,:),ISIStats.ISIhist.logbins,'k','linewidth',2)
ylim(logISIbins([1 end]))
box off
xlabel('P[ISI_j]')


subplot(3,3,2)
plot(logISIbins,ISIcomod.(state).i_ISIhist_marj(:,excell),'k','linewidth',2)
hold on
plot(logISIbins,(ISIcomod.(state).anyspike(:,pair(1),excell)),'b')
plot(logISIbins,(ISIcomod.(state).nospike(:,pair(1),excell)),'r')
axis tight
legend('i ISI dist','j Spike','No Spike','location','northoutside')
xlabel('ISI');ylabel('P[ISI_i]')
LogScale('x',10,'nohalf',true)
box off
title(['Cell j: ',num2str(pair)])

subplot(6,3,1)
imagesc(logISIbins,[1 spikes.numcells],ISIcomod.(state).anyspike(:,ISIStats.sorts.(state).ratepE,excell)'-...
    ISIcomod.(state).nospike(:,ISIStats.sorts.(state).ratepE,excell)')
crameri('berlin','pivot',0)
title(['Cell i: ',num2str(excell)])


NiceSave('ExCell',figfolder,baseName)

%%





