function [ output_args ] = HPCSpikeIrregularityAnalysis( basePath,figfolder )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

basePath = '/mnt/packrat/userdirs/david/zpool2/DT2/DT2_rPPC_rCCG_218um_218um_20160208_160208_142910'; %Parietal/Cingulate
basePath = '/mnt/packrat/userdirs/david/zpool2/DT2/DT2_rPPC_rCCG_3540um_1288um_20160227_160227_121226';


figfolder = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/AnalysisScripts/AnalysisFigs/HPCSpikeIrregularityAnalysis';
%%
sessionInfo = bz_getSessionInfo(basePath);
baseName = bz_BasenameFromBasepath(basePath);

%%
spikes = bz_GetSpikes('basepath',basePath);
%%
%Needs work.....
SleepState = bz_LoadStates(basePath,'SleepState');
%%
CellClass = bz_LoadCellinfo(basePath,'CellClass');
%CellClass = CellClass.CellClass;
%%
CellClass.label = cell(size(CellClass.Pyr));
CellClass.label(CellClass.Pyr) = {'Pyr'};
CellClass.label(CellClass.Int) = {'Int'};

%%
hasfield = cellfun(@(X) X.has_field,Tuning.placeFields);
CellClass.PyrField = CellClass.Pyr & hasfield;
CellClass.PyrNo = CellClass.Pyr & ~hasfield;

%% Calculate pop CV2, rate over time


dt = 0.25;
binsize = 2;
overlap = binsize./dt;
spikemat = bz_SpktToSpkmat(spikes,'binsize',binsize,'overlap',overlap);
%%
cellclasses = {'Pyr','Int','PyrField','PyrNo'};
classcolors = {'k','r','b','g'};
for pp = 1:length(cellclasses)
    spikemat.poprate.(cellclasses{pp}) = ...
        sum(spikemat.data(:,CellClass.(cellclasses{pp})),2)... %Spikes/s/cell
        ./binsize...
        ./sum(CellClass.(cellclasses{pp}));
end

%% Calculate the spike stats
[ ISIStats ] = bz_ISIStats( spikes,'cellclass',CellClass.label,'showfig',true,...
    'ints',SleepState.ints);

%% Mean binned CV2...
clear CV2mat
CV2mat.winsize = binsize;
CV2mat.timestamps = spikemat.timestamps;
CV2mat.binedges = bsxfun(@(X,Y) X+Y,spikemat.timestamps,[-0.5 0.5].*CV2mat.winsize);
for pp = 1:length(cellclasses)
    allspikes.CV2.(cellclasses{pp}) = cat(1,ISIStats.allspikes.CV2{CellClass.(cellclasses{pp})});
    allspikes.times.(cellclasses{pp}) = cat(1,ISIStats.allspikes.times{CellClass.(cellclasses{pp})});
    [CV2mat.timestamps,CV2mat.(cellclasses{pp})] = ...
        BinDataTimes(allspikes.CV2.(cellclasses{pp}),allspikes.times.(cellclasses{pp}),CV2mat.binedges);
    CV2mat.rate.(cellclasses{pp}) = interp1(spikemat.timestamps,spikemat.poprate.(cellclasses{pp}),CV2mat.timestamps);
end


%% Get the Behavior
linear = bz_LoadBehavior(basePath,'linear');

%% 
lfp = bz_GetLFP(sessionInfo.thetaChans(2),'basepath',basePath);

%% PSS
% dt = 0.2;
% binsize = 2;
[specslope,specgram] = bz_PowerSpectrumSlope(lfp,binsize,dt,'showfig',true);

%% mean over trials
aroundtrialwindow = 4;
linear.numtrials = length(linear.events.trialIntervals);
buffer = 2;

numbins = 90;
trialdata.linearizedposition = linspace(-1,2,numbins);
for pp = 1:length(cellclasses)
    trialdata.popCV2.(cellclasses{pp}) = nan(numbins,linear.numtrials);
    trialdata.poprate.(cellclasses{pp}) = nan(numbins,linear.numtrials);
end
trialdata.PSS = nan(numbins,linear.numtrials);
trialdata.specgram = nan(numbins,length(specgram.freqs),linear.numtrials);

[~,sorttrialtimes] = sort(linear.events.trialIntervals(:,1));
sortedtrialtimes = linear.events.trialIntervals(sorttrialtimes,:);

for bb = 2:linear.numtrials-1

    pretrialwindow= sortedtrialtimes(bb,1) + aroundtrialwindow.*[-1 0];
    trialwindow= sortedtrialtimes(bb,:);
    posttrialwindow= sortedtrialtimes(bb,2) + aroundtrialwindow.*[0 1];
    betweenothertrials = [sortedtrialtimes(bb-1,2)+buffer sortedtrialtimes(bb+1,1)-buffer];
    
    %Spike Rate Stuff
    pretrialtimes = InIntervals(spikemat.timestamps,pretrialwindow);
    trialtimes = InIntervals(spikemat.timestamps,trialwindow);
    posttrialtimes = InIntervals(spikemat.timestamps,posttrialwindow);
    linearizedposition = [linspace(-1,-0.01,sum(pretrialtimes)),...
        linspace(0.01,0.99,sum(trialtimes)),linspace(1.01,2,sum(posttrialtimes))]';
    
    %othertrialtimes
    othertrialtimes = ~InIntervals(spikemat.timestamps,betweenothertrials);
    pretrialnans = sum(pretrialtimes&othertrialtimes);
    postrialnans = sum(posttrialtimes&othertrialtimes);
    
%     linearizedposition(1:pretrialnans) = nan;
%     linearizedposition(end-postrialnans:end) = nan;

    for pp = 1:length(cellclasses)
        
        CV2trial = [CV2mat.(cellclasses{pp})(pretrialtimes);...
            CV2mat.(cellclasses{pp})(trialtimes);...
            CV2mat.(cellclasses{pp})(posttrialtimes)];
        CV2trial(1:pretrialnans)=nan;
        if postrialnans>1
            CV2trial(end-postrialnans:end)=nan;
        end
        
        popratetrial = [spikemat.poprate.(cellclasses{pp})(pretrialtimes);...
            spikemat.poprate.(cellclasses{pp})(trialtimes);...
            spikemat.poprate.(cellclasses{pp})(posttrialtimes)];
        popratetrial(1:pretrialnans)=nan;
        if postrialnans>1
            popratetrial(end-postrialnans:end)=nan;
        end
        
        trialdata.popCV2.(cellclasses{pp})(:,bb) = ...
            interp1(linearizedposition,CV2trial,trialdata.linearizedposition)';
        
        trialdata.poprate.(cellclasses{pp})(:,bb) = ...
            interp1(linearizedposition,popratetrial,trialdata.linearizedposition)';
    end
 
    %LFP Spectrum Stuff
    pretrialtimes = InIntervals(specslope.timestamps,pretrialwindow);
    trialtimes = InIntervals(specslope.timestamps,trialwindow);
    posttrialtimes = InIntervals(specslope.timestamps,posttrialwindow);
    linearizedposition = [linspace(-1,-0.01,sum(pretrialtimes)),...
        linspace(0,1,sum(trialtimes)),linspace(1.01,2,sum(posttrialtimes))]';
    
        %othertrialtimes
    othertrialtimes = ~InIntervals(specslope.timestamps,betweenothertrials);
    pretrialnans = sum(pretrialtimes&othertrialtimes);
    postrialnans = sum(posttrialtimes&othertrialtimes);
    
%     linearizedposition(1:pretrialnans) = nan;
%     linearizedposition(end-postrialnans:end) = nan;
    
    specslopetrial = [specslope.data(pretrialtimes);...
        specslope.data(trialtimes);...
        specslope.data(posttrialtimes)];
    specslopetrial(1:pretrialnans)=nan;
    if postrialnans>1
        specslopetrial(end-postrialnans:end)=nan;
    end
    
    specgramtrial = [specgram.amp(:,pretrialtimes)';...
        specgram.amp(:,trialtimes)';...
        specgram.amp(:,posttrialtimes)'];
    specgramtrial(1:pretrialnans,:)=nan;
    if postrialnans>1
        specgramtrial(end+1-postrialnans:end,:)=nan;
    end
    
    trialdata.PSS(:,bb) = ...
        interp1(linearizedposition,specslopetrial,trialdata.linearizedposition)';
    
    trialdata.specgram(:,:,bb) = ...
        interp1(linearizedposition,specgramtrial,trialdata.linearizedposition);
    
    
    %Individual spike stuff
    normtime = cellfun(@(X) IntervalTimeNormalize(X,trialwindow,aroundtrialwindow),ISIStats.allspikes.times,'UniformOutput',false);
    thistrial = cellfun(@(X) X>betweenothertrials(1) & X<betweenothertrials(2),ISIStats.allspikes.times,'UniformOutput',false);
    %Is it within the window and in the domain of THIS SW
    trialrelidx = cellfun(@(X,Y) X>-1 & X<2 & Y,normtime,thistrial,'UniformOutput',false);
    
    ISIn(bb,:) = cellfun(@(X,Y) X(Y),ISIStats.allspikes.ISIs,trialrelidx,'UniformOutput',false);
    ISInp1(bb,:) = cellfun(@(X,Y) X([false; Y(1:end-1)]),ISIStats.allspikes.ISIs,trialrelidx,'UniformOutput',false);
    trialrelCV2(bb,:) = cellfun(@(X,Y) X(Y),ISIStats.allspikes.CV2,trialrelidx,'UniformOutput',false);
    trialnormtime(bb,:) = cellfun(@(X,Y) X(Y),normtime,trialrelidx,'UniformOutput',false);
    
    
    
end

for pp = 1:length(cellclasses)
    grouptrialdata.popCV2.mean.(cellclasses{pp}) = nanmean(trialdata.popCV2.(cellclasses{pp}),2);
    grouptrialdata.popCV2.std.(cellclasses{pp}) = nanstd(trialdata.popCV2.(cellclasses{pp}),[],2);
    grouptrialdata.poprate.mean.(cellclasses{pp}) = nanmean(log2(trialdata.poprate.(cellclasses{pp})),2);
    grouptrialdata.poprate.std.(cellclasses{pp}) = nanstd(log2(trialdata.poprate.(cellclasses{pp})),[],2);
end
grouptrialdata.PSS.mean = nanmean(trialdata.PSS,2);
grouptrialdata.PSS.std = nanstd(trialdata.PSS,[],2);
grouptrialdata.specgram = nanmean(trialdata.specgram,3);

%%
for cc = 1:spikes.numcells
    trialrelCV2_all{cc} = cat(1,trialrelCV2{:,cc});
    trialnormtime_all{cc} = cat(1,trialnormtime{:,cc});
    ISIn_all{cc} = cat(1,ISIn{:,cc});
    ISInp1_all{cc} = cat(1,ISInp1{:,cc});
    
    if length(ISIn_all{cc})<100
        ISIn_all{cc} = [nan;nan;nan];
        ISInp1_all{cc} = [nan;nan;nan];
        trialrelCV2_all{cc} = [nan;nan;nan];
        trialnormtime_all{cc} = [-0.5 0.5 1.5];
    end
end

%%
binsize = 0.1;
binedges = -1:binsize:2;
[ bincenters,binmeans,binstd,binnum,binneddata ] = cellfun(@(X,Y) BinDataTimes(X,Y,binedges),trialrelCV2_all,trialnormtime_all,'UniformOutput',false);
bincenters = bincenters{1};
binmeans = cat(2,binmeans{:});
binstd = cat(2,binstd{:});

for tt = 1:length(cellclasses)
trialCV2bypop.(cellclasses{tt}).mean = nanmean(binmeans(:,CellClass.(cellclasses{tt})),2);
trialCV2bypop.(cellclasses{tt}).std = nanstd(binmeans(:,CellClass.(cellclasses{tt})),[],2);
end

%%
figure
subplot(4,2,5)
plot([0 0],[0.7 1.5],'k')
hold on
for tt = 1:length(cellclasses)
    plot(bincenters,trialCV2bypop.(cellclasses{tt}).mean,'o-','color',classcolors{tt})
    hold on
    errorshade(bincenters,trialCV2bypop.(cellclasses{tt}).mean,...
        trialCV2bypop.(cellclasses{tt}).std,trialCV2bypop.(cellclasses{tt}).std,classcolors{tt},'scalar')
end
ylim([0.8 1.4])
xlabel('t - relative to Trial');ylabel('<CV2>')

%%
%% ISI return maps around trial
%(add ACG)

win = 0.15;

outtrialreturnmaps = cellfun(@(X,Y,Z) hist3([log10(X(Z<0|Z>1)) log10(Y(Z<0|Z>1))],{ISIStats.ISIhist.logbins,ISIStats.ISIhist.logbins}),...
    ISIn_all,ISInp1_all,trialnormtime_all,'UniformOutput',false);
outtrialreturnmaps = cellfun(@(X) X./sum(X(:)),outtrialreturnmaps,'UniformOutput',false);
outtrialreturnmaps = cat(3,outtrialreturnmaps{:});

intrialreturnmaps = cellfun(@(X,Y,Z) hist3([log10(X(Z<1&Z>0)) log10(Y(Z<1&Z>0))],{ISIStats.ISIhist.logbins,ISIStats.ISIhist.logbins}),...
    ISIn_all,ISInp1_all,trialnormtime_all,'UniformOutput',false);
intrialreturnmaps = cellfun(@(X) X./sum(X(:)),intrialreturnmaps,'UniformOutput',false);
intrialreturnmaps = cat(3,intrialreturnmaps{:});


for tt = 1:length(cellclasses)
    meanouttrialreturn.(cellclasses{tt}) = nanmean(outtrialreturnmaps(:,:,CellClass.(cellclasses{tt})),3);
    meanintrialreturn.(cellclasses{tt}) = nanmean(intrialreturnmaps(:,:,CellClass.(cellclasses{tt})),3);
end

%%
histcolors = flipud(gray);
figure
for tt = 1:length(cellclasses)
    subplot(3,4,0+tt)
    colormap(gca,histcolors)
        imagesc(ISIStats.ISIhist.logbins,ISIStats.ISIhist.logbins,...
            meanouttrialreturn.(cellclasses{tt})')
        LogScale('xy',10)
        %xlabel('ISI_n');ylabel('ISI_n+1')
        axis xy
        caxis([0 5e-3])
        title(cellclasses{tt})
        %colorbar
        
        
    subplot(3,4,8+tt)
    colormap(gca,histcolors)
        imagesc(ISIStats.ISIhist.logbins,ISIStats.ISIhist.logbins,...
            meanintrialreturn.(cellclasses{tt})')
        LogScale('xy',10)
        %xlabel('Preceeding ISI');ylabel('Next ISI')
        axis xy
        caxis([0 5e-3])
        
end

NiceSave('trialReturnmaps',figfolder,baseName)
%%
figure
subplot(2,1,1)
    imagesc(trialdata.linearizedposition,[0 linear.numtrials],trialdata.popCV2.Pyr')
subplot(2,1,2)
    imagesc(trialdata.linearizedposition,[0 linear.numtrials],trialdata.popCV2.Int')
%%
randtrial = randsample(length(linear.events.trialIntervals),1);
%randtrial = bb;
timewin = sortedtrialtimes(randtrial,:) + 10.*[-1 1];

figure

subplot(6,1,1)
    imagesc(specgram.timestamps,log2(specgram.freqs),specgram.amp)
    hold on
   % StateScorePlot({SleepState.ints.NREMstate,SleepState.ints.REMstate,SleepState.ints.WAKEstate},...
   %     {'b','r','k'})
    axis xy
    xlim(timewin)
    ylabel({'Specgram','f (Hz)'})
    LogScale('y',2)
    


subplot(6,1,4)
    plot(specslope.timestamps,specslope.data,'k','linewidth',1)
    axis tight
    xlim(timewin)
    ylabel('PSS')
    box off
    

subplot(6,1,[2:3])
    bz_MultiLFPPlot(lfp,'spikes',spikes,'timewin',timewin,...
      	 'sortmetric',ISIStats.summstats.WAKEstate.meanrate )
     %        'cellgroups',{CellClass.Pyr},...
       % ,...
        %)
    hold on
    plot(linear.events.trialIntervals',ones(size(linear.events.trialIntervals))','r')
    xlim(timewin)
    
subplot(6,1,5)
for pp = 1:length(cellclasses)
    plot(spikemat.timestamps,log2(spikemat.poprate.(cellclasses{pp})),classcolors{pp},'linewidth',1)
    hold on
end
    axis tight
    xlim(timewin)
    box off
    ylabel({'Pop. Rate', '(Spk/Cell/S)'})
    legend(cellclasses{:})
    LogScale('y',2)
    
subplot(6,1,6)
for pp = 1:length(cellclasses)
    plot(CV2mat.timestamps,CV2mat.(cellclasses{pp}),classcolors{pp},'linewidth',1)
    hold on
end
    axis tight
    xlim(timewin)
    box off
    ylabel('<CV2>')
   
NiceSave('exampletrial',figfolder,baseName)  
    %%
figure
subplot(4,2,6)
for pp = 1:length(cellclasses)
    plot(trialdata.linearizedposition,grouptrialdata.popCV2.mean.(cellclasses{pp}),classcolors{pp})
    hold on
    errorshade(trialdata.linearizedposition,grouptrialdata.popCV2.mean.(cellclasses{pp}),...
        grouptrialdata.popCV2.std.(cellclasses{pp}),grouptrialdata.popCV2.std.(cellclasses{pp}),...
        classcolors{pp},'scalar')
end
axis tight
box off
plot([0 0],get(gca,'ylim'),'k')
plot([1 1],get(gca,'ylim'),'k')
ylabel('<CV2>')
set(gca,'XTick',[-1 0 1 2]);
set(gca,'XTickLabels',{['-',num2str(aroundtrialwindow)],'S','E',['+',num2str(aroundtrialwindow)]})

subplot(4,2,8)
for pp = 1:length(cellclasses)
    plot(trialdata.linearizedposition,grouptrialdata.poprate.mean.(cellclasses{pp}),classcolors{pp})
    hold on
    errorshade(trialdata.linearizedposition,grouptrialdata.poprate.mean.(cellclasses{pp}),...
        grouptrialdata.poprate.std.(cellclasses{pp}),grouptrialdata.poprate.std.(cellclasses{pp}),...
        classcolors{pp},'scalar')
end
LogScale('y',2)
axis tight
box off
plot([0 0],get(gca,'ylim'),'k')
plot([1 1],get(gca,'ylim'),'k')
ylabel('Pop Rate')
xlabel('Trial')
set(gca,'XTick',[-1 0 1 2]);
set(gca,'XTickLabels',{['-',num2str(aroundtrialwindow),'s'],'S','E',['+',num2str(aroundtrialwindow),'s']})

subplot(4,2,4)
    plot(trialdata.linearizedposition,grouptrialdata.PSS.mean,'k')
    hold on
    errorshade(trialdata.linearizedposition,grouptrialdata.PSS.mean,...
        grouptrialdata.PSS.std,grouptrialdata.PSS.std,...
        'k','scalar')
    axis tight
    ylabel('PSS')
    plot([0 0],get(gca,'ylim'),'k')
plot([1 1],get(gca,'ylim'),'k')
    set(gca,'XTick',[-1 0 1 2]);
set(gca,'XTickLabels',{['-',num2str(aroundtrialwindow)],'S','E',['+',num2str(aroundtrialwindow)]})
    
subplot(4,2,2)
    imagesc(trialdata.linearizedposition,log2(specgram.freqs),grouptrialdata.specgram')
    LogScale('y',2)
    axis xy    
    hold on
    plot([0 0],get(gca,'ylim'),'k')
plot([1 1],get(gca,'ylim'),'k')
    set(gca,'XTick',[-1 0 1 2]);
set(gca,'XTickLabels',{['-',num2str(aroundtrialwindow)],'S','E',['+',num2str(aroundtrialwindow)]})
    
NiceSave('meantrial',figfolder,baseName)
end

