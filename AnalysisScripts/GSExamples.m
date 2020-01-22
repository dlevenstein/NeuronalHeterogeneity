function [ ] = GSExamples(basePath,figfolder)
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
%figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%Load Stuff
sessionInfo = bz_getSessionInfo(basePath);
spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
states = fieldnames(SleepState.ints);
states{4} = 'ALL';
SleepState.ints.ALL = [0 Inf];
statecolors = {'k','b','r',[0.6 0.6 0.6]};

[celltypes,~,typeidx] = unique(CellClass.label);
cellcolor = {'k','r'};

try
    celltypes = CellClass.celltypes;
catch
    celltypes = unique(CellClass.label);
end
cellcolor = {'k','r'};
%% Load the LFP if needed
% 
% lfpchan = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID;
% downsamplefactor = 1;
% lfp = bz_GetLFP(lfpchan,...
%     'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);
% %Noralize the LFP
% lfp.data = NormToInt(single(lfp.data),'modZ', SleepState.ints.NREMstate,lfp.samplingRate);

%% Restrict to state

% state = states{3};
% ISIStats.allspikes.instate = cellfun(@(X) InIntervals(X,double(SleepState.ints.(state))),...
%     ISIStats.allspikes.times,'UniformOutput',false);

%%
ISIrate = bz_LoadCellinfo(basePath,'ISIrate');

%%
ISIoccupancy.logbins = linspace(-3,3,100);
    for ss = 1:3
        ISIrate.OccupancyStats.(states{ss}).MTORatio = ...
            ISIStats.summstats.(states{ss}).meanrate./...
            (1./ISIrate.OccupancyStats.(states{ss}).median);
        
        ISIrate.instate = ...
            InIntervals(ISIrate.timestamps,SleepState.ints.(states{ss}));   

        ISIoccupancy.(states{ss}).loghist = ...
            hist(log10(ISIrate.ISI(ISIrate.instate,:)),ISIoccupancy.logbins);
        ISIoccupancy.(states{ss}).loghist = ...
            ISIoccupancy.(states{ss}).loghist./length(ISIrate.timestamps(ISIrate.instate));
        
    end
%%
histcolors = flipud(gray);
NREMhistcolors = makeColorMap([1 1 1],[0 0 0.8]);
REMhistcolors = makeColorMap([1 1 1],[0.8 0 0]);
statecolormap = {histcolors,NREMhistcolors,REMhistcolors};

for cc = 1:length(celltypes)
excell = randsample(find(CellClass.(celltypes{cc})),4,true);

for ss = 1:3

% excell(1) = randsample(find(log10(ISIrate.OccupancyStats.(states{ss}).MTORatio)<0.4 & CellClass.pE),1);
% excell(2) = randsample(find(log10(ISIrate.OccupancyStats.(states{ss}).MTORatio)>0.8 & CellClass.pE),1);
% 


[~,sortex] = sort(ISIrate.OccupancyStats.(states{ss}).median(excell));
excell = excell(sortex);
exwin = bz_RandomWindowInIntervals(SleepState.ints.(states{ss}),60);

figure
for ee = 1:4
    subplot(6,3,(ee-1).*3+[1:2])
    
        plotspikes = Restrict(spikes.times{excell(ee)},exwin);
        plot(ISIrate.timestamps,log10(ISIrate.ISI(:,excell(ee))),'k')
        hold on
        plot([plotspikes*[1 1]]',[ones(size(plotspikes))*[1.75 2]]','k')
        plot(exwin,log10(1./ISIStats.summstats.(states{ss}).meanrate(excell(ee))).*[1 1],'k:')
        plot(exwin,log10(ISIrate.OccupancyStats.(states{ss}).median(excell(ee))).*[1 1],'r:')
        xlim(exwin)
        ylim([-3 2.3])
        LogScale('y',10,'exp',true)
        ylabel({['UID: ',num2str(spikes.UID(excell(ee)))],'ISI (s)'})
        box off
        
        bz_ScaleBar('s')
        axis ij
        
%      subplot(6,4,(ee)+16)
%         plot(ISIStats.ISIhist.logbins,ISIStats.ISIhist.(states{ss}).log(excell(ee),:),'k','linewidth',2)
%         hold on
%         plot(ISIoccupancy.logbins,ISIoccupancy.(states{ss}).loghist(:,excell(ee)),'k:')
%         plot(log10(1./ISIStats.summstats.(states{ss}).meanrate(excell(ee))),0,'k+')
%         plot(log10(ISIrate.OccupancyStats.(states{ss}).median(excell(ee))),0,'r+')
%         box off
%         axis tight
%         LogScale('x',10,'exp',true)
%         set(gca,'yticklabels',[])
%         xlabel('ISI (s)');ylabel('P[ISI]')
%         %plot(ISIrate.OccupancyStats.WAKEstate.
        
     subplot(6,4,ee+16)
        imagesc(ISIStats.ISIhist.logbins,ISIStats.CV2hist.bins,...
            squeeze(ISIStats.Jointhist.(states{ss}).log(excell(ee),:,:))')
        axis xy
        hold on
        plot(ISIoccupancy.logbins,bz_NormToRange(ISIoccupancy.(states{ss}).loghist(:,excell(ee)),0.75),':','color',statecolors{ss})
        plot(ISIStats.ISIhist.logbins,bz_NormToRange(ISIStats.ISIhist.(states{ss}).log(excell(ee),:),0.75),...
            'linewidth',2,'color',statecolors{ss})
        
        plot(log10(1./ISIStats.summstats.(states{ss}).meanrate(excell(ee))),0,'k+')
        plot(log10(ISIrate.OccupancyStats.(states{ss}).median(excell(ee))),0,'r+')
        %ISIStats.Jointhist.(states{ss}).log(excell(ee),1,:)
        caxis([0 max([ISIStats.Jointhist.(states{ss}).log(excell(ee),:,1),0])])
        LogScale('x',10,'exp',true)
        xlim([-3 2])
        xlabel('ISI (s)');ylabel('CV2')
    
    subplot(6,4,ee+20)   
    imagesc((ISIStats.ISIhist.logbins),(ISIStats.ISIhist.logbins),(ISIStats.ISIhist.(states{ss}).return(:,:,excell(ee))))
    hold on
    plot(log10(1./ISIStats.summstats.(states{ss}).meanrate(excell(ee))),log10(1./ISIStats.summstats.(states{ss}).meanrate(excell(ee))),'k+')
    axis xy
    LogScale('xy',10,'exp',true)
    %set(gca,'ytick',[]);set(gca,'xtick',[]);
    caxis([0 0.003])
    xlim(ISIStats.ISIhist.logbins([1 end]));ylim(ISIStats.ISIhist.logbins([1 end]))
    xlim([-3 2]);ylim([-3 2])
    colormap(gca,statecolormap{ss})
end

subplot(3,3,3)
plot(log10((1./ISIrate.OccupancyStats.(states{ss}).median(CellClass.(celltypes{cc})))),...
    log10(ISIrate.OccupancyStats.(states{ss}).MTORatio(CellClass.(celltypes{cc}))),'k.')
hold on
plot(log10((1./ISIrate.OccupancyStats.(states{ss}).median(excell))),...
    log10(ISIrate.OccupancyStats.(states{ss}).MTORatio(excell)),'r.')

NiceSave([(celltypes{cc}),'_GSExamples_',(states{ss})],figfolder,baseName)
end
end
