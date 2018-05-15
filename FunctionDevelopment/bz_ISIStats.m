function [ ISIstats ] = bz_ISIStats( spikes,varargin )
%ISIstats = bz_ISIStats(spikes,varargin) calculates the statistics 
%inter-spike intervals for the spiketimes in spikes.
%
%   INPUTS
%       spikes
%
%       (options)
%       'states'        A structure with intervals.
%                       states.statename = [start stop]
%                       Will calculate ISIsstats separately for each state
%       'cellclass'     Cell array of strings - label for each cell
%       'savecellinfo'
%       'basePath'
%       'figfolder'
%       'showfig'
%
%   OUTPUTS
%       ISIstats        cellinfo structure with ISI statistics
%           .summstats  summary statistics
%           .ISIhist    histograms of ISIs etc
%           .sorts      sorting indices
%           .allspikes  ISI/CV2 value for each spike (ISI is PRECEDING interval)
%
%DLevenstein 2018
%% Parse the inputs
defaultstates.ALLtime = [-Inf Inf];

% parse args
p = inputParser;
addParameter(p,'states',defaultstates)
addParameter(p,'savecellinfo',false,@islogical)
addParameter(p,'basePath',pwd,@isstr)
addParameter(p,'figfolder',false)
addParameter(p,'showfig',false,@islogical);
addParameter(p,'cellclass',[]);


parse(p,varargin{:})
states = p.Results.states;
cellclass = p.Results.cellclass;
basePath = p.Results.basePath;
SAVECELLINFO = p.Results.savecellinfo;
SAVEFIG = p.Results.figfolder;
SHOWFIG = p.Results.showfig;

%%
statenames = fieldnames(states);
numstates = length(statenames);


%% ISI and CV2 statistics
numcells = length(spikes.UID);

%Calculate ISI and CV2 for allspikes
allspikes.ISIs = cellfun(@diff,spikes.times,'UniformOutput',false);
allspikes.CV2 = cellfun(@(X) 2.*abs(X(2:end)-X(1:end-1))./(X(2:end)+X(1:end-1)),allspikes.ISIs ,'UniformOutput',false);
%Make sure times line up
allspikes.times = cellfun(@(X) X(2:end-1),spikes.times,'UniformOutput',false);
allspikes.ISIs = cellfun(@(X) X(1:end-1),allspikes.ISIs,'UniformOutput',false);
%%
for ss = 1:numstates
%ss=1;

%Find which spikes are during state of interest
[statespiketimes,statespikes] = cellfun(@(X) RestrictInts(X,states.(statenames{ss})),...
    allspikes.times,'UniformOutput',false);
CV2 = cellfun(@(X,Y) X(Y(2:end-1)),allspikes.CV2,statespikes,'Uniformoutput',false);
ISIs = cellfun(@(X,Y) X(Y(2:end-1)),allspikes.ISIs,statespikes,'Uniformoutput',false);

%Summary Statistics
summstats.(statenames{ss}).meanISI = cellfun(@(X) mean(X),ISIs);
summstats.(statenames{ss}).meanrate = 1./summstats.(statenames{ss}).meanISI;
summstats.(statenames{ss}).ISICV = cellfun(@(X) std(X)./mean(X),ISIs);
summstats.(statenames{ss}).meanCV2 = cellfun(@(X) mean(X),CV2);


%% CV2 for shuffle (shows that CV2 is not much meaningful?)
% numshuffle = 100;
% for sh = 1:numshuffle
%     ISIs_shuffle = cellfun(@(X) shuffle(X),ISIs,'UniformOutput',false);
%     CV2_shuffle = cellfun(@(X) 2.*abs(X(2:end)-X(1:end-1))./(X(2:end)+X(1:end-1)),...
%         ISIs_shuffle ,'UniformOutput',false);
%     meanshuffle(sh,:) = cellfun(@(X) mean(X),CV2_shuffle);
% end
% shufflemean = mean(meanshuffle);
% shufflestd = std(meanshuffle);
% 
% %CV2_reltoshuff = (summstats.(statenames{ss}).meanCV2-shufflemean)./shufflestd;
% 
% figure
% subplot(2,2,1)
% plot([log10(summstats.(statenames{ss}).meanrate);log10(summstats.(statenames{ss}).meanrate)],...
%     [shufflemean;summstats.(statenames{ss}).meanCV2],'color',0.7.*[1 1 1],'linewidth',0.5)
% hold on
% plot([log10(summstats.(statenames{ss}).meanrate);log10(summstats.(statenames{ss}).meanrate)],...
%     [shufflemean-shufflestd;shufflemean+shufflestd],'color',0.7.*[1 1 1],'linewidth',3)
% plot(log10(summstats.(statenames{ss}).meanrate),summstats.(statenames{ss}).meanCV2,'.r','markersize',10)
% xlabel('FR (Hz)');ylabel('<CV2>')
% NiceSave(['CV2_',(statenames{ss})],figfolder,baseName);

%%
%Set up all the bins and matrices
numbins = 60;
ISIhist.linbins = linspace(0,10,numbins);
ISIhist.logbins = linspace(log10(0.001),log10(200),numbins);
ISIhist.(statenames{ss}).lin = zeros(numcells,numbins);
ISIhist.(statenames{ss}).log = zeros(numcells,numbins);
normcv2hist = zeros(numcells,numbins);

ISIhist.(statenames{ss}).return = zeros(numbins,numbins,numcells);

%Calculate all the histograms: ISI, log(ISI), 1/ISI, log(1/ISI)
for cc = 1:numcells
    numspks(cc) = length(ISIs{cc});
    
    %Calculate ISI histograms
    ISIhist.(statenames{ss}).lin(cc,:) = hist(ISIs{cc},ISIhist.linbins);
    ISIhist.(statenames{ss}).log(cc,:) = hist(log10(ISIs{cc}),ISIhist.logbins);
    
    %Normalize histograms to number of spikes
    ISIhist.(statenames{ss}).lin(cc,:) = ISIhist.(statenames{ss}).lin(cc,:)./numspks(cc);
    ISIhist.(statenames{ss}).log(cc,:) = ISIhist.(statenames{ss}).log(cc,:)./numspks(cc);
    
    %Calculate Return maps
    if numspks(cc)>1
    ISIhist.(statenames{ss}).return(:,:,cc) = hist3(log10([ISIs{cc}(1:end-1) ISIs{cc}(2:end)]),{ISIhist.logbins,ISIhist.logbins});
    end
    ISIhist.(statenames{ss}).return(:,:,cc) = ISIhist.(statenames{ss}).return(:,:,cc)./numspks(cc);
  
end

%Sortings
[~,sorts.(statenames{ss}).rate]=sort(summstats.(statenames{ss}).meanrate);
[~,sorts.(statenames{ss}).ISICV]=sort(summstats.(statenames{ss}).ISICV);
[~,sorts.(statenames{ss}).CV2]=sort(summstats.(statenames{ss}).meanCV2);

%Make the cell-type specific sortings
if ~isempty(cellclass)
    classnames = unique(cellclass);
    for cl = 1:length(classnames)
        inclasscells = strcmp(classnames{cl},cellclass);
        sorttypes = {'rate','ISICV','CV2'};
        for tt = 1:length(sorttypes)
        sorts.(statenames{ss}).([sorttypes{tt},classnames{cl}]) = ...
            intersect(sorts.(statenames{ss}).(sorttypes{tt}),find(inclasscells),'stable');
        
        if cl==1
            sorts.(statenames{ss}).([sorttypes{tt},'byclass'])=[];
        end
        sorts.(statenames{ss}).([sorttypes{tt},'byclass']) = ...
            [sorts.(statenames{ss}).([sorttypes{tt},'byclass']) sorts.(statenames{ss}).([sorttypes{tt},classnames{cl}])];
        end
            
    end  
end



%% ACG
%[ccg,t] = CCG(statespiketimes,[],<options>)

%%
if SHOWFIG || SAVEFIG
figure
    subplot(2,2,1)
        plot(log10(summstats.(statenames{ss}).meanrate(CellClass.pE)),...
            log2(summstats.(statenames{ss}).ISICV(CellClass.pE)),'k.')
        hold on
        plot(log10(summstats.(statenames{ss}).meanrate(CellClass.pI)),...
            log2(summstats.(statenames{ss}).ISICV(CellClass.pI)),'r.')
        LogScale('x',10);LogScale('y',2);
        xlabel('Mean Rate (Hz)');ylabel('ISI CV')
        title(statenames{ss})
        box off
        
    subplot(2,2,3)
        plot(log10(summstats.(statenames{ss}).meanrate(CellClass.pE)),...
            (summstats.(statenames{ss}).meanCV2(CellClass.pE)),'k.')
        hold on
        plot(log10(summstats.(statenames{ss}).meanrate(CellClass.pI)),...
            (summstats.(statenames{ss}).meanCV2(CellClass.pI)),'r.')
        plot(get(gca,'xlim'),[1 1])
        LogScale('x',10);
        xlabel('Mean Rate (Hz)');ylabel('ISI <CV2>')
        title(statenames{ss})
        box off
        
    subplot(2,2,2)
        imagesc((ISIhist.logbins),[1 numcells],...
            ISIhist.(statenames{ss}).log(sorts.(statenames{ss}).rateEI,:))
        hold on
        plot(log10(1./(summstats.(statenames{ss}).meanrate(sorts.(statenames{ss}).rateEI))),[1:numcells],'k.','LineWidth',2)
        plot(ISIhist.logbins([1 end]),sum(CellClass.pE).*[1 1]+0.5,'r')
        LogScale('x',10)
        xlabel('ISI (s)')
        xlim(ISIhist.logbins([1 end]))
        colorbar
      %  legend('1/Mean Firing Rate (s)','location','southeast')
        ylabel('Cell (Sorted by FR, Type)')
        %legend('1/Mean Firing Rate (s)','location','southeast')
        caxis([0 0.1])
        title('ISI Distribution (Log Scale)')
        
    subplot(2,2,4)
        imagesc((ISIhist.logbins),[1 numcells],...
            ISIhist.(statenames{ss}).log(sorts.(statenames{ss}).ISICVEI,:))
        hold on
       % plot(log10(1./(summstats.(statenames{ss}).meanrate(sorts.(statenames{ss}).rateEI))),[1:numcells],'k.','LineWidth',2)
        plot(ISIhist.logbins([1 end]),sum(CellClass.pE).*[1 1]+0.5,'r')
        LogScale('x',10)
        xlabel('ISI (s)')
        xlim(ISIhist.logbins([1 end]))
        colorbar
      %  legend('1/Mean Firing Rate (s)','location','southeast')
        ylabel('Cell (Sorted by CV2, Type)')
        %legend('1/Mean Firing Rate (s)','location','southeast')
        caxis([0 0.1])
        title('ISI Distribution (Log Scale)')
        
%     subplot(2,2,3)
%         plot(log2(summstats.(statenames{ss}).meanrate(CellClass.pE)),...
%             log2(summstats.(statenames{ss}).meanCV2(CellClass.pE)),'k.')
%         hold on
%         plot(log2(summstats.(statenames{ss}).meanrate(CellClass.pI)),...
%             log2(summstats.(statenames{ss}).meanCV2(CellClass.pI)),'r.')
%         LogScale('xy',2)
%         xlabel('Mean Rate (Hz)');ylabel('Mean CV2')
%         title(statenames{ss})
%         box off


if SAVEFIG
    NiceSave(['ISIstats_',(statenames{ss})],figfolder,baseName);
end


%exneurons %top/bottom 25 %ile rate/ISICV
%%
% exwindur = 4; %s
% STATEtimepoints = Restrict(lfp.timestamps,double(SleepState.ints.(statenames{ss})));
% samplewin = STATEtimepoints(randi(length(STATEtimepoints))) + [0 exwindur];
% %%
% figure
% bz_MultiLFPPlot( lfp,'spikes',spikes,'timewin',samplewin,...
%     'sortmetric',summstats.(statenames{ss}).meanCV2,...
%     'cellgroups',{CellClass.pI,CellClass.pE})


end

ISIstats.summstats = summstats;
ISIstats.ISIhist = ISIhist;
ISIstats.sorts = sorts;
ISIstats.UID = spikes.UID;
ISIstats.allspikes = allspikes;

if SAVECELLINFO
    save
end
    
end

