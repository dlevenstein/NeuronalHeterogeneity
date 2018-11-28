function [ EIMUARateMap,EIBalRateMap,poprate,spikemat ] = bz_EIMUARateMap( spikes,CellClass,varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
% requires: spikes.times
%%
% parse args
p = inputParser;
addParameter(p,'intervals',[0 Inf],@isnumeric)
addParameter(p,'binsize',0.2)
addParameter(p,'overlap',10)
addParameter(p,'nbins',25)
addParameter(p,'Nspikesthresh',30)
addParameter(p,'timebinthresh',2) %units: seconds
addParameter(p,'metric',[])
addParameter(p,'SHOWFIG',false)
addParameter(p,'spikemat',[])


parse(p,varargin{:})
intervals = p.Results.intervals;
binsize = p.Results.binsize;
overlap = p.Results.overlap;
nbins = p.Results.nbins;
Nspikesthresh = p.Results.Nspikesthresh;
timebinthresh = p.Results.timebinthresh;
metric = p.Results.metric;
SHOWFIG = p.Results.SHOWFIG;
spikemat = p.Results.spikemat;

%% Calculate spike count matrix
%binsize = 0.2; %s
%overlap = 10;

if isempty(spikemat)
    spikemat = bz_SpktToSpkmat(spikes,'binsize',binsize,'overlap',overlap);
end

%% For each cell, calculate E and I pop rates of all OTHER cells
celltypes = [unique(CellClass.label),'ALL'];
CellClass.ALL = true(size(CellClass.(celltypes{1})));
%E/I Pop Rate: all cells
for tt = 1:length(celltypes)
    spikemat.poprate.(celltypes{tt}) = sum(spikemat.data(:,CellClass.(celltypes{tt})),2);%./...
            %sum(CellClass.(celltypes{tt}))./binsize;
end
%E/I Ratio
spikemat.poprate.EIratio = (spikemat.poprate.pE-spikemat.poprate.pI)./...
    (spikemat.poprate.pE+spikemat.poprate.pI);

%E/I Pop Rate: other cells
numcells = length(spikes.times);
for cc = 1:numcells
    thiscell = false(size(CellClass.pE));
    thiscell(cc) = true;
    spikemat.cellrate{cc} = spikemat.data(:,cc);
    for tt = 1:length(celltypes)
        spikemat.bycellpoprate.(celltypes{tt}){cc} = sum(spikemat.data(:,CellClass.(celltypes{tt}) & ~thiscell),2);%./...
            %sum(CellClass.(celltypes{tt}) & ~thiscell)./binsize;
    end
    
    %E/I Ratio
    spikemat.bycellpoprate.EIratio{cc} = (spikemat.bycellpoprate.pE{cc}-spikemat.bycellpoprate.pI{cc})./...
        (spikemat.bycellpoprate.pE{cc}+spikemat.bycellpoprate.pI{cc});
end

%%
% xwin = bz_RandomWindowInIntervals( intervals,5 );
% figure
% subplot(3,1,1)
% plot(spikemat.timestamps,spikemat.poprate.pE,'k')
% hold on
% plot(spikemat.timestamps,spikemat.poprate.pI,'r')
% plot(spikemat.timestamps,spikemat.poprate.ALL,'k','LineWidth',2)
% 
% xlim(xwin)
% 
% subplot(3,1,2)
% plot(spikemat.timestamps,spikemat.poprate.EIratio,'r')
% xlim(xwin)
%% Calculate E and I pop rate (other cells) for each spike
for tt = 1:length(celltypes)
    spikes.poprate.(celltypes{tt}) = ...
        cellfun(@(X,Y) interp1(spikemat.timestamps,X,Y,'nearest'),...
        spikemat.bycellpoprate.(celltypes{tt}),spikes.times,...
        'UniformOutput',false);
    
end

spikes.poprate.EIratio = ...
    cellfun(@(X,Y) interp1(spikemat.timestamps,X,Y,'nearest'),...
    spikemat.bycellpoprate.EIratio,spikes.times,...
    'UniformOutput',false);

%% Restrict to Specified State
instatespiketimes = cellfun(@(X) InIntervals(X,intervals),...
    spikes.times,'UniformOutput',false);
instateratetimes = InIntervals(spikemat.timestamps,intervals);


%% Cell Type Pop rate Histogram
%popratehist.bins = {unique(spikemat.poprate.pE),unique(spikemat.poprate.pI)};

%Get E/I MUA bins from all E/I MUA
[~,EIMUARateMap.bins{1},EIMUARateMap.bins{2}]...
    = histcounts2(spikemat.poprate.pE(instateratetimes),spikemat.poprate.pI(instateratetimes),nbins);

%FOR NORMALIZATION
%Get E/I MUA Bins for all spikes
[~,~,~,spikes.Ebin,spikes.Ibin] = ...
    cellfun(@(X,Y) histcounts2(X,Y,EIMUARateMap.bins{1},EIMUARateMap.bins{2}),...
    spikes.poprate.pE,spikes.poprate.pI,...
    'UniformOutput',false);
%Get per-cell MUA bin counts
[EIMUARateMap.Nbins] = ...
    cellfun(@(X,Y) histcounts2(X(instateratetimes),Y(instateratetimes),...
    EIMUARateMap.bins{1},EIMUARateMap.bins{2}),...
    spikemat.bycellpoprate.pE,spikemat.bycellpoprate.pI,...
    'UniformOutput',false);

%CALCULATE RATE
%Get per-cell MUA spike counts
[EIMUARateMap.Nspikes] = ...
    cellfun(@(X,Y,Z) histcounts2(X(Z),Y(Z),...
    EIMUARateMap.bins{1},EIMUARateMap.bins{2}),...
    spikes.poprate.pE,spikes.poprate.pI,instatespiketimes,...
    'UniformOutput',false);
%Get Rate in each bin
[EIMUARateMap.rate] = ...
    cellfun(@(X,Y) X./(Y.*binsize),EIMUARateMap.Nspikes,EIMUARateMap.Nbins,...
    'UniformOutput',false);

%Remove bins with not enough time in the bin to calculate rate
for cc = 1:numcells
    EIMUARateMap.rate{cc}(EIMUARateMap.Nbins{cc}<(timebinthresh./spikemat.dt)) = nan;
end


%% EI Ratio Pop rate Histogram
%popratehist.bins = {unique(spikemat.poprate.pE),unique(spikemat.poprate.pI)};

%Get MUA/EI bins from all E/I MUA
[~,EIBalRateMap.bins{1},EIBalRateMap.bins{2}]...
    = histcounts2(spikemat.poprate.EIratio(instateratetimes),spikemat.poprate.ALL(instateratetimes),nbins);

%FOR NORMALIZATION
%Get MUA/EI Bins for all spikes
[~,~,~,spikes.Ebin,spikes.Ibin] = ...
    cellfun(@(X,Y) histcounts2(X,Y,EIBalRateMap.bins{1},EIBalRateMap.bins{2}),...
    spikes.poprate.EIratio,spikes.poprate.ALL,...
    'UniformOutput',false);
%Get per-cell MUA bin counts
[EIBalRateMap.Nbins] = ...
    cellfun(@(X,Y) histcounts2(X(instateratetimes),Y(instateratetimes),...
    EIBalRateMap.bins{1},EIBalRateMap.bins{2}),...
    spikemat.bycellpoprate.EIratio,spikemat.bycellpoprate.ALL,...
    'UniformOutput',false);

%CALCULATE RATE
%Get per-cell MUA spike counts
[EIBalRateMap.Nspikes] = ...
    cellfun(@(X,Y,Z) histcounts2(X(Z),Y(Z),...
    EIBalRateMap.bins{1},EIBalRateMap.bins{2}),...
    spikes.poprate.EIratio,spikes.poprate.ALL,instatespiketimes,...
    'UniformOutput',false);
%Get Rate in each bin
[EIBalRateMap.rate] = ...
    cellfun(@(X,Y) X./(Y.*binsize),EIBalRateMap.Nspikes,EIBalRateMap.Nbins,...
    'UniformOutput',false);

%Remove bins with not enough time in the bin to calculate rate
for cc = 1:numcells
    EIBalRateMap.rate{cc}(EIBalRateMap.Nbins{cc}<(timebinthresh./spikemat.dt)) = nan;
end


%%
if ~isempty(metric)
    %Go through all the E/I bins to get metric
    for ee = 1:length(EIMUARateMap.bins{1})-1
        for ii = 1:length(EIMUARateMap.bins{2})-1

            %Find all the spikes in the bin and in the state
            inbinspikes = cellfun(@(X,Y,Z) X==ee & Y==ii & Z,...
                spikes.Ebin,spikes.Ibin,instatespiketimes,...
                'UniformOutput',false);

            %Cell maps: rate
            for cc = 1:numcells
                EIMUARateMap.metric{cc}(ee,ii) = ...
                    nanmean(metric{cc}(inbinspikes{cc}));

                if sum(inbinspikes{cc}) < Nspikesthresh || EIMUARateMap.Nbins{cc}(ee,ii)<(timebinthresh./spikemat.dt)
                    EIMUARateMap.metric{cc}(ee,ii) = nan;
                end
            end

    %         %All Pop Spike Maps
    %         for tt = 1:length(celltypes)
    %             
    %             allspikeCV2s = cellfun(@(X,Y) X(Y),...
    %                 spikes.CV2(CellClass.(celltypes{tt})),...
    %                 inbinspikes(CellClass.(celltypes{tt})),...
    %                 'UniformOutput',false);
    %             allspikeCV2s = cat(1,allspikeCV2s{:});
    %             
    %             allspikeISIs = cellfun(@(X,Y) X(Y),...
    %                 spikes.ISIs(CellClass.(celltypes{tt})),...
    %                 inbinspikes(CellClass.(celltypes{tt})),...
    %                 'UniformOutput',false);
    %             allspikeISIs = cat(1,allspikeISIs{:});
    %             
    %             popratehist.popstats.meanCV2.(celltypes{tt})(ee,ii) = mean(allspikeCV2s);
    %             popratehist.popstats.meanISI.(celltypes{tt})(ee,ii) = mean(allspikeISIs);
    %             popratehist.popstats.numspikes.(celltypes{tt})(ee,ii) = length(allspikeCV2s);
    %             
    %             if length(allspikeCV2s) < Nspikesthresh
    %                 popratehist.popstats.meanCV2.(celltypes{tt})(ee,ii) = nan;
    %                 popratehist.popstats.meanISI.(celltypes{tt})(ee,ii) = nan;
    %             end
    %             %popratehist.popstats.meanrate.(celltypes{tt}) = nanmean(1./cat(3,popratehist.meanISI_bycell{CellClass.(celltypes{tt})}),3);
    %         end
    %     end
        end
    end
end


%% Metric Map: MUA vs EI
if ~isempty(metric)
    %Go through all the E/I bins to get metric
    for ee = 1:length(EIBalRateMap.bins{1})-1
        for ii = 1:length(EIBalRateMap.bins{2})-1

            %Find all the spikes in the bin and in the state
            inbinspikes = cellfun(@(X,Y,Z) X==ee & Y==ii & Z,...
                spikes.Ebin,spikes.Ibin,instatespiketimes,...
                'UniformOutput',false);

            %Cell maps: rate
            for cc = 1:numcells
                EIBalRateMap.metric{cc}(ee,ii) = ...
                    nanmean(metric{cc}(inbinspikes{cc}));

                if sum(inbinspikes{cc}) < Nspikesthresh || EIBalRateMap.Nbins{cc}(ee,ii)<(timebinthresh./spikemat.dt)
                    EIBalRateMap.metric{cc}(ee,ii) = nan;
                end
            end

        end
    end
end
%% Cell type averages
for tt = 1:length(celltypes)
    EIMUARateMap.meanrate.(celltypes{tt}) = nanmean(cat(3,EIMUARateMap.rate{CellClass.(celltypes{tt})}),3);
    EIMUARateMap.meanmetric.(celltypes{tt}) = nanmean(cat(3,EIMUARateMap.metric{CellClass.(celltypes{tt})}),3);
    EIMUARateMap.geomeanrate.(celltypes{tt}) = exp(nanmean(log(cat(3,EIMUARateMap.rate{CellClass.(celltypes{tt})})),3));

    EIBalRateMap.meanrate.(celltypes{tt}) = nanmean(cat(3,EIBalRateMap.rate{CellClass.(celltypes{tt})}),3);
    EIBalRateMap.meanmetric.(celltypes{tt}) = nanmean(cat(3,EIBalRateMap.metric{CellClass.(celltypes{tt})}),3);
    EIBalRateMap.geomeanrate.(celltypes{tt}) = exp(nanmean(log(cat(3,EIBalRateMap.rate{CellClass.(celltypes{tt})})),3));

    
end

%% OUTPUT STUFF
poprate = spikes.poprate;


%%
if SHOWFIG
    figure
    for tt = 1:length(celltypes)
        subplot(3,3,tt)
            h = imagesc(EIMUARateMap.bins{1}./(sum(CellClass.pE)-1)./binsize,...
                EIMUARateMap.bins{2}./(sum(CellClass.pI))./binsize,...
                log10(EIMUARateMap.meanrate.(celltypes{tt}))');
            set(h,'AlphaData',~isnan(EIMUARateMap.meanrate.(celltypes{tt})'));
            xlabel('pE Rate (Hz/cell)');ylabel('pI Rate (Hz/cell)')
            title(['Rate - ',celltypes{tt},' cells'])
            axis xy
            colorbar
            %caxis([-1.5 0])
            LogScale('c',10)

        subplot(3,3,tt+3)
            h = imagesc(EIBalRateMap.bins{1},...
                EIBalRateMap.bins{2}./(sum(CellClass.ALL))./binsize,...
                log10(EIBalRateMap.meanrate.(celltypes{tt}))');
            set(h,'AlphaData',~isnan(EIBalRateMap.meanrate.(celltypes{tt})'));
            xlabel('EI Bal');ylabel('MUA (Hz/cell)')
            title(['Rate - ',celltypes{tt},' cells'])
            axis xy
            colorbar
            %caxis([-1.5 0])
            LogScale('c',10)
    end
    %%
    figure
    for tt = 1:length(celltypes)
        subplot(3,3,tt)
            h = imagesc(EIMUARateMap.bins{1}./(sum(CellClass.pE)-1)./binsize,...
                EIMUARateMap.bins{2}./(sum(CellClass.pI))./binsize,...
                (EIMUARateMap.meanmetric.(celltypes{tt}))');
            set(h,'AlphaData',~isnan(EIMUARateMap.meanmetric.(celltypes{tt})'));
            xlabel('pE Rate (Hz/cell)');ylabel('pI Rate (Hz/cell)')
            title(['Metric - ',celltypes{tt},' cells'])
            axis xy
            colorbar
            %caxis([-1.5 0])
           % LogScale('c',10)

        subplot(3,3,tt+3)
            h = imagesc(EIBalRateMap.bins{1},...
                EIBalRateMap.bins{2}./(sum(CellClass.ALL))./binsize,...
                (EIBalRateMap.meanmetric.(celltypes{tt}))');
            set(h,'AlphaData',~isnan(EIBalRateMap.meanmetric.(celltypes{tt})'));
            xlabel('EI Bal');ylabel('MUA (Hz/cell)')
            title(['Metric - ',celltypes{tt},' cells'])
            axis xy
            colorbar
            %caxis([-1.5 0])
           % LogScale('c',10)
    end

end

%%

% binthreshold = 500;
% for tt = 1:length(celltypes)
%     popratehist.meancellstats.meanCV2.(celltypes{tt}) = nanmean(cat(3,popratehist.meanCV2_bycell{CellClass.(celltypes{tt})}),3);
%     popratehist.meancellstats.meanrate.(celltypes{tt}) = nanmean(1./cat(3,popratehist.meanISI_bycell{CellClass.(celltypes{tt})}),3);
%     
%     popratehist.Nbins_all = sum(cat(3,popratehist.Nbins{CellClass.(celltypes{tt})}),3);
%     popratehist.popstats.meanrate.(celltypes{tt}) = popratehist.popstats.numspikes.(celltypes{tt})./popratehist.Nbins_all./spikemat.dt;
%     popratehist.popstats.meanrate.(celltypes{tt})(popratehist.Nbins_all<binthreshold) = nan;
% end
% %popratehist.popstats.meanCV2.pI = mean(cat(3,popratehist.meanCV2{CellClass.pI}),3);

%%
%EIMUARateMap.Ecellbins = 
 %%
% figure
% for tt = 1:length(celltypes)
% subplot(2,2,tt)
%     imagesc(log10(popratehist.meanrate.(celltypes{tt}))')
%     colorbar
%     axis xy
%     title(celltypes{tt})
%     LogScale('c',10)
% end
% 
% for tt = 1:length(celltypes)
% subplot(2,2,tt+2)
%     imagesc((popratehist.meanmetric.(celltypes{tt}))')
%     colorbar
%     axis xy
%     title(celltypes{tt})
% end
%     %%
%     figure
%     subplot(2,2,1)
%         imagesc(popratehist.metric{12}')
%         axis xy
%         colorbar
%         
%     subplot(2,2,2)
%         imagesc(popratehist.bins{1},popratehist.bins{2},log10(popratehist.rate{10})')
%         axis xy
%         colorbar
end

