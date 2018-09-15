function [PowerPhaseRatemap,spikebinIDs] = bz_PowerPhaseRatemap(spikes,filteredLFP,varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%
%INPUTS
%   spikes          structure containing spikes.times, from bz_getSpikes
%   filteredLFP     structure containing filteredLFP.timestamps, 
%                       filteredLFP.amp, filteredLFP.phase,
%                       filteredLFP.samplingRate. from bz_Filter(lfp).
% 
%  (options)
%   'ints'
%   'powernorm'     (default: 'modZ')
%   'ratenorm'      (default: 'none')
%   'numbins'       (default: 20)
%
%OUTPUTS
%   PowerPhaseRatemap
%   spikebinIDs
%
%DLevenstein 2018
%% DEV
p = inputParser;
addParameter(p,'ints',[0 Inf],@isnumeric)
addParameter(p,'powernorm','Zlog')

parse(p,varargin{:})
ints = p.Results.ints;
powernorm = p.Results.powernorm;

%% Find lfp and spikes in the intervals
instatespiketimes = cellfun(@(X) InIntervals(X,ints),...
    spikes.times,'UniformOutput',false);
instateLFPtimes = InIntervals(filteredLFP.timestamps,ints);


%%

%Normalize Power
switch powernorm
    case 'Zlog'
        filteredLFP.amp = NormToInt(log10(filteredLFP.amp),'Z',ints,filteredLFP.samplingRate);
end

%%

%Get Power/Phase at each spike
spikes.amp = cellfun(@(X) interp1(filteredLFP.timestamps,filteredLFP.amp,X,'nearest'),...
    spikes.times,'uniformoutput',false);
spikes.phase = cellfun(@(X) interp1(filteredLFP.timestamps,filteredLFP.phase,X,'nearest'),...
    spikes.times,'uniformoutput',false);

%%
numbins = 20;

%Power/Phase Bins
phaseedges = linspace(-pi,pi,numbins+1);
PowerPhaseRatemap.phasebins=phaseedges(1:end-1)+0.5.*diff(phaseedges(1:2));
poweredges = linspace(-1.75,1.75,numbins+1);
PowerPhaseRatemap.powerbins=poweredges(1:end-1)+0.5.*diff(poweredges(1:2));
poweredges(1) = -Inf; poweredges(end) = Inf;

%Calculate Power/Phase Spike Histogram
%phaseamphist = cellfun(@(X,Y) hist3([X,Y],{powerbins,phasebins}),spikeamp,spikephase,'UniformOutput',false);

[PowerPhaseRatemap.Nspikes,~,~,spikebinIDs.powerbin,spikebinIDs.phasebin] = ...
    cellfun(@(X,Y,Z) histcounts2(X(Z),Y(Z),poweredges,phaseedges),...
    spikes.amp,spikes.phase,instatespiketimes,...
    'UniformOutput',false);

%Calculate Power/Phase Time Histogram
%phaseamphist_t = hist3([t_amp,t_phase],{powerbins,phasebins});
%phaseamphist_t = phaseamphist_t./sf_LFP;

PowerPhaseRatemap.occupancy = ...
    histcounts2(filteredLFP.amp(instateLFPtimes),...
    filteredLFP.phase(instateLFPtimes),poweredges,phaseedges);
PowerPhaseRatemap.occupancy = ...
    PowerPhaseRatemap.occupancy./filteredLFP.samplingRate;
PowerPhaseRatemap.powerdist = sum(PowerPhaseRatemap.occupancy,2);

%Normalize Rate
%totaltime = sum(diff(ints,1,2));
%meanrate = cellfun(@(X) length(X)./totaltime,spiketimes,'UniformOutput',false);
PowerPhaseRatemap.ratemap = cellfun(@(X) X./PowerPhaseRatemap.occupancy,...
    PowerPhaseRatemap.Nspikes,'UniformOutput',false);
%phaseamprate = cellfun(@(X,Y) X./Y,phaseamprate,meanrate,'UniformOutput',false);

PowerPhaseRatemap.meanrate = nanmean(cat(3,PowerPhaseRatemap.ratemap{:}),3);

%%
%excell = randsample(spikes.numcells,1);
figure
% subplot(2,2,1)
%     imagesc(PowerPhaseRatemap.phasebins,PowerPhaseRatemap.powerbins,...
%         PowerPhaseRatemap.ratemap{excell})
%     hold on
%     imagesc(PowerPhaseRatemap.phasebins+2*pi,PowerPhaseRatemap.powerbins,...
%         PowerPhaseRatemap.ratemap{excell})
%     xlim([-pi 3*pi])
%     axis xy
%     colorbar
    
subplot(2,2,1)
    imagesc(PowerPhaseRatemap.phasebins,PowerPhaseRatemap.powerbins,...
        PowerPhaseRatemap.meanrate)
    hold on
    imagesc(PowerPhaseRatemap.phasebins+2*pi,PowerPhaseRatemap.powerbins,...
        PowerPhaseRatemap.meanrate)
    plot(linspace(-pi,3*pi,100),cos(linspace(-pi,3*pi,100)),'k')
    xlim([-pi 3*pi])
    axis xy
    colorbar  
    xlabel('Phase');
    ylabel('Power (Z(log))')
    
subplot(2,2,3)
    bar(PowerPhaseRatemap.powerbins,PowerPhaseRatemap.powerdist)
    xlabel('Power (Z(log))');ylabel('Time (s)')
    box off
    axis tight


end

