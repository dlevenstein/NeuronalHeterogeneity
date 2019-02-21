function [ output_args ] = bz_ConditionalLFPCoupling( spikes,condition,filtLFP )
%bz_ConditionalLFPCoupling( spikes,condition,filtLFP ) calculates the
%phase-coupling of spikes to the LFP conditioned on some property of the
%spike (for example, ISI or CV2)
%
%INPUTS
%   spikes      structure with 
%                   spikes.times: cell array with spiketimes for each cell
%                   (from bz_GetSpikes)
%   condition   cell array with conditional property of each spike
%               -or- (add this)
%               structure with 
%                   condition.data
%                   condition.timestamps
%   filtLFP     structure with
%                   filtLFP.data: complex valued LFP (i.e. hilbert or wavelet transorm)
%                   filtLFP.timestamps
%
%   (options)
%   'intervals'
%
%
%OUTPUT
%
%
%DLevenstein 2019
%%
Xbounds = [-2.5 1];
numXbins = 70;
minX = 50;
intervals = 

%% Restrict spikes and lfp to the interval
filtLFP.inint = InIntervals(filtLFP.timestamps,ints);
filtLFP.timestamps = filtLFP.timestamps(filtLFP.inint);
filtLFP.data = filtLFP.data(filtLFP.inint,:);

spikes.inint = cellfun(@(X) InIntervals(X,ints),spikes.times,'UniformOutput',false);
spikes.times = cellfun(@(X,Y) X(Y),spikes.times,spikes.inint,'UniformOutput',false);
spikes.condition = cellfun(@(X,Y) X(Y),condition,spikes.inint,'UniformOutput',false);
clear condition
%%
%Mean-Normalize power
filtLFP.meanpower = mean(abs(filtLFP.data),1);
filtLFP.data = bsxfun(@(X,Y) X./Y,filtLFP.data,filtLFP.meanpower);

%% Get Power/Phase at each spike

%Get complex-valued filtered LFP at each spike time
if spikes.numcells>50
    display('Interpolating LFP at each spike... if you have a lot of cells this can take a few minutes')
end
for cc = 1:spikes.numcells
    spikes.LFP{cc} = interp1(filtLFP.timestamps,filtLFP.data,spikes.times{cc},'nearest');
end

%%

%Divide condition into bin (as in ConditionalHist)
Xedges = linspace(Xbounds(1),Xbounds(2),numXbins+1);
Xbins = Xedges(1:end-1)+ 0.5.*diff(Xedges([1 2]));
Xedges(1) = -inf;Xedges(end) = inf;

%First calculate the marginal probability of the conditional (X), 
%get the bin in which each spike falls 
[Xhist,~,spikes.XbinID] = cellfun(@(X) histcounts(X,Xedges),spikes.condition,'UniformOutput',false);

%Xhist4norm = Xhist;Xhist4norm(Xhist4norm<=minX) = nan;

%%

%%
%For each bin calculate - mean power of spikes, power-weighted MRL (in
%state)
%Mean Y given X
clear meanpower
clear mrl
clear mrlangle
for xx = 1:length(Xbins)
    %Mean power at spike
    meanpowertemp = cellfun(@(lfp,binID) ...
        nanmean(abs(lfp(binID==xx,:))),...
        spikes.LFP,spikes.XbinID,'UniformOutput',false);
    
    %Mean resultant vector
    pMRVtemp = cellfun(@(lfp,binID) ...
        nanmean(abs(lfp(binID==xx,:)).*exp(1i.*angle(lfp(binID==xx,:)))),...
        ISIStats.allspikes.LFP,ISIStats.allspikes.XbinID,'UniformOutput',false);
   
   for cc = 1:length(ISIStats.allspikes.times)
       meanpower(xx,:,cc)=meanpowertemp{cc};
       mrl(xx,:,cc)=abs(pMRVtemp{cc});
       mrlangle(xx,:,cc)=angle(pMRVtemp{cc});
       if Xhist{cc}(xx)<minX
           meanpower(xx,:,cc) = nan;
           mrl(xx,:,cc) = nan;
           mrlangle(xx,:,cc) = nan;
       end
           
   end
end
%%
meanpower(:,:,CellClass.(celltypes{tt}))
mrl(:,:,CellClass.(celltypes{tt}))

%%
%%
powermap = makeColorMap([0 0 0.8],[1 1 1],[0.8 0 0]);
figure
for tt = 1:length(celltypes)
subplot(3,3,tt*3)
colormap(gca,powermap)
    imagesc(Xbins,log2(filtLFP.freqs), conditionalpower.(celltypes{tt})')
    colorbar
    caxis([0.4 1.6])
    LogScale('x',10);
    LogScale('y',2)
    axis xy
    xlabel('ISI (s)');ylabel('freq (Hz)')
    title((celltypes{tt}))
end   
    
%%
%ADD: condition distribution given power, power disirbution given condition
%Mutual information? Conditional Entropy? Pwoer distribution or even
%power/phase distribution?
end

