function [ ConditionalLFPCoupling ] = bz_ConditionalLFPCoupling( spikes,condition,filtLFP,varargin )
%bz_ConditionalLFPCoupling( spikes,condition,filtLFP ) calculates the
%phase-coupling of spikes to the LFP conditioned on some property of the
%spike (for example, ISI or CV2)
%
%INPUTS
%   spikes      structure (i.e. from bz_GetSpikes) 
%                   spikes.times: cell array with spiketimes for each cell  
%   condition   cell array with conditional property of each spike
%               -or- (add this)
%               structure with 
%                   condition.data
%                   condition.timestamps
%   filtLFP     structure (i.e. from bz_WaveSpec)
%                   filtLFP.data: complex valued LFP (i.e. hilbert or wavelet transorm)
%                   filtLFP.freqs (if multiple)
%                   filtLFP.timestamps
%
%   (options)
%   'intervals' only calculate coupling in some time intervals
%   'numXbins'  number of bins for your conditional variable (default 60)
%   'Xbounds'   bounds of your conditional variable
%   'minX'      minumum number of spikes to calculate coupling (default 25)
%   'spikeLim'  limit number of spikes to look at for each cell 
%               (randomly omits spikes, default: Inf)
%   'showFig'   true/false
%   'saveFig'   folder in which to save the figure
%   'figName'   default: 'CondLFPCouping'
%   'baseName'  (for figure saving)
%
%
%OUTPUT
%
%
%DLevenstein 2019
%%
p = inputParser;
addParameter(p,'numXbins',75)
addParameter(p,'Xbounds',[])
addParameter(p,'minX',25)
addParameter(p,'intervals',[0 Inf])
addParameter(p,'showFig',false)
addParameter(p,'CellClass',[])
addParameter(p,'saveFig',false)
addParameter(p,'figName','CondLFPCouping')
addParameter(p,'baseName',[])
addParameter(p,'spikeLim',Inf)
parse(p,varargin{:})
numXbins = p.Results.numXbins;
Xbounds = p.Results.Xbounds;
minX = p.Results.minX;
ints = p.Results.intervals;
SHOWFIG = p.Results.showFig;
CellClass = p.Results.CellClass;
saveFig = p.Results.saveFig;
figName = p.Results.figName;
baseName = p.Results.baseName;
spikeLim = p.Results.spikeLim;


%% Restrict spikes and lfp to the interval
filtLFP.inint = InIntervals(filtLFP.timestamps,ints);
filtLFP.timestamps = filtLFP.timestamps(filtLFP.inint);
filtLFP.data = filtLFP.data(filtLFP.inint,:);

spikes.inint = cellfun(@(X) InIntervals(X,ints),spikes.times,'UniformOutput',false);

%Apply spikelimit here
spikes.toomany = cellfun(@(X) (sum(X)-spikeLim).*((sum(X)-spikeLim)>0),spikes.inint);
spikes.numcells = length(spikes.times);
for cc = 1:spikes.numcells
    spikes.inint{cc}(randsample(find(spikes.inint{cc}),spikes.toomany(cc)))=false;
end

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
    disp('Interpolating LFP at each spike...')
    disp('If you have a lot of cells/spikes/freqs this can take a few minutes.')
    disp('If this is prohibitive (time or RAM), try using ''spikeLim''')
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

%%
%For each bin calculate - mean power of spikes, power-weighted MRL (in
%state)
%Mean Y given X
filtLFP.nfreqs = size(filtLFP.data,2);
meanpower = nan(numXbins,filtLFP.nfreqs,spikes.numcells);
mrl = nan(numXbins,filtLFP.nfreqs,spikes.numcells);
mrlangle = nan(numXbins,filtLFP.nfreqs,spikes.numcells);
for xx = 1:length(Xbins)
    %Mean power at spike
    meanpowertemp = cellfun(@(lfp,binID) ...
        nanmean(abs(lfp(binID==xx,:))),...
        spikes.LFP,spikes.XbinID,'UniformOutput',false);
    
    %Mean resultant vector
    pMRVtemp = cellfun(@(lfp,binID) ...
        nanmean(abs(lfp(binID==xx,:)).*exp(1i.*angle(lfp(binID==xx,:)))),...
        spikes.LFP,spikes.XbinID,'UniformOutput',false);
   
   for cc = 1:spikes.numcells
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




%% Mutual Information
powerbins = linspace(-0.5,0.5,10);
for cc = 1:spikes.numcells
    cc
    for ff = 1:filtLFP.nfreqs
        

        joint = hist3([spikes.condition{cc} log10(abs(spikes.LFP{cc}(:,ff)))],{Xbins,powerbins});
        joint = joint./sum(joint(:));

        margX = hist(spikes.condition{cc},Xbins);
        margX = margX./sum(margX);
        margPower = hist(log10(abs(spikes.LFP{cc}(:,ff))),powerbins);
        margPower = margPower./sum(margPower);
        jointindependent =  bsxfun(@times, margX.', margPower);


        mutXPow = joint.*log2(joint./jointindependent); % mutual info at each bin
        totmutXPow(cc,ff) = nansum(mutXPow(:)); % sum of all mutual information 
    end
end



if ~isempty(CellClass)
    celltypes = unique(CellClass.label);
    for tt = 1:length(celltypes)
        groupmutinf.(celltypes{tt}) = nanmean(totmutXPow(CellClass.(celltypes{tt}),:),1);
    end
else
    groupmutinf.ALL = nanmean(totmutXPow,1);
    celltypes={'ALL'};
end
%%
figure
subplot(2,2,1)
    imagesc(log2(filtLFP.freqs),[1 spikes.numcells],totmutXPow)
    LogScale('x',2)
    xlabel('Freq (Hz)');ylabel('Cell')
subplot(2,2,2)
hold on
for tt = 1:length(celltypes)
    plot(log2(filtLFP.freqs),groupmutinf.(celltypes{tt}))
end
    LogScale('x',2)
%%
figure
subplot(2,2,1)
imagesc(Xbins,powerbins,joint')
axis xy
subplot(2,2,2)
imagesc(Xbins,powerbins,jointindependent')
axis xy

subplot(2,2,3)
imagesc(Xbins,powerbins,mutXPow')
axis xy
colorbar

%%
if SHOWFIG
    
    % Separate cell classes
    if ~isempty(CellClass)
        celltypes = unique(CellClass.label);
        for tt = 1:length(celltypes)
            allmeanpower.(celltypes{tt}) = nanmean(meanpower(:,:,CellClass.(celltypes{tt})),3);
            almeanpMRL.(celltypes{tt}) = nanmean(mrl(:,:,CellClass.(celltypes{tt})),3);
        end
    else
        allmeanpower.ALL = nanmean(meanpower,3);
        almeanpMRL.ALL = nanmean(mrl,3);
        celltypes={'ALL'};
    end

    
    powermap = makeColorMap([0 0 0.8],[1 1 1],[0.8 0 0]);
    figure
    for tt = 1:length(celltypes)
    subplot(3,2,tt*2)
    colormap(gca,powermap)
        imagesc(Xbins,log2(filtLFP.freqs), (allmeanpower.(celltypes{tt}))')
        colorbar
        ColorbarWithAxis([0.5 1.5],'Power (mean^-^1)')
        LogScale('x',10);
        LogScale('y',2)
        axis xy
        xlabel('ISI (s)');ylabel('freq (Hz)')
        title((celltypes{tt}))
    end  
    
    
    for tt = 1:length(celltypes)
    subplot(3,2,tt*2-1)
        imagesc(Xbins,log2(filtLFP.freqs), almeanpMRL.(celltypes{tt})')
        colorbar
        hold on
        %caxis([0.5 1.5])
        LogScale('x',10);
        LogScale('y',2)
        ColorbarWithAxis([0 0.4],'Phase Coupling (pMRL)')

        axis xy
        xlabel('ISI (s)');ylabel('freq (Hz)')
        title((celltypes{tt}))
    end 
    
    if saveFig
        try
            NiceSave(figName,saveFig,baseName,'includeDate',true)
        catch
            disp('Sorry, I wasn''t able to save your figure :''(')
        end
    end
end


%% Output
ConditionalLFPCoupling.Xbins = Xbins;
ConditionalLFPCoupling.freqs = filtLFP.freqs;
ConditionalLFPCoupling.meanpower = meanpower;
ConditionalLFPCoupling.mrl = mrl;
ConditionalLFPCoupling.mrlangle = mrlangle;

%%
clear spikes
clear filtLFP
%%
%ADD: condition distribution given power, power disirbution given condition
%Mutual information? Conditional Entropy? Pwoer distribution or even
%power/phase distribution?
end

