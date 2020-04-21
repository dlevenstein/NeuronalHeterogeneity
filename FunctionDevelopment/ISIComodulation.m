function [outputArg1,outputArg2] = ISIComodulation(spiketimes)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%%
state = states{1};
ISIStats.allspikes.instate = cellfun(@(X) InIntervals(X,double(SleepState.ints.(state))),...
    ISIStats.allspikes.times,'UniformOutput',false);

%%
usespikes.times = cellfun(@(X,Y) X(Y), ISIStats.allspikes.times,ISIStats.allspikes.instate,'UniformOutput',false);
usespikes.ISIs = cellfun(@(X,Y) X(Y), ISIStats.allspikes.ISIs,ISIStats.allspikes.instate,'UniformOutput',false);
usespikes.ISInp1 = cellfun(@(X,Y) X(Y), ISIStats.allspikes.ISInp1,ISIStats.allspikes.instate,'UniformOutput',false);


twin = 0.1; %(s)
ISItol = 0.1; %log(s)
numspikethresh = 50; %Need this many spikes to calculate hist (i or j), otherwise nan.

numjbins = 200;
jbins = linspace(ISIStats.ISIhist.logbins(1),ISIStats.ISIhist.logbins(end),numjbins);


i_ISIhist = nan(length(jbins),length(ISIStats.ISIhist.logbins),spikes.numcells,spikes.numcells);
i_ISIhist_0bin = nan(length(ISIStats.ISIhist.logbins),spikes.numcells,spikes.numcells); 
i_ISIhist_any = nan(length(ISIStats.ISIhist.logbins),spikes.numcells,spikes.numcells); 

%DEV
for ii = 1
    %bz_Counter(ii,spikes.numcells,'Cell i')
    %
    for jj = 1:spikes.numcells
        bz_Counter(jj,spikes.numcells,['Cell i: ',num2str(ii)])

        %Get the marginal ISI distribution for cell i
        [i_ISIhist_marj] = hist(log10([usespikes.ISIs{ii};usespikes.ISInp1{ii}]),...
            ISIStats.ISIhist.logbins);
        i_ISIhist_marj(i_ISIhist_marj<numspikethresh)=nan;
        i_ISIhist_marj = i_ISIhist_marj./nansum(i_ISIhist_marj); %Normalize
        
        parfor bb = 1:numjbins %For each ISI bin,

            %find cell j spikes in that bin.   
            thisbinspikes = ismembertol(log10(usespikes.ISIs{jj}),jbins(bb),ISItol,'DataScale', 1)...
                 |ismembertol(log10(usespikes.ISInp1{jj}),jbins(bb),ISItol,'DataScale', 1);

            if sum(thisbinspikes)<numspikethresh
               continue
            end

            %Find all cell j spikes within twin of those spikes
            cellispikes = ismembertol(usespikes.times{ii},...
                usespikes.times{jj}(thisbinspikes),twin,'DataScale', 1);
            if sum(cellispikes)<numspikethresh
               continue
            end

            %and histogram of their ISIs (prev/next)
            [i_ISIhist(bb,:,ii,jj)] = ...
                hist(log10([usespikes.ISIs{ii}(cellispikes);usespikes.ISInp1{ii}(cellispikes)]),...
                ISIStats.ISIhist.logbins);
            i_ISIhist(bb,:,ii,jj) = i_ISIhist(bb,:,ii,jj)./sum(i_ISIhist(bb,:,ii,jj)); %Normalize
            i_ISIhist(bb,:,ii,jj) = i_ISIhist(bb,:,ii,jj)./i_ISIhist_marj;
        end

        % Here: find i spikes that aren't within tolerance of any j spikes
        cellispikes_0bin = ~ismembertol(usespikes.times{ii},...
            usespikes.times{jj},twin,'DataScale', 1);
        [i_ISIhist_0bin(:,ii,jj)] = hist(log10([usespikes.ISIs{ii}(cellispikes_0bin);...
            usespikes.ISInp1{ii}(cellispikes_0bin)]),...
            ISIStats.ISIhist.logbins);
        i_ISIhist_0bin(:,ii,jj) = i_ISIhist_0bin(:,ii,jj)./sum(i_ISIhist_0bin(:,ii,jj)); %Normalize
        i_ISIhist_0bin(:,ii,jj) = i_ISIhist_0bin(:,ii,jj)./i_ISIhist_marj';
        
        [i_ISIhist_any(:,ii,jj)] = hist(log10([usespikes.ISIs{ii}(~cellispikes_0bin);...
            usespikes.ISInp1{ii}(~cellispikes_0bin)]),...
            ISIStats.ISIhist.logbins);
        i_ISIhist_any(:,ii,jj) = i_ISIhist_any(:,ii,jj)./sum(i_ISIhist_any(:,ii,jj)); %Normalize
        i_ISIhist_any(:,ii,jj) = i_ISIhist_any(:,ii,jj)./i_ISIhist_marj';

    end
end
end

