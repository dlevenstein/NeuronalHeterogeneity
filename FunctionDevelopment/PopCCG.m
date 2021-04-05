function [popCCG] = PopCCG(spikes,varargin)
%popCCG = PopCCG(spikes,varargin) calculates the statistics 
%inter-spike intervals for the spiketimes in spikes.
%
%   INPUTS
%       spikes      Structure with spikes.times and spikes.UID
%
%       (options)
%       'ints'        
%       'cellclass'     Cell array of strings - label for each cell. 
%       'classnames'  (To expect)
%       'ignorepairs' [ncells x ncells] matrix where entry (i,j)=true means to
%                       ignore ccg(i,j): events from cell j with respect to
%                       reference events from cell i



%return: (sort by pop coupling (chorister/soloist) and value of coupling and
%timescale of coupling?
%%
%defaultstates.ALL = [-Inf Inf];

% parse args
p = inputParser;
addParameter(p,'ints',[-Inf Inf])
addParameter(p,'showfig',false,@islogical);
addParameter(p,'cellclass',[]);
addParameter(p,'binsize',0.001);
addParameter(p,'duration',0.4);
addParameter(p,'sortcells',[]);
addParameter(p,'classnames',[]);
addParameter(p,'figfolder',false)
addParameter(p,'figname',[])
addParameter(p,'ignorepairs',[])
addParameter(p,'minspikes',[])


parse(p,varargin{:})
ints = p.Results.ints;
cellclass = p.Results.cellclass;
SHOWFIG = p.Results.showfig;
binsize = p.Results.binsize;
duration = p.Results.duration;
sortcells = p.Results.sortcells;
classnames = p.Results.classnames;
figfolder = p.Results.figfolder;
figname = p.Results.figname;
ignorepairs = p.Results.ignorepairs;
minSpikes = p.Results.minspikes;

%%
spikes.instate = cellfun(@(X) InIntervals(X,ints),spikes.times,'UniformOutput',false);
spikes.numcells = length(spikes.times);
%% CCGs
ccgspikes = cellfun(@(X,Y) X(Y),spikes.times,spikes.instate,'UniformOutput',false);
numspikes = cellfun(@length,ccgspikes);
[allccg,popCCG.t_ccg] = CCG(ccgspikes,[],'binSize',binsize,'duration',duration,'norm','rate'); 

%%
for cc = 1:spikes.numcells
    othercells = true(size(spikes.times));
    othercells(cc) = false;
    othercells(ignorepairs(cc,:)) = false;
    %othercells(numspikes<minSpikes) = false;
    popCCG.cells.ALL(:,cc) = mean(allccg(:,cc,othercells),3);
end
%%
figure
imagesc(popCCG.t_ccg,[0 spikes.numcells],popCCG.cells.ALL')
colorbar
caxis([0 1])
%% Cell Class Populations

%Make the cell-type specific sortings and average distributions
if ~isempty(cellclass)
    %Check for empty cell class entries
    noclass = cellfun(@isempty,cellclass);
    sorts.numclassycells = sum(~noclass);
    %cellclass(noclass)={'none'};
    if isempty(classnames)
        classnames = unique(cellclass(~noclass));
    end
    numclasses = length(classnames);
    for tt = 1:numclasses
        inclasscells{tt} = strcmp(classnames{tt},cellclass);
        %CellClass.(classnames{tt}) = inclasscells{tt};
        popspikes{tt} = cat(1,ccgspikes{inclasscells{tt}});
    end
end

%%
clear meanCCG
for cc = 1:spikes.numcells
    for tt = 1:numclasses
        popothercells = inclasscells{tt};
        popothercells(cc) = false;
        popothercells(ignorepairs(cc,:)) = false;
        %popothercells(numspikes<minSpikes) = false;
        popCCG.cells.(classnames{tt})(:,cc) = nanmean(allccg(:,cc,popothercells),3);
    end
end

for tt = 1:numclasses
    for tt2 = 1:numclasses
        popothercells = inclasscells{tt2};
        popothercells(numspikes<minSpikes) = false;
        popCCG.cellsmean.(classnames{tt})(:,tt2) = nanmean(popCCG.cells.(classnames{tt})(:,popothercells),2);
    end
end
%%
%popspikes = {cat(1,ccgspikes{CellClass.pE}),cat(1,ccgspikes{CellClass.pI})};
[popccg,~] = CCG(popspikes,[],'binSize',binsize,'duration',duration,'norm','rate'); 
for tt = 1:numclasses
    popCCG.pop.(classnames{tt}) = popccg(:,:,tt)./sum(inclasscells{tt});
end
%%
if SHOWFIG || figfolder

if isempty(sortcells)
   sortcells = [1:length(spikes.times)]; 
end
figure
for tt = 1:numclasses
    subplot(3,3,tt)
        imagesc(popCCG.t_ccg,[0 spikes.numcells],(popCCG.cells.(classnames{tt})(:,sortcells))')
        hold on
        
        plot(popCCG.t_ccg,bz_NormToRange(-popccg(:,1,tt),[1 sum(inclasscells{1})-1]),'linewidth',2);%,'color',cellcolor{tt})

        try
        for tt2 = 2:numclasses
            plot(xlim(gca),sum(inclasscells{tt2-1}).*[1 1],'w')
            plot(popCCG.t_ccg,bz_NormToRange(-popccg(:,tt2,tt),sum(inclasscells{tt2-1})+[1 sum(inclasscells{tt2})-1]),'linewidth',2);%,'color',cellcolor{tt})
        end
        catch
           disp('Only one class') 
        end
        %plot(t_ccg,bz_NormToRange(-popccg(:,tt,2),sum(CellClass.pE)+[1 sum(CellClass.pI)]),'color',cellcolor{tt})
        %plot(xlim(gca),sum(CellClass.pE).*[1 1],'w')
        colorbar
        title(classnames{tt})
        xlabel('t lag')
        %caxis([0 2])
end

if figfolder
    if ~isempty(figname)
        baseName = figname;
    end
    NiceSave(['PopCCG'],figfolder,baseName);
end
end


end

