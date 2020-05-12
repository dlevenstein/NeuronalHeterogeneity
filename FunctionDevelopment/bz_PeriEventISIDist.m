function [PeriEventISIDist] = bz_PeriEventISIDist(spikes,eventtimes,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%
p = inputParser;
addParameter(p,'savecellinfo',false,@islogical)
addParameter(p,'basePath',pwd,@isstr)
addParameter(p,'figfolder',false)
addParameter(p,'figname',[])
addParameter(p,'showfig',false,@islogical);
addParameter(p,'cellclass',[]);
addParameter(p,'forceRedetect',false,@islogical);
addParameter(p,'shuffleCV2',false,@islogical);
addParameter(p,'numISIbins',100);
addParameter(p,'numCV2bins',50);
addParameter(p,'ISIbounds',[0.001 100]);
addParameter(p,'minX',50);
addParameter(p,'numXbins',50);
addParameter(p,'winsize',[-1 1]);
addParameter(p,'whichISIs','both');
%addParameter(p,'fitGammas',false);


parse(p,varargin{:})
cellclass = p.Results.cellclass;
basePath = p.Results.basePath;
SAVECELLINFO = p.Results.savecellinfo;
figfolder = p.Results.figfolder;
figname = p.Results.figname;
SHOWFIG = p.Results.showfig;
forceRedetect = p.Results.forceRedetect;
SHUFFLECV2 = p.Results.shuffleCV2;
numISIbins = p.Results.numISIbins;
numCV2bins = p.Results.numCV2bins;
%fitGammas = p.Results.fitGammas;
logISIbounds = log10(p.Results.ISIbounds);
minX = (p.Results.minX);
numXbins = (p.Results.numXbins);
winsize = (p.Results.winsize);
whichISIs = (p.Results.whichISIs);
%%
if strcmp(cellclass,'load')
    cellclass = bz_LoadCellinfo(basePath,'CellClass');
    cellclass = cellclass.label;
end

%%
if iscell(spikes)
    temp = cellfun(@(X) bz_PeriEventISIDist(X,eventtimes,varargin{:},'showfig',false,'cellclass',[]),...
        spikes,'UniformOutput',false);
    PeriEventISIDist.cells = cat(1,temp{:});
    
    if ~isempty(cellclass)
        %Check for empty cell class entries
        noclass = cellfun(@isempty,cellclass);
        classnames = unique(cellclass(~noclass));
        numclasses = length(classnames);
        for cl = 1:numclasses
            inclasscells{cl} = strcmp(classnames{cl},cellclass);
            PeriEventISIDist.pop.(classnames{cl}) = bz_CollapseStruct(PeriEventISIDist.cells(inclasscells{cl}),...
                3,'mean',true);
        end
    end
    
    return
end


%%
ISIs = diff(spikes);
spiketimes = spikes(2:end-1);
ISInp1 = ISIs(2:end);
ISIs = ISIs(1:end-2);

ISIs_rel = [];
ISInp1_rel = [];
spiketimes_rel = [];
numevents =length(eventtimes);
for ee = 1:numevents
    reltime = spiketimes-eventtimes(ee);
    inwin = reltime >= winsize(1) & reltime <= winsize(2);
    
    switch whichISIs
        case 'both'
            ISIs_rel = [ISIs_rel;ISIs(inwin)];
            ISInp1_rel = [ISInp1_rel;ISInp1(inwin)];
            spiketimes_rel = [spiketimes_rel;reltime(inwin)];
        case 'prev'
            ISIs_rel = [ISIs_rel;ISIs(inwin)];
            ISInp1_rel = [ISInp1_rel;ISIs(inwin)];
            spiketimes_rel = [spiketimes_rel;reltime(inwin)]; 
        case 'next'
            ISIs_rel = [ISIs_rel;ISInp1(inwin)];
            ISInp1_rel = [ISInp1_rel;ISInp1(inwin)];
            spiketimes_rel = [spiketimes_rel;reltime(inwin)];
    end
end

%%
[ PeriEventISIDist ] = ConditionalHist( [spiketimes_rel;spiketimes_rel],log10([ISIs_rel;ISInp1_rel]),...
    'Xbounds',winsize,'numXbins',numXbins,'Ybounds',logISIbounds,'numYbins',numISIbins,'minX',minX);

PeriEventISIDist.rate = PeriEventISIDist.Xhist./(diff(PeriEventISIDist.Xbins([1 2])).*numevents*2); 
%factor of 2 for double counting n and n+1
%%
if SHOWFIG
    figure
        imagesc(PeriEventISIDist.Xbins,(PeriEventISIDist.Ybins),PeriEventISIDist.pYX')
        hold on
        plot(PeriEventISIDist.Xbins,log10(1./PeriEventISIDist.rate),'r')
        LogScale('y',10)
end



end

