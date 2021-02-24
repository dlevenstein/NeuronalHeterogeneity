function [] = bz_PlotISIDistModes(GammaFits,UID,varargin)
% Function for plotting mixture of gamma modes for the ISI distribtion
%Does best with subplot(2,3,x)
%
%INPUT
%   GammaFits   mixture of gamma model fit from bz_FitISISharedGammaModes
%               via SharedGammaModeFitAnalysis.m. Should have the following
%               fields:
%                   .logtimebins .ISIdists .taubins .cellstats
%                   .sharedfit or .singlecell
%   UID         which cell to plot?

%Could also load from GammaFit = bz_LoadCellinfo(basePath,'GammaFit');
%%
% parse args
p = inputParser;
addParameter(p,'sharORsing','sharedfit')
parse(p,varargin{:})
sharORsing = p.Results.sharORsing;

% add: Mode color...
%%
% Find the cell that matches the UID
plotcell = find(GammaFits.cellstats.UID==UID);
%Option: random
%plotcell = randi(GammaFits.numcells,1);
%UID
switch sharORsing
    case 'sharedfit'
        GFmodel.ASlogrates = GammaFits.sharedfit.ASlogrates;
        GFmodel.ASCVs = GammaFits.sharedfit.ASCVs;
        GFmodel.ASweights = GammaFits.sharedfit.ASweights(plotcell,:);
        GFmodel.GSlogrates = GammaFits.sharedfit.GSlogrates(plotcell);
        GFmodel.GSCVs = GammaFits.sharedfit.GSCVs(plotcell);
        GFmodel.GSweights = GammaFits.sharedfit.GSweights(plotcell);
    case 'singlecell'
        GFmodel = GammaFits.singlecell(plotcell);
end

numAS = length(GFmodel.ASlogrates);
fitcolor = 'k';
GScolor = [0.6 0.4 0];

scaleDist = 3.5;
offset = 0.7;
dotscale = 200;
%ISI Distribution
plot(GammaFits.logtimebins,...
    GammaFits.ISIdists(:,plotcell).*scaleDist+offset,...
    'color',[0.5 0.5 0.5],'linewidth',2)
hold on
%Full Model
plot(GammaFits.logtimebins,...
    GSASmodel(GFmodel,...
    GammaFits.taubins).*scaleDist+offset,...
    fitcolor,'linewidth',1)

%GS and AS modes
plot(GammaFits.logtimebins,...
    LogGamma(GFmodel.GSlogrates,...
    GFmodel.GSCVs,...
    GFmodel.GSweights',...
    GammaFits.taubins').*scaleDist+offset,'color',GScolor,'linewidth',0.25);
for aa = 1:numAS
    plot(GammaFits.logtimebins,...
        LogGamma(GFmodel.ASlogrates(aa),...
        GFmodel.ASCVs(aa),...
        GFmodel.ASweights(aa)',...
        GammaFits.taubins').*scaleDist+offset,'k','linewidth',0.25);
end
%box off
axis tight

    ylabel(['UID: ',num2str(GammaFits.cellstats.UID(plotcell))])

xlim([-3 2])


scatter(-GFmodel.ASlogrates(:),...
    log10(GFmodel.ASCVs(:)),...
    dotscale*GFmodel.ASweights(:)+0.00001,'k','filled')
hold on
scatter(-GFmodel.GSlogrates,...
    log10(GFmodel.GSCVs),...
    dotscale*GFmodel.GSweights+0.00001,GScolor,'filled')
plot(GammaFits.logtimebins([1 end]),[0 0],'k--')
ylabel('CV');xlabel('mean ISI (s)')
xlim([-2.75 1.75])
% ylim([-2 0.75])
LogScale('x',10,'exp',true,'nohalf',true)
LogScale('y',10)
box off





end

