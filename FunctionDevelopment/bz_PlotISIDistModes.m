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
addParameter(p,'whichShare',1)
addParameter(p,'dotscale',200)
addParameter(p,'dotscaleAS',[])
parse(p,varargin{:})
sharORsing = p.Results.sharORsing;
ws = p.Results.whichShare;
dotscale = p.Results.dotscale;
dotscaleAS = p.Results.dotscaleAS;
if isempty(dotscaleAS)
    dotscaleAS = dotscale;
end
% add: Mode color...
%%
% Find the cell that matches the UID
plotcell = find(GammaFits.cellstats.UID==UID);
%Option: random
%plotcell = randi(GammaFits.numcells,1);
%UID
ASdots = 'filled';
switch sharORsing
    case 'sharedfit'

            GFmodel.ASlogrates = GammaFits.sharedfit(ws).ASlogrates;
            GFmodel.ASCVs = GammaFits.sharedfit(ws).ASCVs;
            GFmodel.ASweights = GammaFits.sharedfit(ws).ASweights(plotcell,:);
            GFmodel.GSlogrates = GammaFits.sharedfit(ws).GSlogrates(plotcell);
            GFmodel.GSCVs = GammaFits.sharedfit(ws).GSCVs(plotcell);
            GFmodel.GSweights = GammaFits.sharedfit(ws).GSweights(plotcell);
            
        if length(plotcell)>1
            %Percentiles for Range Bars
            ASstds = prctile(GammaFits.sharedfit(ws).ASweights(plotcell,:),[20 80],1);
            GSstds = prctile(GammaFits.sharedfit(ws).GSweights(plotcell),[20 80],2);
            GSstds_rate = prctile(GammaFits.sharedfit(ws).GSlogrates(plotcell),[20 80],2);
            ASdots = 'o';
        elseif length(plotcell)< 1
            error('No Cells with that UID')
        end
    case 'singlecell'
        GFmodel = GammaFits.singlecell(plotcell);
end

numAS = length(GFmodel.ASlogrates);
fitcolor = 'k';
GScolor = [0.6 0.4 0];

scaleDist = 3.5;
offset = 0.7;
%ISI Distribution
plot(GammaFits.logtimebins,...
    mean(GammaFits.ISIdists(:,plotcell),2).*scaleDist+offset,...
    'color',[0.5 0.5 0.5],'linewidth',2)
hold on
%Full Model
plot(GammaFits.logtimebins,...
    mean(GSASmodel(GFmodel,...
    GammaFits.taubins),2).*scaleDist+offset,...
    fitcolor,'linewidth',1)

%GS and AS modes
plot(GammaFits.logtimebins,...
    LogGamma(mean(GFmodel.GSlogrates),...
    mean(GFmodel.GSCVs),...
    mean(GFmodel.GSweights)',...
    GammaFits.taubins').*scaleDist+offset,'color',GScolor,'linewidth',0.25);

for aa = 1:numAS
    plot(GammaFits.logtimebins,...
        LogGamma(GFmodel.ASlogrates(aa),...
        GFmodel.ASCVs(aa),...
        mean(GFmodel.ASweights(:,aa),1)',...
        GammaFits.taubins').*scaleDist+offset,'k','linewidth',0.25);
    
    %For multiple cells: bounds
    if length(plotcell)>1
        plot(-GFmodel.ASlogrates(aa).*[1 1],...
            LogGamma(GFmodel.ASlogrates(aa),...
            GFmodel.ASCVs(aa),...
            ASstds(:,aa)',...
            -GFmodel.ASlogrates(aa)./log10(exp(1))').*scaleDist+offset,'-k','linewidth',0.25);
    end
end
    
    %For Multiple cells: bounds
    if length(plotcell)>1
        plot(-mean(GFmodel.GSlogrates).*[1 1],...
            LogGamma(mean(GFmodel.GSlogrates),...
            mean(GFmodel.GSCVs),...
            GSstds(:)',...
            -mean(GFmodel.GSlogrates)./log10(exp(1))').*scaleDist+offset,'-','color',GScolor,'linewidth',0.25);
        
        plot(-GSstds_rate,...
            LogGamma(mean(GFmodel.GSlogrates),...
            mean(GFmodel.GSCVs),...
            mean(GFmodel.GSweights)',...
            -[1 1].*mean(GFmodel.GSlogrates)'./log10(exp(1))').*scaleDist+offset,'-','color',GScolor,'linewidth',0.25);
    end
box on
axis tight

    ylabel(['UID: ',num2str(GammaFits.cellstats.UID(plotcell))])

%xlim([-3 2])


scatter(-GFmodel.ASlogrates,...
    log10(GFmodel.ASCVs),...
    dotscaleAS*mean(GFmodel.ASweights,1)+0.00001,[0 0 0],ASdots)
hold on
scatter(-GFmodel.GSlogrates,...
    log10(GFmodel.GSCVs),...
    dotscale*GFmodel.GSweights+0.00001,GScolor,'filled')
plot(GammaFits.logtimebins([1 end]),[0 0],'k--')
ylabel('CV');xlabel('mean ISI (s)')
xlim([-2.75 1.75])
%ylim([-2 0.75])
LogScale('x',10,'exp',true,'nohalf',true)
LogScale('y',10)
%box off





end

