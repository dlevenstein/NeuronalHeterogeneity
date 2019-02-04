function [ meanZ,N ] = ConditionalHist3( X,Y,Z,varargin )
%[ jointXYZ ] = ConditionalHist3( X,Y,Z ) for a set of observations [X,Y,Z] this 
%function calculates the statistics of Z given X and Y.
%%
p = inputParser;
addParameter(p,'numXbins',50)
addParameter(p,'numYbins',50)
addParameter(p,'Xbounds',[])
addParameter(p,'Ybounds',[])
addParameter(p,'minXY',25)
addParameter(p,'sigma',[])
parse(p,varargin{:})
numXbins = p.Results.numXbins;
numYbins = p.Results.numYbins;
Xbounds = p.Results.Xbounds;
Ybounds = p.Results.Ybounds;
minXY = p.Results.minXY;
sig = p.Results.sigma;

if ~isempty(sig)
    bintype = 'gaussian';
else
    bintype = 'bins';
end

%%
if isempty(Xbounds)
    Xbounds(1) = min(X); Xbounds(2) = max(X);
end
if isempty(Ybounds)
    Ybounds(1) = min(Y(~isinf(Y))); Ybounds(2) = max(Y(~isinf(Y)));
end

Xedges = linspace(Xbounds(1),Xbounds(2),numXbins+1);
Xbins = Xedges(1:end-1)+ 0.5.*diff(Xedges([1 2]));
Xedges(1) = -inf;Xedges(end) = inf;

Yedges = linspace(Ybounds(1),Ybounds(2),numYbins+1);
Ybins = Yedges(1:end-1)+ 0.5.*diff(Yedges([1 2]));
Yedges(1) = -inf;Yedges(end) = inf;


%%
[N,~,~,BINX,BINY] = histcounts2(X,Y,Xedges,Yedges);
meanZ = zeros(numXbins,numYbins);


        
switch bintype
    case 'bins'
        for xx = 1:length(Xbins)
            for yy = 1:length(Ybins)
                meanZ(xx,yy) = nanmean(Z(BINX==xx & BINY==yy));
            end
        end

    case 'gaussian'
        for xx = 1:length(Xbins)
            for yy = 1:length(Ybins)
                pointdist = sqrt((X-Xbins(xx)).^2 + (Y-Ybins(yy)).^2);
                weight = exp(-.5 * (pointdist/sig) .^ 2) ./ (sig * sqrt(2*pi));   %Weight by gaussian

                N(xx,yy) = sum(weight);
                meanZ(xx,yy) = sum(Z.*weight)./N(xx,yy);
            end
        end
end

meanZ(N<minXY)=nan;
N = N./length(Z);
%weightmean(totweight<minXY)=nan;

%%
% figure
% subplot(2,2,1)
% imagesc(Xbins,Ybins,weightmean')
% hold on
% %plot(X,Y,'k.')
% 
% subplot(2,2,2)
% imagesc(Xbins,Ybins,totweight')
% hold on
% plot(X,Y,'k.')
% colorbar
% caxis([0 50])
% subplot(2,2,1)
% plot(Xbins(xx),Ybins(yy),'+')
% hold on
% scatter(X,Y,1,pointdist)
% subplot(2,2,2)
% %plot(Xbins(xx),Ybins(yy),'+')
% hold on
% scatter(X,Y,1,weight)
% colorbar
%%
% figure
% subplot(2,2,1)
% plot(X,Y,'.')
% subplot(2,2,2)
% imagesc(Xbins,Ybins,N)
% axis xy
% 
% subplot(2,2,3)
% imagesc(Xbins,Ybins,meanZ)
% colorbar
% axis xy

