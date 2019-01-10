function [ CONDXY ] = ConditionalHist( X,Y,varargin)
%[CONDXY] = ConditionalHist(X,Y) Calculates the conditional probabilty 
%of Y given X
%
%
%OUTPUT
%   CONDXY
%       .pYX    a [numXbins x numYbins] matrix in which the [x,y]th element is P(y|x) 
%       .XYhist joint histogram of X and Y
%       .Xhist  histogram of X
%       .Xbins  X bins
%       .Ybins  Y bins
%
%DLevenstein 2019
%% DEV
% X = CV2mat.PSS;
% Y = CV2mat.pE;

%%
p = inputParser;
addParameter(p,'numXbins',50)
addParameter(p,'numYbins',50)
addParameter(p,'Xbounds',[])
addParameter(p,'Ybounds',[])
addParameter(p,'minX',25)
parse(p,varargin{:})
numXbins = p.Results.numXbins;
numYbins = p.Results.numYbins;
Xbounds = p.Results.Xbounds;
Ybounds = p.Results.Ybounds;
minX = p.Results.minX;

%%
if isempty(Xbounds)
    Xbounds(1) = min(X); Xbounds(2) = max(X);
end
if isempty(Ybounds)
    Ybounds(1) = min(Y); Ybounds(2) = max(Y);
end

Xbins = linspace(Xbounds(1),Xbounds(2),numXbins+1);
Xbins = Xbins(1:end-1)+ 0.5.*diff(Xbins([1 2]));

Ybins = linspace(Ybounds(1),Ybounds(2),numYbins+1);
Ybins = Ybins(1:end-1)+ 0.5.*diff(Ybins([1 2]));

%First calculate the marginal probability of X
Xhist = hist(X,Xbins);
Xhist4norm = Xhist;Xhist4norm(Xhist4norm<=minX) = nan;

%Then calculate the joint probabilty of X and Y
[XYhist] = hist3([X,Y],{Xbins,Ybins});

%
pYX = bsxfun(@(x,y) x./y,XYhist,Xhist4norm');

CONDXY.pYX = pYX;
CONDXY.XYhist = XYhist;
CONDXY.Xhist = Xhist;
CONDXY.Xbins = Xbins;
CONDXY.Ybins = Ybins;


%%
% figure
% imagesc(CONDXY.Xbins,CONDXY.Ybins,CONDXY.pYX')
% axis xy
end

