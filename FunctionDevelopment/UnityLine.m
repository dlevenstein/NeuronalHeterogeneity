function [  ] = UnityLine( varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%
p = inputParser;
addParameter(p,'linetype',':')
addParameter(p,'linecolor','k')
parse(p,varargin{:})
linetype = p.Results.linetype;
linecolor = p.Results.linecolor;

%%
xrange = get(gca,'xlim');
yrange = get(gca,'ylim');

plot([min([xrange yrange]) max([xrange yrange])],...
    [min([xrange yrange]) max([xrange yrange])],linetype,'color',linecolor)

end

