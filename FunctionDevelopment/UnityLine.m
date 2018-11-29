function [  ] = UnityLine(  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

xrange = get(gca,'xlim');
yrange = get(gca,'ylim');

plot([min([xrange yrange]) max([xrange yrange])],...
    [min([xrange yrange]) max([xrange yrange])],'k:')

end

