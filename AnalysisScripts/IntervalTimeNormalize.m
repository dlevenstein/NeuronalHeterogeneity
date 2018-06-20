function [ normt ] = IntervalTimeNormalize( t,interval,outintscaling )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Normalizes time in the interval to 0-1. Before the interval time will be
%negative. After the interval time will be +1. Outside interval is scaled
%by outintscaling factor (i.e. if outintscaling=4, -1 will be 4s before
%interval start
% <-scale------S---E-------+scale>
%      -1      0   1      +1
%%

trelstart = t-interval(1);
trelend = t-interval(2);

normt = zeros(size(t));
normt(trelstart<0) = trelstart(trelstart<0)./outintscaling;
normt(trelend>0) = trelend(trelend>0)./outintscaling+1;

inintt = trelstart>=0 & trelend<=0;

normt(inintt) = interp1(interval,[0 1],t(inintt));

end

