function [ logpdf,logbinedges,logbincenters ] = LinToLogPDF( linpdf,linbinedges,numlogbins )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%   Make option to do hist as well...
%% DEV
linpdf = ccg.hists(:,1,1);
linbinedges = [ccg.t(1)-0.5.*diff(ccg.t([1 2])); ccg.t+0.5.*diff(ccg.t([1 2]))];
numlogbins

end

