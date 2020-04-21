function [cdfX,cdfY] = PlotCDF(data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% 
cdfX = sort(data);
cdfY = [1:length(data)]./length(data);

end

