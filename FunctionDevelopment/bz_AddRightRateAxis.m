function [] = bz_AddRightRateAxis()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    bounds = ylim(gca);
    yyaxis right
    ylim(-fliplr(bounds))
    LogScale('y',10,'nohalf',true,'exp',true)
    ylabel('Rate (Hz)')
    
end

