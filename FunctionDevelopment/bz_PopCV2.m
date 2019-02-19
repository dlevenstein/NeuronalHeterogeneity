function [ PopCV2 ] = bz_PopCV2( input_args )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% If Exists, load


%%
CV2mat.winsize = spkwinsize;
CV2mat.timestamps = spikemat.timestamps;
CV2mat.binedges = bsxfun(@(X,Y) X+Y,spikemat.timestamps,[-0.5 0.5].*CV2mat.winsize);
for tt = 1:length(cellclasses)
    allspikes.CV2.(cellclasses{tt}) = cat(1,ISIStats.allspikes.CV2{CellClass.(cellclasses{tt})});
    allspikes.times.(cellclasses{tt}) = cat(1,ISIStats.allspikes.times{CellClass.(cellclasses{tt})});
    [CV2mat.timestamps,CV2mat.(cellclasses{tt})] = ...
        BinDataTimes(allspikes.CV2.(cellclasses{tt}),allspikes.times.(cellclasses{tt}),CV2mat.binedges);
    CV2mat.rate.(cellclasses{tt}) = interp1(spikemat.timestamps,spikemat.poprate.(cellclasses{tt}),CV2mat.timestamps);
end

%% Save


end

