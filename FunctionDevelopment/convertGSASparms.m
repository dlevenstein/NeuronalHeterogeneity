function [parmstype2] = convertGSASparms(parmstype1,numcells,numAS)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% If structure, returns vector (for use in fitting algorithms).
% If vector, returns structure, (needs numcells, numAS)

if isstruct(parmstype1)
    parmstype2 = [parmstype1.GSlogrates'; parmstype1.GSCVs'; parmstype1.GSweights';,...
        parmstype1.ASlogrates'; parmstype1.ASCVs'; parmstype1.ASweights(:)];
else
    parmstype2.GSlogrates = parmstype1(1:numcells)';
    parmstype2.GSCVs = parmstype1(numcells+1 : 2*numcells)';
    parmstype2.GSweights = parmstype1(2*numcells+1 : 3*numcells)';
    
    parmstype2.ASlogrates = parmstype1(3*numcells+1 : 3*numcells+numAS)';
    parmstype2.ASCVs = parmstype1(3*numcells+numAS+1 : 3*numcells+2*numAS)';
    parmstype2.ASweights = parmstype1(3*numcells+2*numAS+1 : end)';
    parmstype2.ASweights = reshape(parmstype2.ASweights,numcells,numAS);
end

end

