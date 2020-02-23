function [allISIdist] = GSASmodel(GSASparms,logtbins,numcells,numAS)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%
if ~isstruct(GSASparms)
    GSASparms = convertGSASparms(GSASparms,numcells,numAS);
end

GSISI = LogGamma(GSASparms.GSlogrates,GSASparms.GSCVs,GSASparms.GSweights,logtbins');
%%
for aa = 1:length(GSASparms.ASlogrates)
    ASISI(:,:,aa) = LogGamma(GSASparms.ASlogrates(aa),GSASparms.ASCVs(aa),GSASparms.ASweights(:,aa)',logtbins');
end
%%
allISIdist = sum(ASISI,3)+GSISI;

end

