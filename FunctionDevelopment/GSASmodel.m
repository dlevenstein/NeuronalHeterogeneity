function [allISIdist] = GSASmodel(GSASparms,logtbins,numcells,numAS)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%
%Note: vectr form uses logCV... for fitting
%logtbins should be base e

if ~isstruct(GSASparms)
    GSASparms = convertGSASparms(GSASparms,numcells,numAS);
end

%Here: if logtbins = 'sample'. Put in a very large vector of possible times
%to sample from, save that we're sampling.
if strcmp(logtbins,'sample')
    logtbins = 
end

GSISI = LogGamma(GSASparms.GSlogrates,GSASparms.GSCVs,GSASparms.GSweights,logtbins');
%%
for aa = 1:length(GSASparms.ASlogrates)
    ASISI(:,:,aa) = LogGamma(GSASparms.ASlogrates(aa),GSASparms.ASCVs(aa),GSASparms.ASweights(:,aa)',logtbins');
end
%%
allISIdist = sum(ASISI,3)+GSISI;

%Here: multiply allISIdist by some (large) number to get counts and return
%a sample with those counts... (for KS test). use CUMsum
end

