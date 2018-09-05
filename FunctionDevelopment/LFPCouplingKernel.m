function [ lambda ] = LFPCouplingKernel( binnedpowers,phases,prefphase,powerdependence )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%phasepowergram = spkmat.specgram;
%% Nonparametric
% npowerbins = 15;
% nfreqs = size(phasepowergram,2);
% 
% %Get the phase and normalized power from the spectrogram
% abX = NormToInt(log10(abs(phasepowergram)),'modZ');
% angX = angle(phasepowergram);
% 
% %Bin the powers
% poweredges = linspace(-1.75,1.75,npowerbins+1);
% powercenters = poweredges(1:end-1)+0.5.*diff(poweredges(1:2));
% poweredges(1) = -Inf; poweredges(end) = Inf;
% [~,~,powerBIN] = histcounts(abX,poweredges);
% 
% binnedpowers = zeros(size(powerBIN,1),npowerbins);
% for ff = 1:nfreqs
%     for tt = 1:length(powerBIN)
%     %binnedphasepowers(tt,phaseBIN(tt),powerBIN(tt)) = 1;
%     binnedpowers(tt,powerBIN(tt,ff),ff) = 1;
%     end
% end

%%
%keep = powerdependence;
%nfreqs = size(phases,2);
%%
%powerdependence = reshape(powerdependence,1,[],nfreqs);
%powerdependence = permute(powerdependence,[1,3,2]);
%%
ntsteps = size(phases,1);
TH = phases;
POW = reshape(binnedpowers,[ntsteps,length(powerdependence)]);

%A = @(a,POW,TH) (POW*a(2:npowerbins+1)'.*cos(TH+a(npowerbins+2)) + log(a(1)));
%lambda = zeros(
lambda = POW*powerdependence'.*sum(cos(bsxfun(@plus,TH,prefphase)),2);

end

