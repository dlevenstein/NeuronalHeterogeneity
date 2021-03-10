function [pS_x,pAS_x,pS_AS,pS_GS] = EncodingModal_VMPoisson(s,x,x0,k,pAS_0,pAS_pi,rAS,rGS,dt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%
%Tuning Curve: P(AS|x)
%pAS_x = exp(k.*cos(x-x0))./(2.*pi.*besseli(0,k)); (Von misses)
%Note: normalization removes bessel function!
pAS_x = (exp(cos(x-x0)./k)-exp(-1./k))./(exp(1./k)-exp(-1./k));
pAS_x = pAS_x.*(pAS_0-pAS_pi)+pAS_pi;

%AS Likelihood Function
pS_AS = (exp(-rAS.*dt).*(rAS*dt).^s)./factorial(s);

%GS Likelihood Function
pS_GS = (exp(-rGS.*dt).*(rGS*dt).^s)./factorial(s);

%Total likelihood
pS_x = pAS_x.*(pS_AS - pS_GS) + pS_GS;

end

