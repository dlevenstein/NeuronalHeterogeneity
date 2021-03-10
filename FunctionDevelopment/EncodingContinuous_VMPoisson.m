function [pS_x,R_x] = EncodingContinuous_VMPoisson(s,x,x0,k,R_0,R_pi,dt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%
%Tuning Curve: P(AS|x)
%pAS_x = exp(k.*cos(x-x0))./(2.*pi.*besseli(0,k)); (Von misses)
%Note: normalization removes bessel function!
R_x = (exp(cos(x-x0)./k)-exp(-1./k))./(exp(1./k)-exp(-1./k));
R_x = R_x.*(R_0-R_pi)+R_pi;

%AS Likelihood Function
pS_R = (exp(-R_x.*dt).*(R_x*dt).^s)./factorial(s);

%Total likelihood
pS_x = pS_R;

%Add: option to show full X,S space and curves...
end

