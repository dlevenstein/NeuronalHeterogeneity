function [ R_e,R_i ] = CondLIFReparm( V_inf,Gamma,cellparams,synparams )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


g_e = cellparams.g_L.*(V_inf.*Gamma - cellparams.E_L - (Gamma-1).*cellparams.E_i) ./ ...
        (cellparams.E_e-cellparams.E_i);

g_i =  cellparams.g_L.*(Gamma-1)-g_e;

R_e = g_e./(synparams.w_e.*synparams.tau_se.*synparams.K_e);
R_i = g_i./(synparams.w_i.*synparams.tau_si.*synparams.K_i);
R_e = R_e.*1000;
R_i = R_i.*1000;


end