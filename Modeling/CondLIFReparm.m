function [ Re,Ri ] = CondLIFReparm( V_inf,Gamma,cellparams,synparams )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Re = cellparams.g_L.*(V_inf.*Gamma - cellparams.E_L - (Gamma-1).*cellparams.E_i) ./ ...
        synparams.w_e.*synparams.tau_se.*(cellparams.E_e-cellparams.E_i);
    
Ri = (cellparams.g_L.*(Gamma-1) - Re.*synparams.w_e.*synparams.tau_se) ./ ...
    (synparams.w_i.*synparams.tau_si);

Re = Re./1000;
Ri = Ri./1000;

end

