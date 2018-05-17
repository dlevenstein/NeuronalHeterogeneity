function [ P_T ] = PoissonBasis(r,P_r,T)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%NEXT:gamma distribution basis with constant t_ref?
%% DEV
r = 0.01:0.1:20;
P_r = LogNorm(r,-2,1);
P_r = P_r./sum(P_r);
T = 0:0.001:10;

figure
plot(r,P_r)


%%
[gridT,gridR] = meshgrid(T,r);
gridP = gridR.*exp(-gridT.*gridR);
P_T = sum(gridP,1);
P_T = P_T./sum(P_T);
%%
figure
    subplot(2,2,1)
        imagesc(T,r,gridP)
        colorbar
        axis xy
        xlabel('T');ylabel('R')
    subplot(2,2,2)
        plot(T,gridP)
    subplot(2,2,3)
        plot(T,P_T,'k')

end

