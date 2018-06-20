function [ C ] = ConnectSpatialNetwork(N,p_mean,width)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%
%N  number of neurons, can be [Nin Nout] if connecting two different
%   populations
%Size spatial extent of the network (periodic) default: 1
%
%C output connectivity matrix = 1 indicates connection, 0 no connection
%% Dev

N = 5000;
p_mean = 0.1;
Size = 2.*pi;
width = 1;
%%
locations = Size.*rand(N,2);
[~,sortlocations] = sort(locations(:,1));
locations = locations(sortlocations,:);

xdistances = pdist(locations(:,1),'cityblock');
xdistances = squareform(xdistances);

ydistances = pdist(locations(:,2),'cityblock');
ydistances = squareform(ydistances);

distances = pdist(locations);
distances = squareform(distances);
%%

k = 1/(width.^2);
gx = exp(k.*cos(xdistances))./(2.*pi*besselj(0,k));
gy = exp(k.*cos(ydistances))./(2.*pi*besselj(0,k));

pconnect = gx.*gy;
pconnect = pconnect.* (p_mean./mean(pconnect(:))); %Ad hoc normalization...

connect = rand(size(pconnect))<pconnect;
%%
figure
imagesc(distances)
%%



excell = randsample(N,1);
th = 0:pi/500:2*pi;
xunit = width * cos(th) + locations(excell,1);
yunit = width * sin(th) + locations(excell,2);

figure
plot(locations(:,1),locations(:,2),'k.','markersize',1)
hold on
plot(locations(connect(:,excell),1),locations(connect(:,excell),2),'r.','markersize',10)
%scatter(locations(:,1),locations(:,2),2,pconnect(:,excell))
plot(mod(xunit,2.*pi), mod(yunit,2.*pi),'k.');
%colorbar
axis tight

%%
xtest = linspace(-10,10,500);
figure
plot(xtest,exp(k.*cos(xtest).*(2*pi))./(2.*pi*besselj(0,k)))
end

