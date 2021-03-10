function [model_m,model_c,tuningcurve] = FitEncodingModel_HD(s,x,dt)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Possible here: make a guess that's close! (Use continuous model fit?....)

maxspikes = 15;
    CONDXY = ConditionalHist(x,s,...
         'numXbins',20,'Xbounds',[0 2*pi],'numYbins',maxspikes+1,'Ybounds',[0 maxspikes],...
         'Xbinoverlap',2);
     
[peakrate,peakx]=max(CONDXY.meanYX);
peakrate = peakrate./dt;
peakx = CONDXY.Xbins(peakx);
%%

init = zeros(6,1);
init(1) = peakx;    %x0
init(2) = 1;    %k
init(3) = 0.8;  %pAS_0
init(4) = 0.1;  %pAS_pi
init(5) = peakrate;   %rAS
init(6) = 0.5;    %rGS

lb(1) = -5*pi;    %x0
lb(2) = 0.1;    %k
lb(3) = 0;   %pAS_0
lb(4) = 0;  %pAS_pi
lb(5) = 1;   %rAS
lb(6) = 1e-2;    %rGS

ub(1) = 5*pi;    %x0
ub(2) = 10;    %k
ub(3) = 1;   %pAS_0
ub(4) = 1;  %pAS_pi
ub(5) = 200;   %rAS
ub(6) = 20;    %rGS


likelihood_c = @(a,x,s) (EncodingContinuous_VMPoisson(double(s),x,a(1),a(2),a(5),a(6),dt));
nlogL_c = @(k) -(nansum(log(likelihood_c(k,x,s))));

likelihood_m = @(a,x,s) (EncodingModal_VMPoisson(double(s),x,a(1),a(2),a(3),a(4),a(5),a(6),dt));
nlogL_m = @(k) -(nansum(log(likelihood_m(k,x,s))));
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');%,'UseParallel',true);
%try also: 'Algorithm','active-set','sqp'
options.MaxFunctionEvaluations = 2e4;

%%
kernelPredict_c = fmincon(nlogL_c,init,[],[],[],[],lb,ub,[],options);

init(5)=kernelPredict_c(5); init(1)=kernelPredict_c(1);
kernelPredict_m = fmincon(nlogL_m,init,[],[],[],[],lb,ub,[],options);

%%
nll_m = nlogL_m(kernelPredict_m);
nll_c = nlogL_c(kernelPredict_c);

AIC_m = 2.*nll_m + 2.*6;
AIC_c= 2.*nll_c + 2.*4;

BIC_m = 2.*nll_m + log(sum(~isnan(x))).*6;
BIC_c= 2.*nll_c + log(sum(~isnan(x))).*4;
%%
%likelihood(init,x,s);
model_m.nll = nll_m;
model_m.AIC = AIC_m;
model_m.BIC = BIC_m;
model_m.parms.x0 = kernelPredict_m(1);
model_m.parms.k=kernelPredict_m(2);
model_m.parms.pAS_0 = kernelPredict_m(3);
model_m.parms.pAS_pi = kernelPredict_m(4);
model_m.parms.rAS = kernelPredict_m(5);
model_m.parms.rGS = kernelPredict_m(6);

model_c.nll = nll_c;
model_c.AIC = AIC_c;
model_c.BIC = BIC_c;
model_c.parms.x0 = kernelPredict_c(1);
model_c.parms.k=kernelPredict_c(2);
model_c.parms.R_0 = kernelPredict_c(5);
model_c.parms.R_pi = kernelPredict_c(6);

tuningcurve.Xbins = CONDXY.Xbins
tuningcurve.meanRate = CONDXY.meanYX./dt;
%%

%Check out best nlogL vs one you pick by eye...

x_plot = linspace(0,2.*pi,100);
s_plot = 0:maxspikes;
[X,S] = meshgrid(x_plot,s_plot);

x0=kernelPredict_m(1);
k=kernelPredict_m(2);
pAS_0 = kernelPredict_m(3);
pAS_pi = kernelPredict_m(4);

%dt = 0.5;
rAS = kernelPredict_m(5);
rGS = kernelPredict_m(6);


[pS_x_modal,pAS_x,pS_AS,pS_GS] = EncodingModal_VMPoisson(S,X,model_m.parms.x0,...
    model_m.parms.k,model_m.parms.pAS_0,model_m.parms.pAS_pi,model_m.parms.rAS,...
    model_m.parms.rGS,dt);
x0=kernelPredict_c(1);
k=kernelPredict_c(2);
R_0 = kernelPredict_c(5);
R_pi = kernelPredict_c(6);
[pS_x_cont,R_x] = EncodingContinuous_VMPoisson(S,X,x0,k,R_0,R_pi,dt);
%% Plot: fitted vs observed tuning curve

%Tuning Curves
%AS/GS observation distirbution
%Comapre to ISI dist and fit!


figure

    subplot(3,3,1)
        imagesc(CONDXY.Xbins,CONDXY.Ybins,CONDXY.pYX')
        hold on
        %imagesc(CONDXY(1).Xbins+2.*pi,CONDXY(1).Ybins,CONDXY(HDcells(cc)).pYX')
        axis xy
        %xlim([0 4*pi])
        bz_piTickLabel('x')
        ColorbarWithAxis([0 0.4],'P[s|HD]')
        xlabel('Head Direction');ylabel('Spike Count')
        title('Observed')

    subplot(3,3,4)
        plot(AIC_c,AIC_m,'.')
        hold on
        UnityLine
    
    subplot(3,3,2)
        imagesc(x_plot,s_plot,pS_x_modal)
        xlabel('x');ylabel('s')
        ColorbarWithAxis([0 0.5] ,'P[s|x]')
        axis xy
        title('Probabilistic Modal Tuning')
        bz_piTickLabel('x')
        
        subplot(3,3,3)
            imagesc(x_plot,s_plot,pS_x_cont)
            xlabel('x');ylabel('s')
            ColorbarWithAxis([0 0.5] ,'P[s|x]')
            axis xy
            title('Continuous Rate Tuning')
            bz_piTickLabel('x')
            
    subplot(3,3,5)
        plot(s_plot,pS_AS(:,1))
        hold on
        plot(s_plot,pS_GS(:,1))
        xlabel(['# Spikes in dt=',num2str(dt),' bin (s)'])
        ylabel('P(s|state)')
        box off
        legend(['AS: ',num2str(rAS),'Hz'],['GS: ',num2str(rGS),'Hz'])
        
subplot(3,3,6)
    %plot(x,pAS_x(1,:))
    plot(x_plot,pAS_x(1,:))
    xlabel('x');ylabel('P[AS|x]') 
    ylim([0 1])
    %xlim([-pi pi]+peakx)
    xlim([0 2.*pi])
    bz_piTickLabel('x')
    yyaxis right
    plot (x_plot,R_x(1,:),'--')
    ylim([rGS rAS])
    ylabel('R[x]') 
end

