%CompareEncodingModels.m

figfolder = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/Modeling/Figures/EncodingModels';
%%

x = linspace(-pi,pi,100);
s = 0:15;
[X,S] = meshgrid(x,s);

x0=0;
k=1;
pAS_0 = 0.7;
pAS_pi = 0.05;

dt = 0.5;
rAS = 10;
rGS = 1;


[pS_x_modal,pAS_x,pS_AS,pS_GS] = EncodingModal_VMPoisson(S,X,x0,k,pAS_0,pAS_pi,rAS,rGS,dt);

%[pS_x_modal,pAS_x,pS_AS,pS_GS] = EncodingContinuous_VMPoisson(S,X,x0,k,pAS_0,pAS_pi,rAS,rGS,dt);
R_0 = pAS_0.*rAS + (1-pAS_0).*rGS;
R_pi = pAS_pi.*rAS + (1-pAS_pi).*rGS;
[pS_x_cont,R_x] = EncodingContinuous_VMPoisson(S,X,x0,k,R_0,R_pi,dt);


%%


figure
subplot(2,2,1)
    %plot(x,pAS_x(1,:))
    plot(x,pAS_x(1,:))
    xlabel('x');ylabel('P[AS|x]') 
    ylim([0 1])
    xlim([-pi pi])
    bz_piTickLabel('x')
    yyaxis right
    plot (x,R_x(1,:),'--')
    ylim([rGS rAS])
    ylabel('R[x]') 
    
subplot(3,2,2)
    plot(s,pS_AS(:,1))
    hold on
    plot(s,pS_GS(:,1))
    xlabel(['# Spikes in dt=',num2str(dt),' bin (s)'])
    ylabel('P(s|state)')
    legend(['AS: ',num2str(rAS),'Hz'],['GS: ',num2str(rGS),'Hz'])

subplot(2,2,3)
imagesc(x,s,pS_x_modal)
xlabel('x');ylabel('s')
ColorbarWithAxis([0 0.5] ,'P[s|x]')
axis xy
title('Probabilistic Modal Tuning')
bz_piTickLabel('x')

subplot(2,2,4)
imagesc(x,s,pS_x_cont)
xlabel('x');ylabel('s')
ColorbarWithAxis([0 0.5] ,'P[s|x]')
axis xy
title('Continuous Rate Tuning')
bz_piTickLabel('x')

NiceSave('TuningCurves',figfolder,'Modal_VMPoisson')


%% Compare parameters: dt

k=1;
pAS_0 = 0.7;
pAS_pi = 0.05;

dt = 0.05;
rAS = 10;
rGS = 1;


[pS_x_modal,pAS_x,pS_AS,pS_GS] = EncodingModal_VMPoisson(S,X,x0,k,pAS_0,pAS_pi,rAS,rGS,dt);

%[pS_x_modal,pAS_x,pS_AS,pS_GS] = EncodingContinuous_VMPoisson(S,X,x0,k,pAS_0,pAS_pi,rAS,rGS,dt);
R_0 = pAS_0.*rAS + (1-pAS_0).*rGS;
R_pi = pAS_pi.*rAS + (1-pAS_pi).*rGS;
[pS_x_cont,R_x] = EncodingContinuous_VMPoisson(S,X,x0,k,R_0,R_pi,dt);


%%

figure
subplot(2,2,1)
    %plot(x,pAS_x(1,:))
    plot(x,pAS_x(1,:))
    xlabel('x');ylabel('P[AS|x]') 
    ylim([0 1])
    xlim([-pi pi])
    bz_piTickLabel('x')
    yyaxis right
    plot (x,R_x(1,:),'--')
    ylim([rGS rAS])
    ylabel('R[x]') 
    
subplot(3,2,2)
    plot(s,pS_AS(:,1))
    hold on
    plot(s,pS_GS(:,1))
    xlabel(['# Spikes in dt=',num2str(dt),' bin (s)'])
    ylabel('P(s|state)')
    legend(['AS: ',num2str(rAS),'Hz'],['GS: ',num2str(rGS),'Hz'])

subplot(2,2,3)
imagesc(x,s,pS_x_modal)
xlabel('x');ylabel('s')
ColorbarWithAxis([0 0.5] ,'P[s|x]')
axis xy
title('Probabilistic Modal Tuning')
bz_piTickLabel('x')

subplot(2,2,4)
imagesc(x,s,pS_x_cont)
xlabel('x');ylabel('s')
ColorbarWithAxis([0 0.5] ,'P[s|x]')
axis xy
title('Continuous Rate Tuning')
bz_piTickLabel('x')

NiceSave('TuningCurves_Lowdt',figfolder,'Modal_VMPoisson')

%%
k=1;
pAS_0 = 0.8;
pAS_pi = 0.1;

dt = 0.5;
rAS = 5;
rGS = 0.1;


[pS_x_modal,pAS_x,pS_AS,pS_GS] = EncodingModal_VMPoisson(S,X,x0,k,pAS_0,pAS_pi,rAS,rGS,dt);

%[pS_x_modal,pAS_x,pS_AS,pS_GS] = EncodingContinuous_VMPoisson(S,X,x0,k,pAS_0,pAS_pi,rAS,rGS,dt);
R_0 = pAS_0.*rAS + (1-pAS_0).*rGS;
R_pi = pAS_pi.*rAS + (1-pAS_pi).*rGS;
[pS_x_cont,R_x] = EncodingContinuous_VMPoisson(S,X,x0,k,R_0,R_pi,dt);


%%

figure
subplot(2,2,1)
    %plot(x,pAS_x(1,:))
    plot(x,pAS_x(1,:))
    xlabel('x');ylabel('P[AS|x]') 
    ylim([0 1])
    xlim([-pi pi])
    bz_piTickLabel('x')
    yyaxis right
    plot (x,R_x(1,:),'--')
    ylim([rGS rAS])
    ylabel('R[x]') 
    
subplot(3,2,2)
    plot(s,pS_AS(:,1))
    hold on
    plot(s,pS_GS(:,1))
    xlabel(['# Spikes in dt=',num2str(dt),' bin (s)'])
    ylabel('P(s|state)')
    legend(['AS: ',num2str(rAS),'Hz'],['GS: ',num2str(rGS),'Hz'])

subplot(2,2,3)
imagesc(x,s,pS_x_modal)
xlabel('x');ylabel('s')
ColorbarWithAxis([0 0.5] ,'P[s|x]')
axis xy
title('Probabilistic Modal Tuning')
bz_piTickLabel('x')

subplot(2,2,4)
imagesc(x,s,pS_x_cont)
xlabel('x');ylabel('s')
ColorbarWithAxis([0 0.5] ,'P[s|x]')
axis xy
title('Continuous Rate Tuning')
bz_piTickLabel('x')

%NiceSave('TuningCurves_Lowdt',figfolder,'Modal_VMPoisson')