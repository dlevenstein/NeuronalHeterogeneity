reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/'; %Laptop
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/HeadDirectionEncodingAnalysis'];


[HDALL,baseNames] = GetMatResults(figfolder,'HeadDirectionEncodingAnalysis');
HDALL = bz_CollapseStruct(HDALL);

%%
allcells_c = bz_CollapseStruct(HDALL.model_c,'match','justcat',true);
allcells_m = bz_CollapseStruct(HDALL.model_m,'match','justcat',true);
%ISIbyHD_alignGam = bz_CollapseStruct(HDALL.ISIbyHD_alignGam,3,'justcat',true);
% MutInfo = bz_CollapseStruct(HDALL.MutInfo,'match','justcat',true);
% cellISIStats = bz_CollapseStruct(HDALL.cellISIStats,3,'justcat',true);
%%
figure
subplot(3,4,1)
BoxAndScatterPlot(allcells_m.BIC-allcells_c.BIC)
hold on
box off
UnityLine
ylabel('Modal-Continuous Model BIC')
subplot(3,3,2)
plot(allcells_m.parms.pAS_0,allcells_m.parms.pAS_pi,'.')
hold on
UnityLine
xlabel('p[AS|x=x_0]');ylabel('p[AS|x=x_0+pi]')

subplot(3,3,3)
plot(allcells_m.parms.rAS,allcells_m.parms.pAS_0,'.')
hold on
%UnityLine
xlabel('R_A_S');ylabel('p[AS|x=x_0]')

subplot(3,3,4)
plot(allcells_m.parms.rGS,allcells_m.parms.pAS_pi,'.')
hold on
%UnityLine
box off
xlabel('R_G_S');ylabel('p[AS|x=x_0+pi]')

% subplot(3,3,5)
% plot(log10(allcells_m.parms.rGS),GammaFit.WAKEstate.sharedfit.GSlogrates(HDcells_GammaIDX),'.')
% LogScale('xy',10)
% box off
% xlabel('GS Rate: Encoding Model');ylabel('GS Rate: Gamma ISI Model')
NiceSave('EncodingModelStats',figfolder,[])
