reporoot = '/Users/dlevenstein/Project Repos/NeuronalHeterogeneity/';
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/GroundStateAnalysis'];


[GroundStateAll,baseNames] = GetMatResults(figfolder,'GroundStateAnalysis','select',true);
GroundStateAll = bz_CollapseStruct(GroundStateAll);
thisregion = 'CA1';


%%
