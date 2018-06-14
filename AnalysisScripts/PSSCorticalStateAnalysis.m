function [ output_args ] = PSSCorticalStateAnalysis( basePath,figfolder )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% DEV
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
basePath = '/mnt/proraidDL/Database/BWCRCNS/JennBuzsaki22/20140526_277um';
figfolder = '/Users/dlevenstein/Dropbox/Research/Current Projects/FRHetAndDynamics/AnalysisScripts/AnalysisFigs';
figfolder = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/AnalysisScripts/AnalysisFigs';
%%
baseName = bz_BasenameFromBasepath(basePath);


spikes = bz_GetSpikes('basePath',basePath);
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
lfp = bz_GetLFP(SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID,...
    'basepath',basePath);
%%
downsamplefactor = 5;
lfp_down = bz_DownsampleLFP(lfp,downsamplefactor);
%%
dt = 0.5;

numwins = 10;
winsizes = logspace(-0.301,1.5,numwins);
for ww = numwins:-1:1
    ww
    [specslope_temp] = bz_PowerSpectrumSlope(lfp_down,winsizes(ww),dt,'showfig',false);
    if ww == numwins
        specslope=specslope_temp;
        specslope.winsizes = winsizes;
    else
        specslope.data(:,ww)= interp1(specslope_temp.timestamps,specslope_temp.data,specslope.timestamps,'nearest');
        specslope.rsq(:,ww) = interp1(specslope_temp.timestamps,specslope_temp.rsq,specslope.timestamps,'nearest');
        specslope.intercept(:,ww) = interp1(specslope_temp.timestamps,specslope_temp.intercept,specslope.timestamps,'nearest');
    end
end

%%
figure
hist(specslope.rsq)
