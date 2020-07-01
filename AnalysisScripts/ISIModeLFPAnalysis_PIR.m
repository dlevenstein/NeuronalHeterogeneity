function [ModalLFPModulation] = ISIModeLFPAnalysis_PIR(basePath,figfolder)
% Date XX/XX/20XX
%
%Question: 
%
%Plots
%-
%-
%
%% Load Header
%Initiate Paths
reporoot = '/gpfs/data/buzsakilab/DL/NeuronalHeterogeneity/';
%reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/';
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/Cicero_09102014';
% %basePath = pwd;
% %basePath = fullfile(reporoot,'Datasets/onProbox/AG_HPC/Achilles_11012013');
%figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%Load Stuff
sessionInfo = bz_getSessionInfo(basePath);
spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
states = fieldnames(SleepState.ints);
states{4} = 'ALL';
SleepState.ints.ALL = [0 Inf];
statecolors = {[0 0 0],[0 0 1],[1 0 0],[0.6 0.6 0.6]};

try
celltypes = CellClass.celltypes;
catch
    celltypes = unique(CellClass.label);
end
cellcolor = {'k','r'};


%%

LFPMapFolder = [reporoot,'AnalysisScripts/AnalysisFigs/ISILFPMap'];

%Check for an LFP Map
try
    [ISILFPMap] = GetMatResults(LFPMapFolder,'ISILFPMap','baseNames',baseName);
catch
    error('No Channel selected')
end

region = 'pir';
%If One exists: bz_tagChannel
lfpchannel = ISILFPMap.MIMap.selectedchans.(region).channel;
usecells = ISILFPMap.MIMap.(ISILFPMap.MIMap.selectedchans.(region).regname).UIDs;
%If it doesnt, skip all the mode stuff and just run bz_ISILFPMap
%Question: what to do about pir/bla
%Note: need to only use UIDs from the right region

    %%
    %lfpchannel = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.THchanID; 
    %Load from the analysisresults. and tag channel.
    %Then run [ISILFPMap] =
    %bz_ISILFPMap(basePath,varargin) to save again and save with channel
    %marking? (also 
    %In the future, here, tagChannel for each recording!
%% Load the LFP
downsamplefactor = 1;
lfp = bz_GetLFP(lfpchannel,...
    'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);

%%
  
%% Calculate Residuals
dt = 0.01;
winsize = 10;
frange = [1 312];
nfreqs = 150;
[specslope] = bz_PowerSpectrumSlope(lfp,winsize,dt,'spectype','wavelet',...
    'nfreqs',nfreqs,'showfig',true,'showprogress',true,'frange',frange,...
    'saveMat',basePath,'saveName',['Chan',num2str(lfpchannel)],...
    'saveFolder','WavPSS','Redetect',true);
    

%% Zoom on Theta (load?)
    %Load Gamma Modes
    GammaFit = bz_LoadCellinfo(basePath,'GammaFit');



%%

pc = parcluster('local');

    % sto
% % store temporary files in the 'scratch' drive on the cluster, labeled by job ID
pc.JobStorageLocation = strcat(getenv('SCRATCH'), '/', getenv('SLURM_JOB_ID'));
% % enable MATLAB to utilize the multiple cores allocated in the job script
% % SLURM_NTASKS_PER_NODE is a variable set in the job script by the flag --tasks-per-node
% % we use SLURM_NTASKS_PER_NODE - 1, because one of these tasks is the original MATLAB script itself
parpool(pc, str2num(getenv('SLURM_NTASKS_PER_NODE'))-1);

%% Ex cell: all frequencies
specslope = rmfield(specslope,'specgram');
specslope = rmfield(specslope,'phase');

%%
%clear fPower
clear FConditionalISIDist

numAS = GammaFit.WAKEstate.parms.numAS;



%numCells = length(usecells);
numstates = 2;
MutInf = nan(length(specslope.freqs),spikes.numcells,numstates);
GSModulation = nan(length(specslope.freqs),spikes.numcells,numstates);
ASModulation = nan(length(specslope.freqs),spikes.numcells,numAS,numstates);
ASlogRates = nan(length(specslope.freqs),spikes.numcells,numAS,numstates);

fPower.timestamps = specslope.timestamps;
%fPower(length(specslope.freqs)+1).data = [];
%FConditionalISIDist = struct
for ff = 1:length(specslope.freqs)
    bz_Counter(ff,length(specslope.freqs),'Freq')
    fPower.data = specslope.resid(:,ff);

    parfor cc = 1:spikes.numcells
        %cc = excell;
        cellUID(cc) = spikes.UID(cc);
        
       % state = 'WAKEstate';
        for ss =1:2
        %Find the UID of the cell in the Gamma fit so match...
        %Put the gamma fit parms to conditional dist in as initial parms
        GFIDX = find(GammaFit.(states{ss}).cellstats.UID==cellUID(cc));
        cellGamma = GammaFit.(states{ss}).singlecell(GFIDX);
        if length(GFIDX)==1 && ismember(cellUID(cc),usecells)

        
            [FConditionalISIDist] = bz_ConditionalISI(spikes.times{cc},fPower,...
                'ints',SleepState.ints.(states{ss}),'GammaFitParms',cellGamma,...
                'showfig',false,'GammaFit',true,'minX',0);

            MutInf(ff,cc,ss) = FConditionalISIDist.MutInf;
            %GammaModes(ff,cc) = FConditionalISIDist.GammaModes;

            GSModulation(ff,cc,ss) = FConditionalISIDist.GammaModes.GS_R;
            ASModulation(ff,cc,:,ss) = FConditionalISIDist.GammaModes.AS_R;
            ASlogRates(ff,cc,:,ss) = FConditionalISIDist.GammaModes.ASlogrates;
            ASweight(1,cc,:,ss) = cellGamma.ASweights;
        end
        end
    end
end
%%
for ss = 1:2
    ModalLFPModulation.(states{ss}).MutInf = MutInf(:,:,ss);
    ModalLFPModulation.(states{ss}).GSModulation = GSModulation(:,:,ss);
    ModalLFPModulation.(states{ss}).ASModulation = ASModulation(:,:,:,ss);
    ModalLFPModulation.(states{ss}).ASlogRates = ASlogRates(:,:,:,ss);
    ModalLFPModulation.(states{ss}).ASweight = ASweight(:,:,:,ss);

end
ModalLFPModulation.UID = cellUID;
ModalLFPModulation.freq = specslope.freqs;
%save('FConditionalISIDist.mat','FConditionalISIDist')
%%
% for cc = 1:spikes.numcells
%     %AllFConditionalISIDist(cc) = bz_CollapseStruct(FConditionalISIDist(:,cc));
%     AllFConditionalISIModes(cc) = bz_CollapseStruct(GammaModes(:,cc),2);
%     AllFConditionalISIModes_rates(cc) = bz_CollapseStruct(GammaModes(:,cc),1);
%     AllFConditionalISIModes_ASweights(cc) = bz_CollapseStruct(GammaModes(:,cc),3);
% 
%     %AllFConditionalMINF(cc) = bz_CollapseStruct(AllFConditionalISIDist(cc).MutInf,2);
%     if length(AllFConditionalISIModes(cc).GS_R)~=150
%         AllFConditionalISIModes(cc).GS_R = [];
%         AllFConditionalISIDist(cc).MutInf = [];
%     end
% end
% 
% AllFConditionalISIModes_mean = bz_CollapseStruct(AllFConditionalISIModes,1,'mean');
% %%
% AllFConditionalMInf = mean(MutInf,2);

%%
%%

figure
for ss = 1:2
subplot(2,2,ss)
hold on
for cc = 1:spikes.numcells
for rr = 1:5
    modeoccupancy = ModalLFPModulation.(states{ss}).ASweight(1,cc,rr);
    if (modeoccupancy<0.02)
        continue
    end
    %showwhichmodes = (AllFConditionalISIModes(cc).AScorr_p(rr,:)<1); %& ...
        %showwhichmodes = (modeoccupancy>0.05);

        
    scatter(log10(ModalLFPModulation.freq),...
        ModalLFPModulation.(states{ss}).ASlogRates(:,cc,rr),...
        modeoccupancy*30,ModalLFPModulation.(states{ss}).ASModulation(:,cc,rr),'.')
end
end
UnityLine
axis tight
LogScale('xy',10)
ColorbarWithAxis([-0.0 0.15],'Mode-Power Correlation','inclusive',{'<','>'})
crameri('vik','pivot',0)
title(states{ss})
xlabel('LFP frequency (Hz)');ylabel('ISI Mode Rate (Hz)')
end

subplot(4,2,7)
hold on
for ss = 1:2
plot(log10(ModalLFPModulation.freq),nanmean(ModalLFPModulation.(states{ss}).MutInf,2),'color',statecolors{ss},'linewidth',1)
end
LogScale('x',10)
xlabel('f (Hz)')
ylabel('MI[ISI;Power]')
colorbar
axis tight


subplot(4,2,8)
hold on
for ss = 1:2
plot(log10(ModalLFPModulation.freq),-nanmean(ModalLFPModulation.(states{ss}).GSModulation,2),'color',statecolors{ss},'linewidth',2)
end
plot(xlim(gca),[0 0],'k--')
LogScale('x',10)
xlabel('f (Hz)')
ylabel('AR-Power Corr')
colorbar
axis tight


NiceSave('ModeLFPFreqMod',figfolder,baseName)

%%
% figure
% subplot(2,2,1)
% plot(log10(specslope.freqs),-AllFConditionalISIModes.GS_R,'k')
% hold on
% %plot(log10(specslope.freqs(AllFConditionalISIModes.GScorr_p<0.05)),AllFConditionalISIModes.GS_R(AllFConditionalISIModes.GScorr_p<0.05),'o')
% hold on
% plot(xlim(gca),[0 0],'k--')
% LogScale('x',10)
% ylabel('AR-Power Corr')
% xlabel('f (Hz)')
% box off 
% 
% subplot(2,2,2)
% plot(log10(specslope.freqs),AllFConditionalISIDist.MutInf,'k')
% hold on
% LogScale('x',10)
% ylabel('MutInf')
% xlabel('f (Hz)')
% box off
% 
% 
% NiceSave(['PowerModUID',num2str(excellUID)],figfolder,baseName)





    

end
