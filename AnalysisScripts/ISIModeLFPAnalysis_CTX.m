function [PopConditional,AllFConditionalISIDist] = ISIModeLFPAnalysis_CTX(basePath,figfolder)
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
    lfpchannel = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID; 

%% Load the LFP
downsamplefactor = 2;
lfp = bz_GetLFP(lfpchannel,...
    'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);

%%

% ss = 1;
% state = states{ss};
% %Take only subset of time (random intervals) so wavelets doesn't break
% %computer (total 625s)
% usetime = 7200;%2500
% winsize = 25;
% if sum(diff(SleepState.ints.(state),1,2))>usetime
%     nwin = round(usetime./winsize);
%     %winsize = 30; %s
%     try
%         windows = bz_RandomWindowInIntervals( SleepState.ints.(state),winsize,nwin );
%     catch
%         windows = SleepState.ints.(state);
%     end
% else
%     windows = SleepState.ints.(state);
% end
        
%% Calculate Residuals
dt = 0.01;
winsize = 10;
frange = [1 312];
nfreqs = 150;
[specslope] = bz_PowerSpectrumSlope(lfp,winsize,dt,'spectype','wavelet',...
    'nfreqs',nfreqs,'showfig',true,'showprogress',true,'frange',frange,...
    'saveMat',basePath,'saveName',['Chan',num2str(lfpchannel)],...
    'saveFolder','WavPSS');
    
clear lfp




%% %Load Gamma Modes
    GammaFit = bz_LoadCellinfo(basePath,'GammaFit');
    %Normalize the brain state metrics
%     ThetaPower.timestamps = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.t_clus;
%     ThetaPower.data = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.thratio;


%% Ex cell: all frequencies
specslope = rmfield(specslope,'specgram');
specslope = rmfield(specslope,'phase');

%%
%clear fPower
clear FConditionalISIDist
fPower.timestamps = specslope.timestamps;
%fPower(length(specslope.freqs)+1).data = [];
for ff = 1:length(specslope.freqs)
    bz_Counter(ff,length(specslope.freqs),'Freq')
    fPower.data = specslope.resid(:,ff);

    parfor cc = 1:spikes.numcells
        %cc = excell;
    excellUID = spikes.UID(cc);
    %Find the UID of the cell in the Gamma fit so match...
    %Put the gamma fit parms to conditional dist in as initial parms
    GFIDX = find(GammaFit.WAKEstate.cellstats.UID==excellUID);
    cellGamma = GammaFit.WAKEstate.singlecell(GFIDX);
    if ~isempty(cellGamma)
        [FConditionalISIDist(ff,cc)] = ConditionalISI(spikes.times{cc},fPower,...
            'ints',SleepState.ints.WAKEstate,'GammaFitParms',cellGamma,...
            'showfig',false,'GammaFit',true);
    end
   % catch
    %    continue
   % end
    end
%catch
%    continue
%end
end
%%
%save('FConditionalISIDist.mat','FConditionalISIDist')
%%
for cc = 1:spikes.numcells
AllFConditionalISIDist(cc) = bz_CollapseStruct(FConditionalISIDist(:,cc));
AllFConditionalISIModes(cc) = bz_CollapseStruct(AllFConditionalISIDist(cc).GammaModes,2);
AllFConditionalISIModes_rates(cc) = bz_CollapseStruct(AllFConditionalISIDist(cc).GammaModes,1);
AllFConditionalISIModes_ASweights(cc) = bz_CollapseStruct(AllFConditionalISIDist(cc).GammaModes,3);

%AllFConditionalMINF(cc) = bz_CollapseStruct(AllFConditionalISIDist(cc).MutInf,2);
    if length(AllFConditionalISIModes(cc).GS_R)~=150
        AllFConditionalISIModes(cc).GS_R = [];
        AllFConditionalISIDist(cc).MutInf = [];
    end
end

AllFConditionalISIModes_mean = bz_CollapseStruct(AllFConditionalISIModes,1,'mean');
%%
AllFConditionalMInf = mean(cat(1,AllFConditionalISIDist(:).MutInf));

%%
%AllFConditionalISIDist = bz_CollapseStruct(AllFConditionalISIDist,2);
%%
figure
subplot(2,2,1)
hold on
for cc = 1:spikes.numcells
for rr = 1:5
    modeoccupancy = squeeze(mean(AllFConditionalISIModes_ASweights(cc).ASweights(:,rr,:),1));
    if mean(modeoccupancy<0.02)
        continue
    end
    showwhichmodes = (AllFConditionalISIModes(cc).AScorr_p(rr,:)<1); %& ...
        %showwhichmodes = (modeoccupancy>0.05);

        
    scatter(log10(specslope.freqs(showwhichmodes)),...
        AllFConditionalISIModes_rates(cc).ASlogrates(showwhichmodes,rr),...
        modeoccupancy(showwhichmodes)*30,AllFConditionalISIModes(cc).AS_R(rr,showwhichmodes),'.')
end
end
UnityLine
axis tight
LogScale('xy',10)
ColorbarWithAxis([-0.0 0.15],'Mode-Power Correlation','inclusive',{'<','>'})
crameri('vik','pivot',0)
xlabel('LFP frequency (Hz)');ylabel('ISI Mode Rate (Hz)')

subplot(4,2,7)
hold on
for ss = 1
plot(log10(specslope.freqs),AllFConditionalMInf,'color',statecolors{ss},'linewidth',1)
end
LogScale('x',10)
xlabel('f (Hz)')
ylabel('MI[ISI;Power]')
colorbar
axis tight


subplot(4,2,8)
hold on

plot(log10(specslope.freqs),-AllFConditionalISIModes_mean.GS_R,'color',statecolors{ss},'linewidth',2)
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




%% Theta ISI modulation - all cells

% parfor cc = 1:spikes.numcells
%     cc
% excellUID = spikes.UID(cc);
% %Find the UID of the cell in the Gamma fit so match...
% %Put the gamma fit parms to conditional dist in as initial parms
% GFIDX = find(GammaFit.WAKEstate.cellstats.UID==excellUID);
% if isempty(GFIDX)
%     continue
% end
% cellGamma = GammaFit.WAKEstate.singlecell(GFIDX);
% try
% [ThetaConditionalISI(cc)] = ConditionalISI(spikes.times{cc},ThetaPower,...
%     'ints',SleepState.ints.WAKEstate,'GammaFitParms',cellGamma,...
%     'showfig',false,'GammaFit',true);
% catch
%     continue
% end
% end
% 
% %%
% AllThetaISIModes = bz_CollapseStruct(ThetaConditionalISI);
% AllThetaISIModes = bz_CollapseStruct(AllThetaISIModes.GammaModes,2);
% %% Figure: Theta ISI modulation all cells
% figure
% subplot(2,2,1)
%     hist(AllThetaISIModes.GSCorr)
% subplot(2,2,2)
%     plot(AllThetaISIModes.ASlogrates(AllThetaISIModes.AScorr_p<0.05),...
%         AllThetaISIModes.AS_R(AllThetaISIModes.AScorr_p<0.05),'k.')
%     hold on
%     plot(AllThetaISIModes.ASlogrates(AllThetaISIModes.AScorr_p>=0.05),...
%         AllThetaISIModes.AS_R(AllThetaISIModes.AScorr_p>=0.05),'.','color',[0.5 0.5 0.5])
% 
%     
%     plot(AllThetaISIModes.GSlogrates(AllThetaISIModes.GScorr_p<0.05),...
%         AllThetaISIModes.GS_R(AllThetaISIModes.GScorr_p<0.05),'.')
%     axis tight
%     box off
%         plot(xlim(gca),[0 0],'k--')
%         LogScale('x',10)
%         xlabel('Mode Rate (Hz)')
%         ylabel('Weight-Power Corr')
%     
%     
% subplot(2,2,3)
%     plot(AllThetaISIModes.GS_R,...
%         log10(AllThetaISIModes.GScorr_p),'o')
%     
%     
%  NiceSave('ThetaMod',figfolder,baseName)
% 
%     

end
