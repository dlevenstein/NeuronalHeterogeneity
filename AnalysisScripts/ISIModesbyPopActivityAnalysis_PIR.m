function [PopCorr,MUAConditionalISIDist,MUAConditionalISIDist_gamma,ISIdists] = ...
    ISIModesbyPopActivityAnalysis_PIR(basePath,figfolder)

%% DEV
%reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/';
reporoot = '/gpfs/data/buzsakilab/DL/NeuronalHeterogeneity/';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/mnt/proraidDL/Database/BWCRCNS/JennBuzsaki22/20140526_277um';
%basePath = '/mnt/proraidDL/Database/AGData/Cicero/Cicero_09012014';
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/Cicero_09102014';
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/Mouse24-131213';
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/Rat08-20130713';
%basePath = pwd;
%figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/ISIModesbyPopActivityAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
CellClass = bz_LoadCellinfo(basePath,'CellClass');
%OccupancyStats = bz_LoadCellinfo(basePath,'OccupancyStats');
SleepState = bz_LoadStates(basePath,'SleepState');
SleepState.ints.ALL = [0 Inf];
% lfp = bz_GetLFP(SleepState.detectorinfo.detectionparms.SleepScoreMetrics.SWchanID,...
%     'basepath',basePath);
%% Getting the right cells hack

LFPMapFolder = [reporoot,'AnalysisScripts/AnalysisFigs/ISILFPMap'];
[ISILFPMap] = GetMatResults(LFPMapFolder,'ISILFPMap','baseNames',baseName);
region = 'pir';
lfpchannel = ISILFPMap.MIMap.selectedchans.(region).channel;
usecells = ISILFPMap.MIMap.(ISILFPMap.MIMap.selectedchans.(region).regname).UIDs;
spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true,'UID',usecells);


CellClass.keep = ismember(CellClass.UID,usecells);
CellClass.pE = CellClass.pE(CellClass.keep);
CellClass.pI = CellClass.pI(CellClass.keep);
CellClass.UID = CellClass.UID(CellClass.keep);
CellClass.label = CellClass.label(CellClass.keep);

%%
cellinfo.CellClass = CellClass;
%cellinfo.ISIStats = ISIStats.summstats;
%cellinfo.OccupancyStats = OccupancyStats;
%% Load The Gamma stuff
GammaFit = bz_LoadCellinfo(basePath,'GammaFit');


%% Cell types and states
% try
% celltypes = CellClass.celltypes;
% catch
%     celltypes = unique(CellClass.label);
% end
celltypes = {'pE','pI'};
cellcolor = {'k','r'};
statenames = fieldnames(SleepState.ints);

%% Calculate spike count matrix
binsize = 0.06; %s
dt = 0.01;
spikemat = bz_SpktToSpkmat(spikes,'binsize',binsize,'dt',dt,'bintype','gaussian','units','rate');
spikemat.isspike = spikemat.data>0.5;
%spikemat.isspike(spikemat.isspike>1) = 1;
for ss = 1:3
    spikemat.instate.(statenames{ss}) = InIntervals(spikemat.timestamps,SleepState.ints.(statenames{ss}));
end
%% For each cell, calculate E and I pop rates of all OTHER cells
ncellthresh = 5;

for tt = 1:length(celltypes)
    Ncells.(celltypes{tt}) = sum(CellClass.(celltypes{tt}));
    spikemat.totpoprate.(celltypes{tt}) = sum(spikemat.data(:,CellClass.(celltypes{tt})),2);
    spikemat.poprate.(celltypes{tt}) = spikemat.totpoprate.(celltypes{tt})./Ncells.(celltypes{tt});
    spikemat.totpopsynch.(celltypes{tt}) = sum(spikemat.isspike(:,CellClass.(celltypes{tt})),2);
    spikemat.popsynch.(celltypes{tt}) = spikemat.totpopsynch.(celltypes{tt})./Ncells.(celltypes{tt});
    
%     if Ncells.(celltypes{tt})==0
%         celltypes(tt) = [];
%     end
end
Ncells.ALL = Ncells.pE + Ncells.pI;
spikemat.totpoprate.ALL = sum(spikemat.data(:,(CellClass.pI|CellClass.pE)),2);
spikemat.poprate.ALL = spikemat.totpoprate.ALL./Ncells.ALL;

for cc = 1:spikes.numcells %weird roundabout way to calculate is much faster
    bz_Counter(cc,spikes.numcells,'Cell');
    thiscell = false(size(CellClass.pE));
    thiscell(cc) = true;
    spikemat.cellrate{cc} = spikemat.data(:,cc);
    spikemat.cellspike{cc} = spikemat.isspike(:,cc); %New Way Below
%     spikemat.cellspike{cc} = spikemat.cellrate{cc};
%     spikemat.cellspike{cc}(spikemat.cellspike{cc}>1) = 1;
%     spikemat.cellspike{cc}(spikemat.cellspike{cc}<0) = 0;
    for tt = 1:length(celltypes)
        PopCorr.cellcount.(celltypes{tt})(cc) = sum(CellClass.(celltypes{tt}) & ~thiscell);
        if CellClass.(celltypes{tt})(cc) %if it's in theclass, subtract off the current cell
            spikemat.bycellpoprate.(celltypes{tt}){cc} = ...
                (spikemat.totpoprate.(celltypes{tt})-spikemat.cellrate{cc})./...
               PopCorr.cellcount.(celltypes{tt})(cc);
           
            spikemat.bycellpopsynch.(celltypes{tt}){cc} = ...
                (spikemat.totpopsynch.(celltypes{tt})-spikemat.cellspike{cc})./...
               PopCorr.cellcount.(celltypes{tt})(cc);
           %Here: add tiny random amount... for percentile...
           spikemat.bycellpopsynch.(celltypes{tt}){cc} = ...
               spikemat.bycellpopsynch.(celltypes{tt}){cc}+0.01*rand(size(spikemat.bycellpopsynch.(celltypes{tt}){cc}));
        else
            spikemat.bycellpoprate.(celltypes{tt}){cc} = ...
                spikemat.totpoprate.(celltypes{tt})./PopCorr.cellcount.(celltypes{tt})(cc);
            
            spikemat.bycellpopsynch.(celltypes{tt}){cc} = ...
                spikemat.totpopsynch.(celltypes{tt})./PopCorr.cellcount.(celltypes{tt})(cc);
            %Here: add tiny random amount... for percentile...
            spikemat.bycellpopsynch.(celltypes{tt}){cc}+0.01*rand(size(spikemat.bycellpopsynch.(celltypes{tt}){cc}));
        end
    end
    
    if CellClass.pI(cc)||CellClass.pE(cc)
        spikemat.bycellpoprate.ALL{cc} = (spikemat.totpoprate.ALL-spikemat.cellrate{cc})./...
            (Ncells.ALL-1);
    else
        spikemat.bycellpoprate.ALL{cc} = spikemat.totpoprate.ALL./Ncells.ALL;
    end
    
    

end

%%
% figure
% plot(spikemat.poprate.(celltypes{1}),spikemat.popsynch.(celltypes{1}),'.')
%% Correlation and Gamma parms

for cc = 1:spikes.numcells
    bz_Counter(cc,spikes.numcells,'Cell');
    for ss = 1:3
        %Calculate Correlation
        for tt = 1:length(celltypes)
            PopCorr.(statenames{ss}).rate.(celltypes{tt})(cc) = ...
                corr(spikemat.cellrate{cc}(spikemat.instate.(statenames{ss})),...
                spikemat.bycellpoprate.(celltypes{tt}){cc}(spikemat.instate.(statenames{ss})),...
                'type','spearman');
            PopCorr.(statenames{ss}).synch.(celltypes{tt})(cc) = ...
                corr(spikemat.cellspike{cc}(spikemat.instate.(statenames{ss})),...
                spikemat.bycellpopsynch.(celltypes{tt}){cc}(spikemat.instate.(statenames{ss})),...
                'type','spearman');
        end
    
        %Get the mean rate
        PopCorr.(statenames{ss}).meanRate(cc) = ISIStats.summstats.(statenames{ss}).meanrate(cc);
        %Get the GS rate
        cellUID(cc) = spikes.UID(cc);
        GFIDX = find(GammaFit.(statenames{ss}).cellstats.UID==cellUID(cc));
        if isempty(GFIDX)
            PopCorr.(statenames{ss}).GSrate(cc) = nan;
            PopCorr.(statenames{ss}).GSCV(cc) = nan;
            PopCorr.(statenames{ss}).GSweight(cc) = nan;
            continue
        end
        cellGamma = GammaFit.(statenames{ss}).singlecell(GFIDX);
        PopCorr.(statenames{ss}).GSrate(cc) = cellGamma.GSlogrates;
        PopCorr.(statenames{ss}).GSCV(cc) = cellGamma.GSCVs;
        PopCorr.(statenames{ss}).GSweight(cc) = cellGamma.GSweights;
    end
end
PopCorr.CellClass = CellClass;

%%
synchrate = {'rate','synch'};
%%
figure
for tt = 1:length(celltypes)
    for ss = 1:3
        for sr = 1:2
        subplot(4,3,(tt-1)*3+(sr-1).*6+ss)
        plot(PopCorr.(statenames{ss}).GSrate,PopCorr.(statenames{ss}).(synchrate{sr}).(celltypes{tt}),'k.')
        hold on
        plot(log10(PopCorr.(statenames{ss}).meanRate(CellClass.pI)),...
            PopCorr.(statenames{ss}).(synchrate{sr}).(celltypes{tt})(CellClass.pI),'r.')
        axis tight
        box off
        plot(xlim(gca),[0 0],'k--')
        xlabel('GS/Mean Rate (Hz)');
        LogScale('x',10,'exp',true,'nohalf',true)
        if ss ==1
            ylabel([(celltypes{tt}),' ',(synchrate{sr}),' Corr'])
        end
        if tt == 1
            title(statenames{ss})
        end
        end
%         subplot(4,3,(tt-1)*3+ss+6)
%         plot(PopCorr.(statenames{ss}).GSCV,PopCorr.(statenames{ss}).rate.(celltypes{tt}),'.')
%         hold on
%         axis tight
%         box off
%         plot(xlim(gca),[0 0],'k--')
%         xlabel('GS CV');ylabel([(celltypes{tt}),' Corr'])
%         %LogScale('x',10,'exp',false,'nohalf',true)
        
    end 
    
end

NiceSave('MUACorrandGSRate',figfolder,baseName)


%%
% for ss = 1:3
%     PopCorr.(statenames{ss}).GSmod = nan(size(MutInf.WAKEstate.PSS));
%     PopCorr.(statenames{ss}).GSmod_p = nan(size(MutInf.WAKEstate.PSS));
% end
%% Conditional ISI distrobution 
synchbins.pE = 10;
synchbins.pI = 10;
synchbounds.pE = [-1.5 -0.5];
synchbounds.pI = [-1 0];
%Should implement number cells threshold...
clear MUAConditionalISIDist_all
for cc = 1:spikes.numcells
    bz_Counter(cc,spikes.numcells,'Cell');    
    for tt = 1:length(celltypes)
%         if Ncells.(celltypes{tt})==0
%             continue
%         end

        %First synchrony
        MUA.timestamps = spikemat.timestamps;
        MUA.data = (spikemat.bycellpopsynch.(celltypes{tt}){cc});
        for ss = 1:3
            [MUAConditionalISIDist_all.(statenames{ss}).synch.(celltypes{tt})(cc)] = ...
                bz_ConditionalISI(spikes.times{cc},MUA,...
                'ints',SleepState.ints.(statenames{ss}),...
                'showfig',false,'GammaFit',false,'minX',50,'numISIbins',100);
                %'normtype','none','Xwin',synchbounds.(celltypes{tt}),'numXbins',synchbins.(celltypes{tt}) );
        end
        
        %Then Rate
        MUA.data = (spikemat.bycellpoprate.(celltypes{tt}){cc});
        for ss = 1:3
            [MUAConditionalISIDist_all.(statenames{ss}).rate.(celltypes{tt})(cc)] = ...
                bz_ConditionalISI(spikes.times{cc},MUA,...
                'ints',SleepState.ints.(statenames{ss}),...
                'showfig',false,'GammaFit',false,'minX',20,'numISIbins',100);
        
        
            %Gamma fit - this is expensive and inefficient to do it twice..
            cellUID(cc) = spikes.UID(cc);
            GFIDX = find(GammaFit.(statenames{ss}).cellstats.UID==cellUID(cc));
            if isempty(GFIDX)
                PopCorr.(statenames{ss}).(celltypes{tt}).GSmod(cc) = nan;
                PopCorr.(statenames{ss}).(celltypes{tt}).GSmod_p(cc) = nan;
                continue
            end

            cellGamma = GammaFit.(statenames{ss}).singlecell(GFIDX);
            try
            [MUAConditionalISIDist_gamma.(statenames{ss}).(celltypes{tt})(cc)] = ...
                bz_ConditionalISI(spikes.times(cc),MUA,...
                'ints',SleepState.ints.(statenames{ss}),...
                'GammaFitParms',cellGamma,'GammaFit',true,...
                'showfig',false,'numISIbins',100);
            catch

                continue
            end
            
            PopCorr.(statenames{ss}).(celltypes{tt}).GSmod(cc) = MUAConditionalISIDist_gamma.(statenames{ss}).(celltypes{tt})(cc).GammaModes.GS_R;
            PopCorr.(statenames{ss}).(celltypes{tt}).GSmod_p(cc) = MUAConditionalISIDist_gamma.(statenames{ss}).(celltypes{tt})(cc).GammaModes.GScorr_p;
        end
    end
end

%%
for ss = 1:3
    for tt = 1:2
    MUAConditionalISIDist_gamma.modes.(statenames{ss}).(celltypes{tt}) = bz_CollapseStruct([MUAConditionalISIDist_gamma.(statenames{ss}).(celltypes{tt}).GammaModes],3,'justcat');
    MUAConditionalISIDist_gamma.dist.(statenames{ss}).(celltypes{tt}) = bz_CollapseStruct([MUAConditionalISIDist_gamma.(statenames{ss}).(celltypes{tt}).Dist],3,'justcat');
    end
end

%%
figure
for tt = 1:2
for ss = 1:3
subplot(3,2,tt+2.*(ss-1))
hold on
keepcells = MUAConditionalISIDist_gamma.modes.(statenames{ss}).(celltypes{tt}).GScorr_p<=0.05;
%keepcells = true(size(keepcells))
for aa = 1:5
    %keepmodes = keepcells&mean(PSSConditionalGamma.modes.(states{ss}).ASweights(:,aa,:),1)>0.02;
    keepmodes = (MUAConditionalISIDist_gamma.modes.(statenames{ss}).(celltypes{tt}).AScorr_p(:,aa,:))<=0.05;
    %keepmodes = true(size(keepmodes))
scatter(-MUAConditionalISIDist_gamma.modes.(statenames{ss}).(celltypes{tt}).ASlogrates(1,aa,keepmodes),...
    log10(MUAConditionalISIDist_gamma.modes.(statenames{ss}).(celltypes{tt}).ASCVs(1,aa,keepmodes)),...
    60*mean(MUAConditionalISIDist_gamma.modes.(statenames{ss}).(celltypes{tt}).ASweights(:,aa,keepmodes),1)+eps,...
    squeeze(MUAConditionalISIDist_gamma.modes.(statenames{ss}).(celltypes{tt}).AS_R(1,aa,keepmodes)),'filled')
end
scatter(-MUAConditionalISIDist_gamma.modes.(statenames{ss}).(celltypes{tt}).GSlogrates(1,1,keepcells),...
    log10(MUAConditionalISIDist_gamma.modes.(statenames{ss}).(celltypes{tt}).GSCVs(1,1,keepcells)),...
    10,...
    squeeze(MUAConditionalISIDist_gamma.modes.(statenames{ss}).(celltypes{tt}).GS_R(1,1,keepcells)))
colorbar
axis tight
caxis([-0.15 0.15])
crameri('vik','pivot',0)
xlabel('Mean');ylabel('CV')
if ss == 1
    title(celltypes{tt})
end
end
end
NiceSave('ASModPopRate',figfolder,baseName)


%% Hi/Low PSS ints
RateThresh = [0.8 0.2];
clear ISIdists
for cc = 1:spikes.numcells
    bz_Counter(cc,spikes.numcells,'Cell');  
    for ss = 1:3
        for tt = 1:2
        instate = spikemat.instate.(statenames{ss});
        percentilenorm = NormToInt((spikemat.bycellpoprate.(celltypes{tt}){cc}(instate)),'percentile');

        IDX.timestamps = spikemat.timestamps(instate);
        IDX.states = zeros(size(IDX.timestamps));
        IDX.states(percentilenorm<=RateThresh(1)) = 1;
        IDX.states(percentilenorm>=RateThresh(2)) = 2;
        IDX.statenames = {'LowPopRate','HighPopRate'};
        INT = bz_IDXtoINT(IDX);

        spikestemp.times = spikes.times(cc);
        spikestemp.UID = spikes.UID(cc);
        tempstruct = bz_ISIStats(spikestemp,'ints',INT,'showfig',false);
        tempstruct.ISIhist.UID = spikes.UID(cc);
        

         tempstruct.ISIhist.summstats = ...
             tempstruct.summstats;
         %Threshold number of spikes for calculating return ap
    %Later - put this as temp and save the things you want (i.e. no
    %allspikes...)
        ISIdists.(statenames{ss}).(celltypes{tt})(cc) = tempstruct.ISIhist;
        end
    end
end
%%
for ss = 1:3
    for tt = 1:length(celltypes)
        ISIdists_temp.(statenames{ss}).(celltypes{tt}) = ...
            bz_CollapseStruct(ISIdists.(statenames{ss}).(celltypes{tt}),3,'justcat',true);
    for sr = 1:2
        MUAConditionalISIDist.(statenames{ss}).(synchrate{sr}).(celltypes{tt}) = ...
            bz_CollapseStruct(MUAConditionalISIDist_all.(statenames{ss}).(synchrate{sr}).(celltypes{tt}),3,'justcat',true);

        for tt2 = 1:length(celltypes)
            if sum(CellClass.(celltypes{tt2}))==0
                continue
            end
            MeanCondISI.(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}) = ...
                bz_CollapseStruct(MUAConditionalISIDist_all.(statenames{ss}).(synchrate{sr}).(celltypes{tt})(CellClass.(celltypes{tt2})),...
                3,'mean',true);
            
            numspksthresh = 300;
            numspks = squeeze(ISIdists_temp.(statenames{ss}).(celltypes{tt}).summstats.LowPopRatestate.numspikes)';
            MeanISIdists.(statenames{ss}).(celltypes{tt}).(celltypes{tt2}) = ...
                bz_CollapseStruct(ISIdists.(statenames{ss}).(celltypes{tt})(...
                CellClass.(celltypes{tt2})&numspks>numspksthresh),3,'mean',true);
        end
    end

    end
end
ISIdists = ISIdists_temp;
%%
lowhi = {'HighPopRatestate','LowPopRatestate'};

%%
for tt2 = 1:length(celltypes)
figure
for ss = 1:3
    for tt = 1:length(celltypes)
        
            for ll = 1:2
                
         subplot(4,6,ss+(tt-1)*3)
         hold on
         plot(MeanISIdists.(statenames{ss}).(celltypes{tt}).(celltypes{tt2}).logbins,...
             MeanISIdists.(statenames{ss}).(celltypes{tt}).(celltypes{tt2}).(lowhi{ll}).log)
         axis tight
         LogScale('x',10,'exp',true) 
                
         subplot(4,6,ss+(tt-1)*3+(ll-1)*6+6)
         imagesc(MeanISIdists.(statenames{ss}).(celltypes{tt}).(celltypes{tt2}).logbins,...
             MeanISIdists.(statenames{ss}).(celltypes{tt}).(celltypes{tt2}).logbins,...
             MeanISIdists.(statenames{ss}).(celltypes{tt}).(celltypes{tt2}).(lowhi{ll}).return)
         axis xy
         LogScale('xy',10,'exp',true)
            end
        
    end
end
NiceSave(['HiLoReturnMaps_',(celltypes{tt2})],figfolder,baseName)
end
%%
% excell = 13;
% figure
% subplot(2,2,1)
% imagesc(MUAConditionalISIDist.(statenames{2}).(synchrate{2}).(celltypes{1}).Dist.pYX(:,:,excell)')
% subplot(2,2,3)
% plot(MUAConditionalISIDist.(statenames{2}).(synchrate{2}).(celltypes{1}).Dist.Xocc(:,:,excell))
%% 
for sr = 1:2 
figure
for ss = 1:3
for tt = 1:length(celltypes) %Pop
    for tt2 = 1:length(celltypes) %Ref

        try
    subplot(4,3,(tt2-1)*6+(tt-1)*3+ss)
        imagesc(MeanCondISI.(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}).Dist.Xbins,...
            MeanCondISI.(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}).Dist.Ybins,...
            MeanCondISI.(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}).Dist.pYX')
        hold on
        plot(MeanCondISI.(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}).Dist.Xbins,...
            -log10(MeanCondISI.(statenames{ss}).(synchrate{sr}).(celltypes{tt}).(celltypes{tt2}).Dist.SpikeRate),...
            'r','LineWidth',2)
           
        xlabel([(celltypes{tt}),' ',(synchrate{sr})]);ylabel([(celltypes{tt2}),' ISI (s)'])
        LogScale('y',10,'exp',true,'nohalf',true)
        bz_AddRightRateAxis
        if tt ==1 & tt2 == 1
            title(statenames{ss})
        end
        catch
        end
    end 
end
end


NiceSave(['ISIbyMUA_',(synchrate{sr})],figfolder,baseName)

end
end