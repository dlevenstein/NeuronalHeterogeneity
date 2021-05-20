function [ ] = AnalysisXXXXXXXX(basePath,figfolder)
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
%reporoot = '/home/dlevenstein/ProjectRepos/NeuronalHeterogeneity/';
reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
basePath = '/Users/dl2820/Dropbox/Research/Datasets/Cicero_09102014';
%basePath = pwd;
%basePath = fullfile(reporoot,'Datasets/onProbox/AG_HPC/Cicero_09102014');
figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);

%Load Stuff
%sessionInfo = bz_getSessionInfo(basePath);
%spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);
%CellClass = bz_LoadCellinfo(basePath,'CellClass');
SleepState = bz_LoadStates(basePath,'SleepState');
ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
states = fieldnames(SleepState.ints);
%states{4} = 'ALL';
%SleepState.ints.ALL = [0 Inf];
%statecolors = {'k','b','r',[0.6 0.6 0.6]};



%%
load([basePath,'/GammaProcessed1/hmm_out.mat'])
GammaFit = bz_LoadCellinfo(basePath,'GammaFit');
%% Prepare the mode structures
ModeHMM.WAKEstate = WAKEall;
ModeHMM.NREMstate = NREMall;

numModes=7;
numcells = length(WAKEall);
spkthresh = 50;
MeanReturn.logbins = linspace(-3,2,50);
%get next ISI (nan for last one in the state)
%Cat all the cells
for ss = 1:2
for cc = 1:numcells

ModeHMM.(states{ss})(cc).next_isi = cellfun(@(prevISI) [prevISI(2:end) nan],ModeHMM.(states{ss})(cc).state_isi,'UniformOutput',false);
ModeHMM.(states{ss})(cc).next_state = cellfun(@(prevState) [prevState(2:end) nan],ModeHMM.(states{ss})(cc).decoded_mode,'UniformOutput',false);
ModeHMM.(states{ss})(cc).next_pMode = cellfun(@(prevProb) [prevProb(2:end,:);nan(1,numModes)],ModeHMM.(states{ss})(cc).prob_mode,'UniformOutput',false);


ModeHMM.(states{ss})(cc).prev_isi = cellfun(@(prevISI) prevISI(1:end),ModeHMM.(states{ss})(cc).state_isi,'UniformOutput',false);
ModeHMM.(states{ss})(cc).prev_state = cellfun(@(prevState) prevState(1:end),ModeHMM.(states{ss})(cc).decoded_mode,'UniformOutput',false);
ModeHMM.(states{ss})(cc).prev_pMode = cellfun(@(prevProb) prevProb(1:end,:),ModeHMM.(states{ss})(cc).prob_mode,'UniformOutput',false);


ModeHMM.(states{ss})(cc).prev_isi = cat(2,ModeHMM.(states{ss})(cc).prev_isi{:});
ModeHMM.(states{ss})(cc).next_isi = cat(2,ModeHMM.(states{ss})(cc).next_isi{:});

ModeHMM.(states{ss})(cc).prev_state = cat(2,ModeHMM.(states{ss})(cc).prev_state{:});
ModeHMM.(states{ss})(cc).next_state = cat(2,ModeHMM.(states{ss})(cc).next_state{:});
ModeHMM.(states{ss})(cc).state_spk = cat(2,ModeHMM.(states{ss})(cc).state_spk{:});

ModeHMM.(states{ss})(cc).prev_pMode = cat(1,ModeHMM.(states{ss})(cc).prev_pMode{:});
ModeHMM.(states{ss})(cc).next_pMode = cat(1,ModeHMM.(states{ss})(cc).next_pMode{:});

for sm = 1:numModes
    instate_both = ModeHMM.(states{ss})(cc).prev_state == sm & ModeHMM.(states{ss})(cc).next_state==sm;
    instate_either = ModeHMM.(states{ss})(cc).prev_state == sm | ModeHMM.(states{ss})(cc).next_state==sm;
    
    MeanReturn.(states{ss}).cells(cc).both(:,:,sm) = hist3([log10(ModeHMM.(states{ss})(cc).prev_isi(instate_both))',log10(ModeHMM.(states{ss})(cc).next_isi(instate_both))'],...
        {MeanReturn.logbins,MeanReturn.logbins});
    numspk = sum(sum(MeanReturn.(states{ss}).cells(cc).both(:,:,sm)));
    if numspk<spkthresh
        MeanReturn.(states{ss}).cells(cc).both(:,:,sm) = nan(length(MeanReturn.logbins));
    else
        MeanReturn.(states{ss}).cells(cc).both(:,:,sm) = MeanReturn.(states{ss}).cells(cc).both(:,:,sm)./numspk;
    end
    
    MeanReturn.(states{ss}).cells(cc).either(:,:,sm) = hist3([log10(ModeHMM.(states{ss})(cc).prev_isi(instate_either))',log10(ModeHMM.(states{ss})(cc).next_isi(instate_either))'],...
        {MeanReturn.logbins,MeanReturn.logbins});
    numspk = sum(sum(MeanReturn.(states{ss}).cells(cc).either(:,:,sm)));
    if numspk<spkthresh
        MeanReturn.(states{ss}).cells(cc).either(:,:,sm) = nan(length(MeanReturn.logbins));
    else
        MeanReturn.(states{ss}).cells(cc).either(:,:,sm) =  MeanReturn.(states{ss}).cells(cc).either(:,:,sm)./numspk;
    end
end

    %All mean
    MeanReturn.(states{ss}).cells(cc).allspikes(:,:) = hist3([log10(ModeHMM.(states{ss})(cc).prev_isi)',log10(ModeHMM.(states{ss})(cc).next_isi)'],...
        {MeanReturn.logbins,MeanReturn.logbins});
    numspk = sum(sum(MeanReturn.(states{ss}).cells(cc).allspikes(:,:)));
    if numspk<spkthresh
        MeanReturn.(states{ss}).cells(cc).allspikes(:,:) = nan(length(MeanReturn.logbins));
    else
        MeanReturn.(states{ss}).cells(cc).allspikes(:,:) =  MeanReturn.(states{ss}).cells(cc).allspikes(:,:)./numspk;
    end

end
MeanReturn.(states{ss}).mean = bz_CollapseStruct(MeanReturn.(states{ss}),4,'mean',true);

MeanTransition.(states{ss}) = bz_CollapseStruct(ModeHMM.(states{ss}),3,'mean',false);
MeanTransition.(states{ss}) = MeanTransition.(states{ss}).trans_mat;
end


%Mode intervals
modenames = {'AS1','AS2','AS3','AS4','AS5','AS6','GS'};
for ss = 1:2
for cc = 1:numcells
        ModeInts.(states{ss}).cells(cc) = bz_IDXtoINT(ModeHMM.(states{ss})(cc).next_state',...
            'statenames',modenames,'jumptol',1);
%         nexttemp = bz_IDXtoINT(ModeHMM.(states{ss})(cc).next_state',...
%             'statenames',modenames);
        
        %Could just use the spike index above....
        for sm = 1:numModes
            ModeInts_time.(states{ss}).cells(cc).(modenames{sm}) = ...
                [ModeHMM.(states{ss})(cc).state_spk(ModeInts.(states{ss}).cells(cc).([modenames{sm},'state'])(:,1))' ...
                ModeHMM.(states{ss})(cc).state_spk(ModeInts.(states{ss}).cells(cc).([modenames{sm},'state'])(:,2)+1)'];

            ModeInts.(states{ss}).chainlength(cc).(modenames{sm}) = diff(ModeInts.(states{ss}).cells(cc).([modenames{sm},'state']),1,2);
            [ModeInts.(states{ss}).chainhist(cc).(modenames{sm}),ModeInts.(states{ss}).chainbins.(modenames{sm})] = ...
                hist(ModeInts.(states{ss}).chainlength(cc).(modenames{sm}),[0:10]);
            ModeInts.(states{ss}).chainhist(cc).(modenames{sm}) = ModeInts.(states{ss}).chainhist(cc).(modenames{sm})./sum(ModeInts.(states{ss}).chainhist(cc).(modenames{sm}));
            
        end
end
end


%% Single cell example
GScolor = [0.6 0.4 0];
modecolors = crameri('bamako',6);
%modecolors = [modecolors;GScolor];
modecolors = {modecolors(1,:),modecolors(2,:),modecolors(3,:),modecolors(4,:),modecolors(5,:),modecolors(6,:),...
    GScolor};

lowthreshcolor = [0.95 0.95 0.95];
numrepeats = 3;
%excell = excells;
histcolors = [repmat([1 1 1],numrepeats,1);makeColorMap(lowthreshcolor,[0 0 0])];
NREMhistcolors = [repmat([1 1 1],numrepeats,1);makeColorMap(lowthreshcolor,[0 0 0.8])];
REMhistcolors = [repmat([1 1 1],numrepeats,1);makeColorMap(lowthreshcolor,[0.8 0 0])];
statecolormap = {histcolors,NREMhistcolors,REMhistcolors};

[ exwin ] = bz_RandomWindowInIntervals( SleepState.ints.WAKEstate,60,1 );

%exwin = [3308 3329];
exwin_long = [3322 3340];
exwin_short = [3335 3336];

%exwin_long = [3500 3600];

linethick = 10;

for cc = 1:numcells
%cc = 6;
checkUID = ModeHMM.WAKEstate(cc).UID;



figure

    subplot(2,3,1)
        bz_PlotISIDistModes(GammaFit.WAKEstate,checkUID)
        title(['Cell UID: ',num2str(checkUID)])

        
    GSints = ModeHMM.WAKEstate(cc).prev_state==numModes;
    clear modeHist
    for sm = 1:numModes
        inints = ModeHMM.WAKEstate(cc).prev_state==sm;
        modeHist(:,numModes+1-sm) = hist(log10(ModeHMM.WAKEstate(cc).prev_isi(inints)),ISIStats.ISIhist.logbins);
    end

    subplot(6,3,10)
        % hold on
        % histogram(log10(ModeHMM.WAKEstate(cc).prev_isi),'BinEdges', ISIStats.ISIhist.logbins,'facecolor','none')
        % histogram(log10(ModeHMM.WAKEstate(cc).prev_isi(GSints)),'facecolor',GScolor,'BinEdges', ISIStats.ISIhist.logbins)
        ba = bar(ISIStats.ISIhist.logbins,modeHist,'Stacked');
        LogScale('x',10,'exp',true,'nohalf',true)
        xlabel('ISI (s)');ylabel('# Intervals')
        for sm = 1:numModes
            ba(sm).FaceColor = modecolors{numModes+1-sm};
            if sm ==1
                ba(sm).EdgeColor = 'k';
            else
                ba(sm).EdgeColor = 'none';
            end
        end

    subplot(6,3,13)
        hold on
        plot(log10(ModeHMM.WAKEstate(cc).prev_isi(GSints)),ModeHMM.WAKEstate(cc).prev_pMode((GSints),numModes),...
            '.','color',GScolor,'markersize',1)

        plot(log10(ModeHMM.WAKEstate(cc).prev_isi(~GSints)),ModeHMM.WAKEstate(cc).prev_pMode((~GSints),numModes),...
            'k.','markersize',1)
        LogScale('x',10,'exp',true,'nohalf',true)
        xlabel('ISI (s)');
        ylabel('P(GS)')


    subplot(6,3,2:3)
        bz_PlotModeRaster(ModeHMM.WAKEstate,ModeInts_time.WAKEstate.cells,cc,modecolors,exwin,'linethick',linethick)

    subplot(6,3,[5:6])
        hold on
        linewidths = [1 1 1 1 1 1 2];
        for mm = numModes:-1:1
        bothspikes = [ModeHMM.WAKEstate(cc).state_spk;ModeHMM.WAKEstate(cc).state_spk];
        bothpmodes = [ModeHMM.WAKEstate(cc).prev_pMode(:,mm) ModeHMM.WAKEstate(cc).next_pMode(:,mm)]';
        %plot(ModeHMM.WAKEstate(checkUID).state_spk,ModeHMM.WAKEstate(checkUID).prev_pMode(:,6))
        plot(bothspikes(:),bothpmodes(:),'color',modecolors{mm},'linewidth',linewidths(mm))
        end
        xlim(exwin)
        %bz_ScaleBar('s')
        set(gca,'xtick',[])
        ylabel('p(state)')

    subplot(6,7,5*7+1)
        plot(log10(ModeHMM.WAKEstate(cc).prev_isi),log10(ModeHMM.WAKEstate(cc).next_isi),'k.','markersize',0.5)
        xlim([-3 2]);ylim([-3 2])
        title('All Spikes')
        xlabel('ISI_n');ylabel('ISI_n_+_1')
        LogScale('xy',10,'exp',true)


    for sm = 1:numModes
        moderate = round(10.^ModeHMM.WAKEstate(cc).logrates(sm),0);
        modeCV = round(ModeHMM.WAKEstate(cc).cvs(sm),1);
        %if ss==6
            instate_both = ModeHMM.WAKEstate(cc).prev_state == sm & ModeHMM.WAKEstate(cc).next_state==sm;
        %else
            instate_either = ModeHMM.WAKEstate(cc).prev_state == sm | ModeHMM.WAKEstate(cc).next_state==sm;
        %end
        subplot(6,numModes+1,5*(numModes+1)+(sm)+1)
            plot(log10(ModeHMM.WAKEstate(cc).prev_isi(instate_either)),log10(ModeHMM.WAKEstate(cc).next_isi(instate_either)),...
                '.','color',[0.5 0.5 0.5],'markersize',0.1)
            hold on
            plot(log10(ModeHMM.WAKEstate(cc).prev_isi(instate_both)),log10(ModeHMM.WAKEstate(cc).next_isi(instate_both)),...
                'k.','markersize',0.1)
             set(gca,'yticklabel',[]);set(gca,'xticklabel',[])
            xlim([-3 2]);ylim([-3 2])
            if sm == numModes
                title('GS','Color',modecolors{sm})
            else
                title(['AS',num2str(sm),' (',num2str(moderate),'Hz, CV: ',num2str(modeCV),')'],'Color',modecolors{sm})
            end
        
        subplot(3,3,6)
        hold on
            plot(ModeInts.WAKEstate.chainbins.(modenames{sm})+1,ModeInts.WAKEstate.chainhist(cc).(modenames{sm}),'color',modecolors{sm})
            box off
            
            if sm == numModes
                xlabel('Chain Length (# ISIs)')
            end
    end

    subplot(3,3,5)
        imagesc(ModeHMM.WAKEstate(cc).trans_mat(1:numModes,1:numModes))
        xlabel('To State');ylabel('From State')
        ColorbarWithAxis([0 0.6],'P(Transition)','inclusive',{'','>'})

    NiceSave(['CellExample_',num2str(cc)],[figfolder,'/HMMAllCells'],baseName)
%close all
end



%%
figure
for ss = 1:2
    for sm = 1:numModes
        subplot(numModes,2,(sm-1)*2+ss)
            bar(ModeInts_all.(states{ss}).chainbins.(modenames{sm})+1,ModeInts_all.(states{ss}).chainhist.(modenames{sm}))
            box off
            
            if ss == 1
                ylabel(modenames{sm})
            end
            if sm == 1
                title(states{ss})
            end
            if sm == numModes
                xlabel('Chain Length (# ISIs)')
            end
    end
end
NiceSave('ModeChainLengths',figfolder,baseName)


%% Check Mode RATE/CV

    ALLcellWAKEmodes = bz_CollapseStruct(ModeHMM.WAKEstate,1);

[~,ordermatch_ISIGamma] = ismember([ModeHMM.WAKEstate(:).UID],GammaFit.WAKEstate.cellstats.UID);
    
    %%
    figure
    subplot(2,2,1)
        imagesc(ALLcellWAKEmodes.logrates)
        xlabel('Mode Number');ylabel('Cell')
        colorbar
        LogScale('c',10)
        title('Mode Rate')
    subplot(2,2,3)
        imagesc(ALLcellWAKEmodes.cvs)
        xlabel('Mode Number');ylabel('Cell')
        colorbar
        %LogScale('c',10)
        title('Mode CV')
    subplot(2,2,2)
        scatter(GammaFit.WAKEstate.sharedfit.GSlogrates(ordermatch_ISIGamma),ALLcellWAKEmodes.logrates(:,numModes),10,...
            ALLcellWAKEmodes.logrates(:,2))
        %plot(GammaFit.WAKEstate.sharedfit.GSlogrates(ordermatch_ISIGamma),ALLcellWAKEmodes.logrates(:,6),'.')
        xlabel('Gamma Fit GS Rate');ylabel('HMM GS Rate')

%% Get rates for sorting
[~,ordermatch] = ismember([ModeHMM.WAKEstate(:).UID],ISIStats.UID);
cellrates.WAKEstate = ISIStats.summstats.WAKEstate.meanrate(ordermatch);
%% Many-neuron figure

plotnumcells  = 20;
plotrates = cellrates.WAKEstate(1:plotnumcells);
[~,cellorder] = sort(plotrates);
%[ exwin ] = bz_RandomWindowInIntervals( SleepState.ints.WAKEstate,60,1 )
exwin = [4980 5025]
linethick = 4
figure
subplot(2,1,1)
bz_PlotModeRaster(ModeHMM.WAKEstate,ModeInts_time.WAKEstate.cells,cellorder,modecolors,exwin,'linethick',linethick)
ylabel('Cell (30), sort by rate')

subplot(6,7,5*7+1)
    imagesc(MeanReturn.WAKEstate.mean.cells.allspikes(:,:))
    colormap(gca,statecolormap{1})
    axis xy
    set(gca,'yticklabel',[]);set(gca,'xticklabel',[])

for sm = 1:numModes
    subplot(6,numModes+1,4*(numModes+1)+(sm)+1)
        imagesc(MeanReturn.WAKEstate.mean.cells.both(:,:,sm))
            colormap(gca,statecolormap{1})
        axis xy
        set(gca,'yticklabel',[]);set(gca,'xticklabel',[])
        
    subplot(6,numModes+1,5*(numModes+1)+(sm)+1)
        imagesc(MeanReturn.WAKEstate.mean.cells.either(:,:,sm))
            colormap(gca,statecolormap{1})
        axis xy
        set(gca,'yticklabel',[]);set(gca,'xticklabel',[])
end
NiceSave('MultiCellExample',figfolder,baseName)

%%

%[ exwin ] = bz_RandomWindowInIntervals( SleepState.ints.WAKEstate,60,1 )
exwin = [5005 5012]
linethick = 4
figure
subplot(2,1,1)
bz_PlotModeRaster(ModeHMM.WAKEstate,ModeInts_time.WAKEstate.cells,cellorder,modecolors,exwin,'linethick',linethick)
ylabel('Cell (30), sort by rate')

exwin = [5009.5 5011]
subplot(2,1,2)
bz_PlotModeRaster(ModeHMM.WAKEstate,ModeInts_time.WAKEstate.cells,cellorder,modecolors,exwin,'linethick',linethick+2)
ylabel('Cell (30), sort by rate')
ylim([15 31])
NiceSave('MultiCellExample_Zoom',figfolder,baseName)
%ylim(
%%

[ exwin ] = bz_RandomWindowInIntervals( SleepState.ints.WAKEstate,60,1 )

%spiketimes


plotnumcells  = 20;
plotrates = cellrates.WAKEstate(1:plotnumcells);
[~,cellorder] = sort(plotrates);
%rates = 
figure
subplot(2,1,1)
hold on
for cc = 1:plotnumcells
    whichcell = cellorder(cc);
    plotints =  structfun(@(modeints) RestrictInts(modeints,exwin,'inclusive',true),ModeInts_time.WAKEstate.cells(whichcell),'UniformOutput',false);
StateScorePlot( plotints,modecolors,'y',cc,'LineWidth',5)

for sm = 1:6
    %if ss==6
        instate_both = ModeHMM.WAKEstate(whichcell).prev_state == sm & ModeHMM.WAKEstate(whichcell).next_state==sm & ...
            InIntervals(ModeHMM.WAKEstate(whichcell).state_spk',exwin)';
    %else
        instate_either = (ModeHMM.WAKEstate(whichcell).prev_state == sm | ModeHMM.WAKEstate(whichcell).next_state==sm) & ...
            InIntervals(ModeHMM.WAKEstate(whichcell).state_spk',exwin)';
    %end
    plot([ModeHMM.WAKEstate(whichcell).state_spk(instate_either);ModeHMM.WAKEstate(whichcell).state_spk(instate_either)],...
        cc+[zeros(size(ModeHMM.WAKEstate(whichcell).state_spk(instate_either)))-0.4;0.4+zeros(size(ModeHMM.WAKEstate(whichcell).state_spk(instate_either)))],...
        'color',[0.5 0.5 0.5],'linewidth',0.5)
    plot([ModeHMM.WAKEstate(whichcell).state_spk(instate_both);ModeHMM.WAKEstate(whichcell).state_spk(instate_both)],...
        cc+[zeros(size(ModeHMM.WAKEstate(whichcell).state_spk(instate_both)))-0.4;0.4+zeros(size(ModeHMM.WAKEstate(whichcell).state_spk(instate_both)))],...
        'color','k','linewidth',1)
    
end
end
xlim(exwin)
ylim([0 plotnumcells+1])
box off
bz_ScaleBar('s')



% subplot(3,3,1)
% plot(log10(ModeHMM.WAKEstate(cc).prev_isi),log10(ModeHMM.WAKEstate(cc).next_isi),'k.','markersize',1)
% xlim([-3 2]);ylim([-3 2])
% title('All Spikes')
% xlabel('ISI_n');ylabel('ISI_n_+_1')
% LogScale('xy',10,'exp',true)
% 
% subplot(3,3,7)
% imagesc(MeanTransition.(states{ss})([2,4,1,5,3,6],[2,4,1,5,3,6],:))
% xlabel('To State');ylabel('From State')
% ColorbarWithAxis([0 0.6],'P(Transition)','inclusive',{'','>'})
cc = 1;
for sm = 1:numModes
    %if ss==6
        instate_both = ModeHMM.WAKEstate(cc).prev_state == sm & ModeHMM.WAKEstate(cc).next_state==sm;
    %else
        instate_either = ModeHMM.WAKEstate(cc).prev_state == sm | ModeHMM.WAKEstate(cc).next_state==sm;
    %end
subplot(6,numModes,3*numModes+(sm))
    plot(log10(ModeHMM.WAKEstate(cc).prev_isi(instate_either)),log10(ModeHMM.WAKEstate(cc).next_isi(instate_either)),...
        '.','color',[0.5 0.5 0.5],'markersize',0.5)
    hold on
    plot(log10(ModeHMM.WAKEstate(cc).prev_isi(instate_both)),log10(ModeHMM.WAKEstate(cc).next_isi(instate_both)),...
        'k.','markersize',0.5)
     set(gca,'yticklabel',[]);set(gca,'xticklabel',[])
    xlim([-3 2]);ylim([-3 2])
    if sm == numModes
        title('GS','Color',modecolors{sm})
    else
        title(['AS',num2str(sm)],'Color',modecolors{sm})
    end
end

for sm = 1:numModes
    subplot(6,numModes,4*numModes+(sm))
        imagesc(MeanReturn.WAKEstate.mean.cells.both(:,:,sm))
        axis xy
        set(gca,'yticklabel',[]);set(gca,'xticklabel',[])
        
    subplot(6,numModes,5*numModes+(sm))
        imagesc(MeanReturn.WAKEstate.mean.cells.either(:,:,sm))
        axis xy
        set(gca,'yticklabel',[]);set(gca,'xticklabel',[])

end

NiceSave('ModeRaster',figfolder,baseName)


%%
position = bz_LoadBehavior( basePath,'position' ); %For finding a good time window
[ exwin ] = bz_RandomWindowInIntervals( position.Epochs.MazeEpoch,30,1 )
%%
%spiketimes


% plotnumcells  = 1;
% plotrates = cellrates.WAKEstate(1:plotnumcells);
% [~,cellorder] = sort(plotrates);
% %rates = 
% figure
% subplot(2,1,1)
% hold on
% for cc = 1:plotnumcells
%     whichcell = cellorder(cc);
%     plotints =  structfun(@(modeints) RestrictInts(modeints,exwin,'inclusive',true),ModeInts_time.WAKEstate.cells(whichcell),'UniformOutput',false);
% StateScorePlot( plotints,modecolors,'y',cc,'LineWidth',5)
% 
% for sm = 1:6
%     %if ss==6
%         instate_both = ModeHMM.WAKEstate(whichcell).prev_state == sm & ModeHMM.WAKEstate(whichcell).next_state==sm & ...
%             InIntervals(ModeHMM.WAKEstate(whichcell).state_spk',exwin)';
%     %else
%         instate_either = (ModeHMM.WAKEstate(whichcell).prev_state == sm | ModeHMM.WAKEstate(whichcell).next_state==sm) & ...
%             InIntervals(ModeHMM.WAKEstate(whichcell).state_spk',exwin)';
%     %end
%     plot([ModeHMM.WAKEstate(whichcell).state_spk(instate_either);ModeHMM.WAKEstate(whichcell).state_spk(instate_either)],...
%         cc+[zeros(size(ModeHMM.WAKEstate(whichcell).state_spk(instate_either)))-0.4;0.4+zeros(size(ModeHMM.WAKEstate(whichcell).state_spk(instate_either)))],...
%         'color',[0.5 0.5 0.5],'linewidth',0.5)
%     plot([ModeHMM.WAKEstate(whichcell).state_spk(instate_both);ModeHMM.WAKEstate(whichcell).state_spk(instate_both)],...
%         cc+[zeros(size(ModeHMM.WAKEstate(whichcell).state_spk(instate_both)))-0.4;0.4+zeros(size(ModeHMM.WAKEstate(whichcell).state_spk(instate_both)))],...
%         'color','k','linewidth',1)
%     
% end
% end
% xlim(exwin)
% ylim([0 plotnumcells+1])
% box off
% bz_ScaleBar('s')



% subplot(3,3,1)
% plot(log10(ModeHMM.WAKEstate(cc).prev_isi),log10(ModeHMM.WAKEstate(cc).next_isi),'k.','markersize',1)
% xlim([-3 2]);ylim([-3 2])
% title('All Spikes')
% xlabel('ISI_n');ylabel('ISI_n_+_1')
% LogScale('xy',10,'exp',true)
% 
% subplot(3,3,7)
% imagesc(MeanTransition.(states{ss})([2,4,1,5,3,6],[2,4,1,5,3,6],:))
% xlabel('To State');ylabel('From State')
% ColorbarWithAxis([0 0.6],'P(Transition)','inclusive',{'','>'})
%cc = 8;
%for cc = 1:numcells
figure
subplot(4,1,1)
bz_PlotModeRaster(ModeHMM.WAKEstate,ModeInts_time.WAKEstate.cells,cc,modecolors,exwin,'linethick',linethick)


for sm = 1:numModes
    moderate = round(10.^ModeHMM.WAKEstate(cc).logrates(sm),0);
    modeCV = round(ModeHMM.WAKEstate(cc).cvs(sm),1);
    %if ss==6
        instate_both = ModeHMM.WAKEstate(cc).prev_state == sm & ModeHMM.WAKEstate(cc).next_state==sm;
    %else
        instate_either = ModeHMM.WAKEstate(cc).prev_state == sm | ModeHMM.WAKEstate(cc).next_state==sm;
    %end
subplot(6,numModes,3*numModes+(sm))
    plot(log10(ModeHMM.WAKEstate(cc).prev_isi(instate_either)),log10(ModeHMM.WAKEstate(cc).next_isi(instate_either)),...
        '.','color',[0.5 0.5 0.5],'markersize',0.5)
    hold on
    plot(log10(ModeHMM.WAKEstate(cc).prev_isi(instate_both)),log10(ModeHMM.WAKEstate(cc).next_isi(instate_both)),...
        'k.','markersize',0.5)
     set(gca,'yticklabel',[]);set(gca,'xticklabel',[])
    xlim([-3 2]);ylim([-3 2])
    if sm == numModes
        title('GS','Color',modecolors{sm})
    else
        title(['AS',num2str(sm),' (',num2str(moderate),'Hz, CV: ',num2str(modeCV),')'],'Color',modecolors{sm})
    end
end

for sm = 1:numModes
    subplot(6,numModes,4*numModes+(sm))
        imagesc(MeanReturn.WAKEstate.mean.cells.both(:,:,sm))
        axis xy
        set(gca,'yticklabel',[]);set(gca,'xticklabel',[])
        
    subplot(6,numModes,5*numModes+(sm))
        imagesc(MeanReturn.WAKEstate.mean.cells.either(:,:,sm))
        axis xy
        set(gca,'yticklabel',[]);set(gca,'xticklabel',[])

end

NiceSave(['ModeRaster_UID',num2str(ModeHMM.WAKEstate(cc).UID)],figfolder,baseName)
close all
%end

end
