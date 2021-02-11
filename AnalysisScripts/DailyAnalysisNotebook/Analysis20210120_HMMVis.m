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
%ISIStats = bz_LoadCellinfo(basePath,'ISIStats');
states = fieldnames(SleepState.ints);
%states{4} = 'ALL';
%SleepState.ints.ALL = [0 Inf];
%statecolors = {'k','b','r',[0.6 0.6 0.6]};

% try
%     celltypes = CellClass.celltypes;
% catch
%     celltypes = unique(CellClass.label);
% end
% cellcolor = {'k','r'};


%%
load([basePath,'/GammaProcessed/hmm_out.mat'])

%%
ModeHMM.WAKEstate = WAKEall;
ModeHMM.NREMstate = NREMall;

numcells = length(WAKEall);
spkthresh = 50;
MeanReturn.logbins = linspace(-3,2,50);
%get next ISI (nan for last one in the state)
%Cat all the cells
for ss = 1:2
for cc = 1:numcells

ModeHMM.(states{ss})(cc).next_isi = cellfun(@(prevISI) [prevISI(2:end) nan],ModeHMM.(states{ss})(cc).state_isi,'UniformOutput',false);
ModeHMM.(states{ss})(cc).next_state = cellfun(@(prevState) [prevState(2:end) nan],ModeHMM.(states{ss})(cc).decoded_mode,'UniformOutput',false);

ModeHMM.(states{ss})(cc).prev_isi = cellfun(@(prevISI) prevISI(1:end),ModeHMM.(states{ss})(cc).state_isi,'UniformOutput',false);
ModeHMM.(states{ss})(cc).prev_state = cellfun(@(prevState) prevState(1:end),ModeHMM.(states{ss})(cc).decoded_mode,'UniformOutput',false);

ModeHMM.(states{ss})(cc).prev_isi = cat(2,ModeHMM.(states{ss})(cc).prev_isi{:});
ModeHMM.(states{ss})(cc).next_isi = cat(2,ModeHMM.(states{ss})(cc).next_isi{:});

ModeHMM.(states{ss})(cc).prev_state = cat(2,ModeHMM.(states{ss})(cc).prev_state{:});
ModeHMM.(states{ss})(cc).next_state = cat(2,ModeHMM.(states{ss})(cc).next_state{:});
ModeHMM.(states{ss})(cc).state_spk = cat(2,ModeHMM.(states{ss})(cc).state_spk{:});

for sm = 1:6
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

end
MeanReturn.(states{ss}).mean = bz_CollapseStruct(MeanReturn.(states{ss}),4,'mean',true);

MeanTransition.(states{ss}) = bz_CollapseStruct(ModeHMM.(states{ss}),3,'mean',false);
MeanTransition.(states{ss}) = MeanTransition.(states{ss}).trans_mat;
end

%%

%%
ss =1
cc = 7;
figure
subplot(3,3,1)
plot(log10(ModeHMM.WAKEstate(cc).prev_isi),log10(ModeHMM.WAKEstate(cc).next_isi),'k.','markersize',1)
xlim([-3 2]);ylim([-3 2])
title('All Spikes')
xlabel('ISI_n');ylabel('ISI_n_+_1')
LogScale('xy',10,'exp',true)

subplot(3,3,7)
imagesc(MeanTransition.(states{ss})([2,4,1,5,3,6],[2,4,1,5,3,6],:))
xlabel('To State');ylabel('From State')
ColorbarWithAxis([0 0.6],'P(Transition)','inclusive',{'','>'})

for sm = 1:6
    %if ss==6
        instate_both = ModeHMM.WAKEstate(cc).prev_state == sm & ModeHMM.WAKEstate(cc).next_state==sm;
    %else
        instate_either = ModeHMM.WAKEstate(cc).prev_state == sm | ModeHMM.WAKEstate(cc).next_state==sm;
    %end
subplot(6,6,(sm-1)*6+3)
    plot(log10(ModeHMM.WAKEstate(cc).prev_isi(instate_either)),log10(ModeHMM.WAKEstate(cc).next_isi(instate_either)),...
        '.','color',[0.5 0.5 0.5],'markersize',0.5)
    hold on
    plot(log10(ModeHMM.WAKEstate(cc).prev_isi(instate_both)),log10(ModeHMM.WAKEstate(cc).next_isi(instate_both)),...
        'k.','markersize',0.5)
     set(gca,'yticklabel',[]);set(gca,'xticklabel',[])
    xlim([-3 2]);ylim([-3 2])
    if sm == 6
        title('GS')
    else
        title(['AS',num2str(sm)])
    end
end

for sm = 1:6
    subplot(6,6,(sm-1)*6+4)
        imagesc(MeanReturn.(states{ss}).mean.cells.both(:,:,sm))
        axis xy
        set(gca,'yticklabel',[]);set(gca,'xticklabel',[])
        
    subplot(6,6,(sm-1)*6+5)
        imagesc(MeanReturn.(states{ss}).mean.cells.either(:,:,sm))
        axis xy
        set(gca,'yticklabel',[]);set(gca,'xticklabel',[])

end
NiceSave('ModeReturnmaps',figfolder,baseName)



%% Chain Length
modenames = {'AS1','AS2','AS3','AS4','AS5','GS'};
for ss = 1:2
for cc = 1:numcells
        ModeInts.(states{ss}).cells(cc) = bz_IDXtoINT(ModeHMM.(states{ss})(cc).next_state',...
            'statenames',modenames,'jumptol',1);
%         nexttemp = bz_IDXtoINT(ModeHMM.(states{ss})(cc).next_state',...
%             'statenames',modenames);
        
        %Could just use the spike index above....
        for sm = 1:6
        ModeInts_time.(states{ss}).cells(cc).(modenames{sm}) = ...
            [ModeHMM.(states{ss})(cc).state_spk(ModeInts.(states{ss}).cells(cc).([modenames{sm},'state'])(:,1))' ...
            ModeHMM.(states{ss})(cc).state_spk(ModeInts.(states{ss}).cells(cc).([modenames{sm},'state'])(:,2)+1)'];
        end
end
end


%%
ModeInts_all = bz_CollapseStruct(ModeInts,1,'justcat',true);

%%
for ss = 1:2
    for sm = 1:6
        ModeInts_all.(states{ss}).chainlength.(modenames{sm}) = diff(ModeInts_all.(states{ss}).cells.([modenames{sm},'state']),1,2);
        [ModeInts_all.(states{ss}).chainhist.(modenames{sm}),ModeInts_all.(states{ss}).chainbins.(modenames{sm})] = ...
            hist(ModeInts_all.(states{ss}).chainlength.(modenames{sm}),[0:10]);
        ModeInts_all.(states{ss}).chainhist.(modenames{sm}) = ModeInts_all.(states{ss}).chainhist.(modenames{sm})./sum(ModeInts_all.(states{ss}).chainhist.(modenames{sm}));
    end
end

%%
figure
for ss = 1:2
    for sm = 1:6
        subplot(6,2,(sm-1)*2+ss)
            bar(ModeInts_all.(states{ss}).chainbins.(modenames{sm})+1,ModeInts_all.(states{ss}).chainhist.(modenames{sm}))
            box off
            
            if ss == 1
                ylabel(modenames{sm})
            end
            if sm == 1
                title(states{ss})
            end
            if sm == 6
                xlabel('Chain Length (# ISIs)')
            end
    end
end
NiceSave('ModeChainLengths',figfolder,baseName)


%%
GScolor = [0.6 0.4 0];
modecolors = crameri('bamako',5);
%modecolors = [modecolors;GScolor];
modecolors = {modecolors(1,:),modecolors(2,:),modecolors(3,:),modecolors(4,:),modecolors(5,:),...
    GScolor};

cc = 10;

[ exwin ] = bz_RandomWindowInIntervals( SleepState.ints.WAKEstate,100,1 )

figure
subplot(5,1,1)
hold on
for cc = 1:21
StateScorePlot( ModeInts_time.WAKEstate.cells(cc),modecolors,'y',cc )

for sm = 1:6
    %if ss==6
        instate_both = ModeHMM.WAKEstate(cc).prev_state == sm & ModeHMM.WAKEstate(cc).next_state==sm;
    %else
        instate_either = ModeHMM.WAKEstate(cc).prev_state == sm | ModeHMM.WAKEstate(cc).next_state==sm;
    %end
    plot(ModeHMM.WAKEstate(cc).state_spk(instate_either),cc.*ones(size(ModeHMM.WAKEstate(cc).state_spk(instate_either))),...
        '+','color',[0.5 0.5 0.5])
    plot(ModeHMM.WAKEstate(cc).state_spk(instate_both),cc.*ones(size(ModeHMM.WAKEstate(cc).state_spk(instate_both))),'k+')
    
end
end
xlim(exwin)
ylim([0 20])
% for sm = 1:5
%     plot(ModeInts_time.(states{ss}).cells(cc).([modenames{sm},'state']),ones(size(ModeInts_time.(states{ss}).cells(cc).([modenames{sm},'state']))),...
%         'color',modecolors(sm,:),'linewidth',2)
% end

%StateScorePlot( ModeInts_time.(states{ss}).cells(cc),modecolors )
NiceSave('ModeRaster',figfolder,baseName)
