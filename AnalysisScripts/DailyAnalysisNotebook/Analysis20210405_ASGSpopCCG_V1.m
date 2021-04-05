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
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/Cicero_09102014';
basePath = '/Users/dl2820/Dropbox/Research/Datasets/YMV11_171208';
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

% try
%     celltypes = CellClass.celltypes;
% catch
%     celltypes = unique(CellClass.label);
% end
% cellcolor = {'k','r'};


%%
load([basePath,'/GammaProcessed1/hmm_out.mat'])
GammaFit = bz_LoadCellinfo(basePath,'GammaFit');
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


    
end

end

%%
ignorepairs = false(numcells.*7);
ss = 2;
for sm = 1:7
    CellClass.celltypes{sm} = ['Mode',num2str(sm)];
    CellClass.(CellClass.celltypes{sm}) = false(numcells.*7);
    for cc = 1:numcells
        cellidx = cc + (sm-1).*numcells;
        
        if sm == 6
            inmode = ModeHMM.(states{ss})(cc).next_state == sm & ModeHMM.(states{ss})(cc).prev_state==sm;
        elseif sm == 7 %All AS spikes
             inmode = ~(ModeHMM.(states{ss})(cc).next_state == 6 & ModeHMM.(states{ss})(cc).prev_state==6);
        else
            inmode = ModeHMM.(states{ss})(cc).next_state == sm | ModeHMM.(states{ss})(cc).prev_state==sm;
         end
        spikes_modes.times{cellidx} = ModeHMM.(states{ss})(cc).state_spk(inmode)';
        CellClass.label{cellidx} = ['Mode',num2str(sm)];
        CellClass.(CellClass.celltypes{sm})(cellidx) = true;
        for sm2 = 1:7
            cellidx2 = cc + (sm2-1).*numcells;
            ignorepairs(cellidx,cellidx2)=true;
        end
    end
    
end

%%
minSpikes = 150;
duration = 0.5;
[popCCG.(states{ss})] = PopCCG(spikes_modes,'showfig',true,'cellclass',CellClass.label,...
    'classnames',CellClass.celltypes,'ignorepairs',ignorepairs,'minspikes',minSpikes,...
    'duration',duration);%,...
    %'sort',ISIStats.sorts.(states{ss}).ratebyclass);

%%
for sm1 = 1:6
    for sm2 = 1:6
        temp(:,sm1,sm2) = popCCG.NREMstate.cellsmean.(CellClass.celltypes{sm1})(:,sm2);
    end
end

GSGS = temp(:,6,6); %GS-GS
GSAS = mean(temp(:,6,1:5),3); %AS spikes referenced to GS
ASGS = mean(temp(:,1:5,6),2); %AS spikes referenced to GS

temp2=[];
temp3=[];
for sm1 = 1:5
    temp2 = [temp2 temp(:,sm1,sm1)];
    for sm2 = 1:5
        if sm1~=sm2

            temp3 = [temp3 temp(:,sm1,sm2)];
        end
    end
end
ASAS_same = mean(temp2,2);
ASAS_diff = mean(temp3,2);

%%
GScolor = [0.6 0.4 0];

figure
subplot(3,3,1)
%plot(popCCG.WAKEstate.t_ccg,GSGS)
plot(popCCG.NREMstate.t_ccg,popCCG.NREMstate.pop.Mode6(:,6),'linewidth',1,'color',GScolor)
hold on
box off
plot(popCCG.NREMstate.t_ccg,popCCG.NREMstate.pop.Mode6(:,7),'linewidth',1,'color',GScolor.*0.7)
%xlabel('Time relative to GS spike')
%plot(popCCG.WAKEstate.t_ccg,GSAS)
legend('GS','AS')
ylabel('GS Rate')
plot([0 0],[0 1],'--','color',GScolor.*0.7)
xlabel('t - rel to ref spike (s)')

subplot(3,3,2)
plot(popCCG.NREMstate.t_ccg,popCCG.NREMstate.pop.Mode7(:,6),'linewidth',1,'color',GScolor.*0.7)
hold on
box off
plot(popCCG.NREMstate.t_ccg,popCCG.NREMstate.pop.Mode7(:,7),'k','linewidth',1)
legend('GS','AS')
ylabel('AS Rate')
plot([0 0],[0 1],'--','color',GScolor.*0.7)
xlabel('t - rel to ref spike (s)')


subplot(3,3,4)
plot(popCCG.NREMstate.t_ccg,GSGS,'linewidth',1,'color',GScolor)
%plot(popCCG.WAKEstate.t_ccg,popCCG.WAKEstate.pop.Mode6(:,6))
hold on
box off
%plot(popCCG.WAKEstate.t_ccg,popCCG.WAKEstate.pop.Mode6(:,7))
plot(popCCG.NREMstate.t_ccg,GSAS,'linewidth',1,'color',GScolor.*0.7)
legend('GS','AS')
ylabel('GS Rate')
plot([0 0],[0 1],'--','color',GScolor.*0.7)
xlabel('t - rel to ref spike (s)')

subplot(3,3,5)
plot(popCCG.NREMstate.t_ccg,ASGS,'linewidth',1,'color',GScolor.*0.7)
hold on
box off
plot(popCCG.NREMstate.t_ccg,ASAS_same,'linewidth',1,'color','k')
plot(popCCG.NREMstate.t_ccg,ASAS_diff,'linewidth',1,'color',[0.5 0.5 0.5])
legend('GS','Same','Diff')
ylabel('AS Rate')
plot([0 0],[0 0.15],'--','color',GScolor.*0.7)
xlabel('t - rel to ref spike (s)')


subplot(3,3,7)
plot(popCCG.NREMstate.t_ccg,popCCG.NREMstate.pop.Mode6(:,6),'linewidth',1,'color',GScolor)
hold on
box off
plot(popCCG.NREMstate.t_ccg,popCCG.NREMstate.pop.Mode7(:,6),'k','linewidth',1)
plot([0 0],[0 1],'--','color',GScolor)
xlabel('t - rel to GS spike (s)')
ylabel('Pop Rate (spk/s)')

subplot(3,3,8)
plot(popCCG.NREMstate.t_ccg,popCCG.NREMstate.pop.Mode6(:,7),'linewidth',1,'color',GScolor)
hold on
box off
plot(popCCG.NREMstate.t_ccg,popCCG.NREMstate.pop.Mode7(:,7),'k','linewidth',1)
plot([0 0],[0 1],'--','color','k')
xlabel('t - rel to AS spike (s)')
ylabel('Pop Rate (spk/s)')

NiceSave('GSASPopGGS',figfolder,baseName);
%%
figure
subplot(3,3,1)
%plot(popCCG.WAKEstate.t_ccg,GSGS)
plot(popCCG.NREMstate.t_ccg,popCCG.NREMstate.pop.Mode6(:,6),'linewidth',0.5,'color',GScolor)
hold on
box off
plot(popCCG.NREMstate.t_ccg,popCCG.NREMstate.pop.Mode6(:,7),'linewidth',0.5,'color',GScolor.*0.7)
%xlabel('Time relative to GS spike')
%plot(popCCG.WAKEstate.t_ccg,GSAS)
%legend('GS','AS')
ylabel('GS Rate')
axis tight
plot([0 0],ylim,'--','color',GScolor.*0.7)
xlabel('t - rel to ref spike (s)')
xlim(0.03.*[-1 1])

subplot(3,3,2)
plot(popCCG.NREMstate.t_ccg,popCCG.NREMstate.pop.Mode7(:,6),'linewidth',0.5,'color',GScolor.*0.7)
hold on
box off
plot(popCCG.NREMstate.t_ccg,popCCG.NREMstate.pop.Mode7(:,7),'k','linewidth',0.5)
%legend('GS','AS')
ylabel('AS Rate')
axis tight
plot([0 0],ylim,'--','color',GScolor.*0.7)
xlabel('t - rel to ref spike (s)')
xlim(0.03.*[-1 1])


subplot(3,3,4)
plot(popCCG.NREMstate.t_ccg,GSGS,'linewidth',0.5,'color',GScolor)
%plot(popCCG.WAKEstate.t_ccg,popCCG.WAKEstate.pop.Mode6(:,6))
hold on
box off
%plot(popCCG.WAKEstate.t_ccg,popCCG.WAKEstate.pop.Mode6(:,7))
plot(popCCG.NREMstate.t_ccg,GSAS,'linewidth',0.5,'color',GScolor.*0.7)
%legend('GS','AS')
ylabel('GS Rate')
axis tight
plot([0 0],ylim,'--','color',GScolor.*0.7)
xlabel('t - rel to ref spike (s)')
xlim(0.03.*[-1 1])

subplot(3,3,5)
plot(popCCG.NREMstate.t_ccg,ASGS,'linewidth',0.5,'color',GScolor.*0.7)
hold on
box off
plot(popCCG.NREMstate.t_ccg,ASAS_same,'linewidth',0.5,'color','k')
plot(popCCG.NREMstate.t_ccg,ASAS_diff,'linewidth',0.5,'color',[0.5 0.5 0.5])
% legend('GS','Same','Diff')
ylabel('AS Rate')
axis tight
plot([0 0],ylim,'--','color',GScolor.*0.7)
xlabel('t - rel to ref spike (s)')
xlim(0.03.*[-1 1])


subplot(3,3,7)
plot(popCCG.NREMstate.t_ccg,popCCG.NREMstate.pop.Mode6(:,6),'linewidth',0.5,'color',GScolor)
hold on
box off
plot(popCCG.NREMstate.t_ccg,popCCG.NREMstate.pop.Mode7(:,6),'k','linewidth',0.5)
axis tight
plot([0 0],ylim,'--','color',GScolor)
xlabel('t - rel to GS spike (s)')
ylabel('Pop Rate (spk/s)')
xlim(0.03.*[-1 1])

subplot(3,3,8)
plot(popCCG.NREMstate.t_ccg,popCCG.NREMstate.pop.Mode6(:,7),'linewidth',0.5,'color',GScolor)
hold on
box off
plot(popCCG.NREMstate.t_ccg,popCCG.NREMstate.pop.Mode7(:,7),'k','linewidth',0.5)
axis tight
plot([0 0],ylim,'--','color','k')
xlabel('t - rel to AS spike (s)')
ylabel('Pop Rate (spk/s)')
xlim(0.03.*[-1 1])


NiceSave('GSASPopGGS_zoom',figfolder,baseName);
%%
%%
modenames = {'AS1','AS2','AS3','AS4','AS5','GS'};

morder = [2 4 1 5 3 6];
GScolor = [0.6 0.4 0];
modecolors = crameri('bamako',5);
%modecolors = [modecolors;GScolor];
modecolors = {modecolors(1,:),modecolors(2,:),modecolors(3,:),modecolors(4,:),modecolors(5,:),...
    GScolor};

figure
for sm1 = 1:6
    for sm2 = 1:6
        subplot(6,6,sm2+(sm1-1).*6)
            ylim([min(popCCG.NREMstate.cellsmean.(CellClass.celltypes{morder(sm1)})(:)) max(popCCG.NREMstate.cellsmean.(CellClass.celltypes{morder(sm1)})(:))])
            plot([0 0],ylim,'--','color',modecolors{morder(sm2)},'linewidth',1)
            hold on
            plot(popCCG.NREMstate.t_ccg,popCCG.NREMstate.cellsmean.(CellClass.celltypes{morder(sm1)})(:,morder(sm2)),'linewidth',0.5,'color',modecolors{morder(sm1)})
            hold on
            box on
            
            if sm2==1
                ylabel(modenames{morder(sm1)})
            else
                set(gca,'ytick',[])
            end
            if sm1==6
                xlabel(modenames{morder(sm2)})
            else
                set(gca,'xtick',[])
            end
    end
    
    
end
NiceSave('AllModePopCCG',figfolder,baseName);



%%
GScouplingallcells = popCCG.NREMstate.cells.Mode6(:,CellClass.Mode6);
tempmerge = cat(1,ModeHMM.NREMstate(:).logrates);
GSrate = tempmerge(:,6);
tempmerge = cat(1,ModeHMM.NREMstate(:).cvs);
GScv = tempmerge(:,6);

[~,sortGSrate] = sort(GSrate);
[~,sortGSCV] = sort(GScv);

figure
subplot(2,2,1)
imagesc(GScouplingallcells(:,sortGSrate)')
colorbar
caxis([0 2])
xlabel('t - rel to GS spike')
ylabel('Sorted by GS Rate')
axis xy

subplot(2,2,2)
imagesc(GScouplingallcells(:,sortGSCV)')
colorbar
caxis([0 2])
xlabel('t - rel to GS spike')
ylabel('Sorted by GS CV')
axis xy

subplot(2,2,3)
plot(GSrate,GScv,'.')
LogScale('x',10)
xlabel('GS Rate (Hz)');ylabel('GS CV')


NiceSave('GSCouplingHeterogeneity',figfolder,baseName);
