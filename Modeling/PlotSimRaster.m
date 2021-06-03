function [spikemat] = PlotSimRaster(SimValues,timewin,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%
p = inputParser;
addParameter(p,'cellsort',[])
addParameter(p,'overlay',[])
addParameter(p,'plotEIweight',false)
addParameter(p,'trainingfigure',false)
addParameter(p,'title',[])
addParameter(p,'ratebin',5) %ms
parse(p,varargin{:})
cellsort = p.Results.cellsort;
overlay = p.Results.overlay;
ratebin = p.Results.ratebin;
PLOTEI = p.Results.plotEIweight;
TRAINFIG = p.Results.trainingfigure;
plottitle = p.Results.title;

%%
if ~exist('timewin','var') || isempty(timewin)
    timewin = [0 SimValues.TimeParams.SimTime];
end

%%
    spikes = SimValues.spikes;
    Ecells = SimValues.EcellIDX;
    Icells = SimValues.IcellIDX;
    PopNum = length(Ecells) + length(Icells);

    spikes = spikes(InIntervals(spikes(:,1),timewin),:);
    
    celltypes = {'E','I'};
    cellcolors = {'k','r'};
    cellspikes.E = ismember(spikes(:,2),Ecells);
    cellspikes.I = ismember(spikes(:,2),Icells);

    for cc = 1:length(celltypes)
        classspikes = spikes(cellspikes.(celltypes{cc}),:);
        classspikes(:,2) = classspikes(:,2)-min(classspikes(:,2))+1;
        spikemat.(celltypes{cc}) = ...
            bz_SpktToSpkmat(classspikes,'dt',1,'binsize',ratebin,'units','rate','win',timewin);
        
        spikemat.(celltypes{cc}).poprate = mean(spikemat.(celltypes{cc}).data,2).*1000;
    end
    
exneuron = randsample(Ecells,1);
exspiketimes = spikes(spikes(:,2)==exneuron,1);

%Relabel cells by sort
if ~isempty(cellsort)
    [~,sortraster] = sort(cellsort);
    spikes(:,2) = sortraster(spikes(:,2));
end


   %%   
figure
subplot(2,1,1)
    hold on
    for cc = 1:length(celltypes)
        plot(spikes(cellspikes.(celltypes{cc}),1),spikes(cellspikes.(celltypes{cc}),2),...
            '.', 'Markersize' , 0.1,'color',cellcolors{cc})
    end
    if ~isempty(overlay)
        plot(overlay(:,1),bz_NormToRange(overlay(:,2),0.5))
    end
    box off
    plot([0 0],[0 PopNum],'r')
    ylabel('Neuron ID');
    title(plottitle);
    xlim(timewin);ylim([0 PopNum+1]);
    %ylim([0 100])
    bz_ScaleBar('ms')
subplot(4,1,3)
    hold on
    for cc = 1:length(celltypes)
        plot(spikemat.(celltypes{cc}).timestamps,spikemat.(celltypes{cc}).poprate,cellcolors{cc})
    end
    ylabel('Pop Rate (Hz)')
    raterange = ylim;
    ylim([0 raterange(2)]);xlim(timewin)
    box off
    
    
    if PLOTEI
        try
        subplot(4,1,4)
            hold on
            plot(SimValues.t,SimValues.EImean,'k')
            errorshade(SimValues.t,SimValues.EImean,...
                SimValues.EImean+SimValues.EIstd,SimValues.EImean-SimValues.EIstd,...
                'k','vector')
            
            plot(SimValues.t,SimValues.IImean,'r')
            errorshade(SimValues.t,SimValues.IImean,...
                SimValues.IImean+SimValues.IIstd,SimValues.IImean-SimValues.IIstd,...
                'r','vector')
            
            ylabel('Inh Weights')
        catch
            display('something wrong with EIweight data. Womp Womp')
        end
    else
       try
            V_th = SimValues.PopParams.V_th(1); %excitatory... assuming first
            V_rest = SimValues.PopParams.V_rest(1); %excitatory... assuming first

            subplot(4,1,4)
                plot(SimValues.t,SimValues.V(exneuron,:),'k')
                hold on
                plot(exspiketimes,V_th.*ones(size(exspiketimes))+2,'k.')
                box off
                plot(xlim,V_th.*[1 1],'k--')
               % plot(spikes

                xlabel('Time (ms)');ylabel('V, example cell')
                xlim(timewin);ylim([V_rest V_th+2])
       catch
           display('something wrong with voltage data. Womp Womp')
       end
    end
    
    
    %% Training figure
    if TRAINFIG
        %KEI over time
        %Rate minus TR over time
        trainmat = ...
            bz_SpktToSpkmat(spikes,'dt',500,'binsize',1000,'units','rate','win',timewin);
        %%
        %Rate vs TR
        %KEI vs TR at end
        KE = sum(SimValues.WeightMat(SimValues.EcellIDX,SimValues.EcellIDX),2);
        KI = sum(SimValues.WeightMat(SimValues.EcellIDX,SimValues.IcellIDX),2);
        KEI = KE./KI;
        %Variance over time, mean I over time
        
        figure
        subplot(2,2,1)
        plot(log10(SimValues.PopParams.TargetRate(SimValues.EcellIDX)),KEI,'.')
        
        subplot(2,1,2)
        plot(bsxfun(@minus,trainmat.data(:,SimValues.EcellIDX)*1000,SimValues.PopParams.TargetRate(SimValues.EcellIDX)')')

    end
end

