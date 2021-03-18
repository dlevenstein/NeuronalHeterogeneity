function [] = bz_PlotModeRaster(spikemodes,modeintervals,plotcells,modecolors,win,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%   spiketimes      (numcells) structure with fields:
%       .prev_state
%       .next_state
%       .state_spk
%   modeintervals   (numcells) structure with field for time intervals of each mode
%   plotcells       which cells to plot (UID OR INDEX????!)
%   modecolors
%   win             time window to plot
%%
p = inputParser;
addParameter(p,'spikewidth',1)
addParameter(p,'linethick',2)
parse(p,varargin{:})
spikewidth = p.Results.spikewidth;
linethick = p.Results.linethick;


%%

plotnumcells = length(plotcells);

hold on
for cc = 1:plotnumcells
    whichcell = plotcells(cc);
    plotints =  structfun(@(modeints) RestrictInts(modeints,win,'inclusive',true),modeintervals(whichcell),'UniformOutput',false);
StateScorePlot( plotints,modecolors,'y',cc,'LineWidth',linethick)

for sm = 1:6
    %if ss==6
        instate_both = spikemodes(whichcell).prev_state == sm & spikemodes(whichcell).next_state==sm & ...
            InIntervals(spikemodes(whichcell).state_spk',win)';
    %else
        instate_either = (spikemodes(whichcell).prev_state == sm | spikemodes(whichcell).next_state==sm) & ...
            InIntervals(spikemodes(whichcell).state_spk',win)';
    %end
    plot([spikemodes(whichcell).state_spk(instate_either);spikemodes(whichcell).state_spk(instate_either)],...
        cc+[zeros(size(spikemodes(whichcell).state_spk(instate_either)))-0.4;0.4+zeros(size(spikemodes(whichcell).state_spk(instate_either)))],...
        'color',[0.5 0.5 0.5],'linewidth',spikewidth)
    plot([spikemodes(whichcell).state_spk(instate_both);spikemodes(whichcell).state_spk(instate_both)],...
        cc+[zeros(size(spikemodes(whichcell).state_spk(instate_both)))-0.4;0.4+zeros(size(spikemodes(whichcell).state_spk(instate_both)))],...
        'color','k','linewidth',spikewidth)
    
end
end
xlim(win)
ylim([0 plotnumcells+1])
box off
bz_ScaleBar('s')

end

