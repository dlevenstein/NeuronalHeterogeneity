function [SimValues] = Run_LIF_iSTDP(PopParams,TimeParams,varargin)
%LIF Model as used in Brunel 2000 and others
%Jonathan Gornet and DLevenstein 2017-2019

%INPUTS
%   PopParams       a structure that gives all the parameters of the population
%       .EPopNum	Number of excitatory neurons
%       .IPopNum	Number of inhibitory neurons
%       .u_0        Input current to the population. Can either be a
%                   constant, input to [E I] populations,
%                   or a function I_e(t) that returns input at time t
%                   time t. Add: 
%       .u_ext(t)   A function for time-varying external input
%       .V_th       Membrane Threshold
%       .V_reset    Reset Potential
%       .V_rest     Resting Potential
%
%       .tau_m      membrane time constant (ms)
%       .t_ref      Refractory period
%       .delay_s    Delay from cell to ALL postsynaptic partners.
%                   To Do: make synapse specific...
%
%                   SYNAPTIC WEIGHT (voltage jump: mV)
%       .J          Extiatory synaptic weight
%       .g          Inhibitory synaptic weight = -g*J
%                   CONNECTION PROBABILITY (Expected In-Degree)
%       .Kee        E->E
%       .Kii        I->I
%       .Kie        E->I
%       .Kei        I->E
%
%       .p0spike    (optional) proportion of neurons spiking at t0 (default:0)
%   TimeParams
%       .dt        timestep (ms)
%       .SimTime   total simulation time (ms)
%   'showfig'       (optional) show the figure? (default:true)
%   'showprogress'  (optional) show time counter of progress (default:false)
%   'onsettime'     (optional) duration of (removed?) onset time (default: 0ms)
%   'save_dt'       (optional) dt for the saved output (default: 0.1ms)
%   'cellout'
%   'estrate'       (default: 10Hz) lowering this saves RAM,
%   'J_mat'         put in a weight matrix to start with
%

%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%Parse optional inputs
p = inputParser;
addParameter(p,'showfig',true,@islogical)
addParameter(p,'showprogress',false)
addParameter(p,'onsettime',0,@isnumeric)
addParameter(p,'save_dt',0.5,@isnumeric)
addParameter(p,'cellout',false,@islogical)
addParameter(p,'estrate',20)
addParameter(p,'J_mat',[])
addParameter(p,'plotEIweight',false)
addParameter(p,'plotoverlay',[])
parse(p,varargin{:})
SHOWFIG = p.Results.showfig;
SHOWPROGRESS = p.Results.showprogress;
onsettime = p.Results.onsettime;
save_dt = p.Results.save_dt;
cellout = p.Results.cellout;
estrate = p.Results.estrate; 
J_mat = p.Results.J_mat; 
PLOTEI = p.Results.plotEIweight; 
plotoverlay = p.Results.plotoverlay; 

%% Default Parameters
DefaultParms.EPopNum = 8000;
DefaultParms.EPopNum = 2000;
DefaultParms.IPopNum = 2000;
DefaultParms.IPopNum = 500;
DefaultParms.u_0 = 24;
DefaultParms.u_ext = @(t) 0;
DefaultParms.V_th =20;
DefaultParms.tau_m = 20;
DefaultParms.V_reset = 10;
DefaultParms.V_rest = 0;
DefaultParms.t_ref = 0.5;
DefaultParms.delay_s = 0.55;
DefaultParms.J = 0.2;
DefaultParms.J = 0.8;
DefaultParms.g = 5;
DefaultParms.Kee = 800;
DefaultParms.Kie = 800;
DefaultParms.Kei = 200;
DefaultParms.Kii = 200;

DefaultParms.LearningRate = 0;

PopParams = EnterDefaultParms(PopParams,DefaultParms);


%%
%--------------------------------------------------------------------------
%Simulation Parameters
EPopNum     = PopParams.EPopNum;    %Number of excitatory neurons
IPopNum     = PopParams.IPopNum;    %Number of excitatory neurons

PopNum      = EPopNum + IPopNum;    %Number of all neurons

SimTime     = TimeParams.SimTime;   %Simulation Time (ms)
dt          = TimeParams.dt;        %differential (ms)

%LIF Parameters
tau_m       = PopParams.tau_m;       %membrane time constant (ms)
u_0         = PopParams.u_0;      %extenral drive (mV)
V_th        = PopParams.V_th.*ones(PopNum,1);    %spike threshhold (mV)
V_reset     = PopParams.V_reset; %reset value (mV)
V_rest     = PopParams.V_rest.*ones(PopNum,1); %reset value (mV)
u_ext = PopParams.u_ext; %function of time gets replaced later if external poisson

t_ref       = PopParams.t_ref;   %refractory period (ms)
delay_s     = PopParams.delay_s.*ones(PopNum,1); %synaptic delay (ms)

%iSTDP Parameters
LearningRate = PopParams.LearningRate;
TargetRate   = PopParams.TargetRate; %Target Rate for Excitatory cells (units of Hz)
tauSTDP      = PopParams.tauSTDP;    %Time Constant for the STDP curve (Units of ms)

alpha = 2.*(TargetRate./1000).*tauSTDP; %Alpha parameter from Vogels eqn5
PlasticPostCells = ~isnan(TargetRate);  %Inhibitory synapses on which postsynaptic cells are plastic?

%--------------------------------------------------------------------------
%Weight Matrices
Ecells = 1:EPopNum;             EcellIDX = ismember(1:PopNum,Ecells)';
Icells = EPopNum+1:PopNum;      IcellIDX = ismember(1:PopNum,Icells)';
EE_mat = zeros(PopNum);
II_mat = zeros(PopNum);
IE_mat = zeros(PopNum);
EI_mat = zeros(PopNum);

if ~isempty(J_mat)
    disp('Using weight matrix provided...')
    EE_mat(Ecells,Ecells) = J_mat(Ecells,Ecells);
    II_mat(Icells,Icells) = J_mat(Icells,Icells);
    IE_mat(Icells,Ecells) = J_mat(Icells,Ecells);
    EI_mat(Ecells,Icells) = J_mat(Ecells,Icells);
else

    %Here we assign four 2x2 matrices of matrix. There are positive values on the locations where there are connections.
    %For example, there are values for the EE connections on the 1x1 matrix, II
    %on the 2x2 matrix, and etc (this is based on the indexing of the neuron population). 

    %NOTE: presynaptic neurons are columns (dim2) and postsynaptic neurons are rows (dim1).

    %E->E Synapses
    Wee = PopParams.J;
    Kee = PopParams.Kee;
    Pee = Kee./(EPopNum-1); %-1 to account for self-connections (which are then removed)

    EE_mat(Ecells,Ecells) = rand(EPopNum)<=Pee;
    EE_mat = EE_mat.*Wee;
    EE_mat(diag(diag(true(size(EE_mat)))))=0; %Remove selfconnections

    %I->I Synapses
    Wii = -PopParams.J.*PopParams.g;
    Kii = PopParams.Kii;
    Pii = Kii./(IPopNum-1);

    II_mat(Icells,Icells) = rand(IPopNum)<=Pii;
    II_mat = II_mat.*Wii;
    II_mat(diag(diag(true(size(II_mat)))))=0; %Remove selfconnections

    %E->I Synapses
    Wie = PopParams.J;
    Kie = PopParams.Kie;
    Pie = Kie./EPopNum;

    IE_mat(Icells,Ecells) = rand(IPopNum,EPopNum)<=Pie;
    IE_mat = IE_mat.*Wie;

    %I->E Synapses
    Wei = -PopParams.J.*PopParams.g;
    Kei = PopParams.Kei;
    Pei = Kei./IPopNum;

    EI_mat(Ecells,Icells) = rand(EPopNum,IPopNum)<=Pei;
    EI_mat = EI_mat.*Wei;

    J_mat = EE_mat + EI_mat + II_mat + IE_mat;
end
isconnected = J_mat~=0;
%% Poisson spiking input

if isfield(PopParams,'ex_rate') && isfield(PopParams,'N_FF')  && ...
        isfield(PopParams,'K_FF') && isfield(PopParams,'J_FF')
    if SHOWPROGRESS
    	disp('Making FF inputs...')
    end
    
    %check that FFweights is the right size matrix (N_FF x PopNum?)
    Pff = PopParams.K_FF./PopParams.N_FF;
    if length(Pff)==2
        Pff = [Pff(1).*ones(PopParams.N_FF,EPopNum),...
            Pff(2).*ones(PopParams.N_FF,IPopNum)];
    end
    FF_mat = rand(PopParams.N_FF,PopNum)<=Pff;
    
	if length(PopParams.J_FF)==2
        Jff = [PopParams.J_FF(1).*ones(PopParams.N_FF,EPopNum),...
            PopParams.J_FF(2).*ones(PopParams.N_FF,IPopNum)];
    else%if length(PopParams.J_FF)==2
        Jff = PopParams.J_FF;
    end
    
    FF_mat = FF_mat.*Jff;
    SimValues.FF_mat = FF_mat;
    %Note: here is where we can have input follow some rate 
    %(PopParams.ex_rate can be N_t x N_FF or a function of t!) note this
    %might slow stuff down... check.
    if ~isa(PopParams.ex_rate, 'function_handle')
        PopParams.ex_rate = @(t) PopParams.ex_rate;
    end
    u_ext = @(t) ((rand(1,PopParams.N_FF)<(PopParams.ex_rate(t).*dt./1000))*(FF_mat))'.*tau_m./dt; 
elseif isfield(PopParams,'ex_rate')
    u_ext = @(t) (mean(PopParams.J(:)).*(rand(PopNum,1)<(PopParams.ex_rate.*dt./1000))).*tau_m./dt; 
end

%% Initialize Variables
disp('Initializing variables...')
%Simulation Variables
V = zeros(PopNum,1);    %Membrane Potential
% xpre = zeros(size(delay_s)); %PreSynaptic trace
% xpost = zeros(size(delay_s)); %PostSynaptic trace
x = zeros(size(delay_s)); %PostSynaptic trace
t_s = zeros(size(delay_s))-dt;  %synapse delay counter
t_r = zeros(PopNum,1); 	%refratory period counter
RI = 0;


%Calculate time vector from time parameters
SimTimeLength  = length([-onsettime:dt:SimTime]);   %Time Steps (simulated)
SaveTimeLength  = length([0:save_dt:SimTime]);      %Time Steps (saved)

%Saved Variables
SimValues.t = nan(1,SaveTimeLength);
SimValues.V = nan(PopNum,SaveTimeLength);
%SimValues.Input = nan(PopNum,SaveTimeLength);
%SimValues.J_mat = nan(PopNum,PopNum,SaveTimeLength);

estnumspikes = PopNum.*(SimTime+onsettime).*estrate./1000;
spikes = nan(estnumspikes,2,'single');
haswarned = false; %to warn the user if their estimate is low.
%% Initial Conditions - random voltages
%Improvement: set # initial spiking neurons instead of hard coding 
%range: E_L-Vth
if isfield(PopParams,'p0spike') 
    p0spike = PopParams.p0spike;
else
    p0spike = 0.0; %5 chance of initial spiking 
end
if isfield(PopParams,'V0')
    V(:,1) = PopParams.V0;
else
    V0range = [V_rest V_th]; %make this neuron vector
    V(:,1) = V0range(:,1) + (1+p0spike).*diff(V0range,[],2).*rand(PopNum,1);
end
%% Time Loop
savecounter = 1;
timecounter = -onsettime-dt;
spikecounter = 0;
if SHOWPROGRESS
    disp('Starting Simulation...')
end
for tt=1:SimTimeLength
    %% Time Counter
    timecounter = round(timecounter+dt,4);  %Round to deal with computational error
    if strcmp(SHOWPROGRESS,'parloop')   
        if mod(tt,round(SimTimeLength./10))==1
            disp([num2str(round(100.*tt./SimTimeLength)),'% Done!'])
        end
    elseif SHOWPROGRESS && mod(tt,round(SimTimeLength./100))==0
        bz_Counter(round(100.*tt./SimTimeLength),100,'Percent Complete')
    end

    %% Dynamics: update noise, V,s,w based on values in previous timestep
    
	%V - Voltage Equation
    dVdt =  (- V + V_rest ...                                    %Leak
             + u_0 ...                                   %Constant Input
             + RI  ...                                   %Synapses
             + u_ext(timecounter))./tau_m;                           %Time varying input
         
	%x - Synaptic Trace for STDP. Note, if only EI plasticity, could make
	%only one trace, where e is presynaptic and i is postsynaptic
    %dxpredt =  - xpre./tauSTDP;
    %dxpostdt =  - xpost./tauSTDP;
    dxdt =  - x./tauSTDP;
    
    V  = V + dVdt.*dt;
%     xpre   = xpre + dxpredt.*dt;
%     xpost   = xpost + dxpostdt.*dt;
    x   = x + dxdt.*dt;

    %% Spiking
    RI = 0;  %reset the synaptic input to 0.
    if any(V > V_th)
        %Find neurons that crossed threshold and record the spiketimes 
        spikeneurons = V > V_th;
        %Register the cell ID/timestamp of the spiked neurons in the spikes vector
        numspikers = sum(spikeneurons);
        spikes(spikecounter+1:spikecounter+numspikers,:) = ...
            [timecounter.*ones(numspikers,1),find(spikeneurons)];
        spikecounter = spikecounter+numspikers;
        if spikecounter > estnumspikes &~ haswarned
            display('WARNING - NUMBER SPIKES IS ABOVE ESTIMATED. GOING MUCH SLOWER NOW')
            haswarned = true;
        end
        
        %Start the synaptic delay counter %BIG PROBLEM HERE! IF cell spikes
        %during delay, delay resets and spike doesn't get through! (e.g. inhibitory cells...). Keep
        %delays short to fix this as a bandaid for now.......
        %Another option is to not jump delay if neuron is mid delay... (implemented)
        %Switch the order of spiking and refractory period countdowns?
        t_s(spikeneurons & t_s<0) = delay_s(spikeneurons & t_s<0);
        %Set spiking neurons refractory period 
        t_r(spikeneurons) = t_ref;
        
        %Jump the postsynaptic trace of the spiking neurons (E only for iSTDP)
        x(spikeneurons & EcellIDX) = x(spikeneurons & EcellIDX) + 1;
        %xpost(spikeneurons) = xpost(spikeneurons) + 1;
        
    else
        numspikers=0;
    end

    %%  Refractory period Countdowns
    if any(t_r > 0) || any(t_s(:) >= 0)
        refractoryneurons = t_r > 0;
        delaysynapses = t_s >= 0;
        
        %Hold voltage at rest
        V(refractoryneurons) = V_reset;
        %Count down the refractory periods
        t_r(refractoryneurons) = t_r(refractoryneurons) - dt;
        t_s(delaysynapses) = t_s(delaysynapses) - dt;
        
        %Pass along the spikes - cells that were in delay, but are now below
        s = delaysynapses & t_s<0;   %activated synapses
        if any(s(:))
            %Jump the presynaptic trace (I only for iSTDP)
            x(s & IcellIDX) = x(s & IcellIDX) + 1;
            %xpre(s) = xpre(s) + 1;

            %Apply the current in the next timestep
            RI = tau_m.*J_mat*s./dt;
        end
    else 
        s = 0;
    end
        
        %HERE: implement matrix delays...
    if (any(s(:)) || numspikers>0) && (LearningRate ~= 0) 
        %Implement STDP (Vogels 2011 SuppEqn 4/5) I->E only
        %Presynaptic I Cells -  adjust all synapses postsynaptic to spiking I cells
        %                       strengthen I->E if recently active, weaken if not
        PreIspikes = s & IcellIDX;
        %PlasticPostCells = EcellIDX;
        %PlasticPostCells = true(size(EcellIDX));
        J_mat(PlasticPostCells,PreIspikes) = J_mat(PlasticPostCells,PreIspikes) - LearningRate.*(x(PlasticPostCells)-alpha(PlasticPostCells)); 

        %Postsynaptic E cells - adjust all (recently active inhibitory) 
        %                       synapses presynaptic to spiking E cells
        PostSpikers = spikeneurons & PlasticPostCells;            
        J_mat(PostSpikers,IcellIDX) = J_mat(PostSpikers,IcellIDX) - LearningRate.*(x(IcellIDX)');
        %(Negative because inhibitory)
        
        %if sum(PreIspikes)>1
         %   keyboard
        %end

        J_mat(~isconnected) = 0; %Get rid of any negative synapses and unconnected pairs...
        J_mat(J_mat(:,IcellIDX)>0) = 0;
    end
    %% Add data to the output variables
    %Question: is accessing structure slower than doubles?
    if mod(timecounter,save_dt)==0 && timecounter>=0
         SimValues.t(savecounter)                 = timecounter;
         SimValues.V(:,savecounter)               = V;
         
         EIconnections = J_mat(Ecells,Icells);
         
         SimValues.EImean(savecounter)   = mean(EIconnections(isconnected(Ecells,Icells))); 
         SimValues.EIstd(savecounter)   = std(EIconnections(isconnected(Ecells,Icells))); 
         
         IIconnections = J_mat(Icells,Icells);
         
         SimValues.IImean(savecounter)   = mean(IIconnections(isconnected(Icells,Icells))); 
         SimValues.IIstd(savecounter)   = std(IIconnections(isconnected(Icells,Icells))); 
         %SimValues.Input(:,savecounter)          = I_e(timecounter);
         %SimValues.J_mat(:,:,savecounter)         = J_mat;
         savecounter = savecounter+1;
    end
    
    %%Idea: add a catch for silent network or excessive firing network?
end
if strcmp(SHOWPROGRESS,'parloop')
    disp('Simulation finished...')
end
%%
%Catch for no spiking in simulation error
spikes(spikecounter+1:end,:)=[];
if isempty(spikes); spikes = [nan nan]; end


%% Output Structure

%Remove onset time

disp(['Used ',num2str(spikecounter./estnumspikes.*100),'% of prediced spikes'])

spikes(spikes(:,1)<=0,:) = [];

if cellout
    for cc = 1:PopNum %This can go very slow with lots of spikes....
        spikesbycell{cc} = spikes(spikes(:,2)==cc,1);
    end
    SimValues.spikesbycell    = spikesbycell;
end


SimValues.spikes          = spikes;
SimValues.EcellIDX        = Ecells;
SimValues.IcellIDX        = Icells;
SimValues.WeightMat_initial       = EE_mat+II_mat+EI_mat+IE_mat;
SimValues.WeightMat       = J_mat;
SimValues.PopParams = PopParams;
SimValues.TimeParams = TimeParams;
SimValues.optionalinputs = p.Results;

if SHOWPROGRESS
    %disp('Saved to Structure')
end

%%
if SHOWFIG
    try
        %disp('Plotting')
        PlotSimRaster(SimValues,[-onsettime SimTime],'plotEIweight',PLOTEI,...
            'overlay',plotoverlay);
        %disp('Plot Success!')
    catch
        disp('Failed to plot...')
    end
end
end