function [SimValues] = Run_BrunelLIF(PopParams,TimeParams,varargin)
%LIF Model as used in Brunel 2000 and others
%by Jonathan Gornet and DLevenstein
%Last update: 4/8/2018

%INPUTS
%   PopParams       a structure that gives all the parameters of the population
%       .EPopNum	Number of excitatory neurons
%       .IPopNum	Number of inhibitory neurons
%       .u_0        Input current to the population. Can either be a
%                   constant, input to [E I] populations,
%                   or a function I_e(t) that returns input at time t
%                   time t. Add: 
%       .V_th       Membrane Threshold
%       .V_reset    Reset Potential
%
%       .tau_m      membrain time constant (ms)
%       .t_ref      Refractory period
%       .delay_s
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
%

%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%Parse optional inputs
p = inputParser;
addParameter(p,'showfig',true,@islogical)
addParameter(p,'showprogress',false,@islogical)
addParameter(p,'onsettime',0,@isnumeric)
addParameter(p,'save_dt',0.5,@isnumeric)
addParameter(p,'cellout',false,@islogical)
parse(p,varargin{:})
SHOWFIG = p.Results.showfig;
SHOWPROGRESS = p.Results.showprogress;
onsettime = p.Results.onsettime;
save_dt = p.Results.save_dt;
cellout = p.Results.cellout;

%--------------------------------------------------------------------------
%Simulation Parameters
EPopNum     = PopParams.EPopNum;    %Number of excitatory neurons
IPopNum     = PopParams.IPopNum;    %Number of excitatory neurons

PopNum      = EPopNum + IPopNum;    %Number of all neurons

SimTime     = TimeParams.SimTime;   %Simulation Time (ms)
dt          = TimeParams.dt;        %differential (ms)

%Calculate time vector from time parameters
SimTimeLength  = length([-onsettime:dt:SimTime]);   %Time Steps (simulated)
SaveTimeLength  = length([0:save_dt:SimTime]);      %Time Steps (saved)

%--------------------------------------------------------------------------
%Weight Matrices
EE_mat = zeros(PopNum);
II_mat = zeros(PopNum);
IE_mat = zeros(PopNum);
EI_mat = zeros(PopNum);

Ecells = 1:EPopNum;             EcellIDX = ismember(1:PopNum,Ecells);
Icells = EPopNum+1:PopNum;      IcellIDX = ismember(1:PopNum,Icells);

%Here we assign four 2x2 matrices of matrix (tensor?). There are positive values on the locations where there are connections.
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
%--------------------------------------------------------------------------
%Simulation Parameters
%LIF Parameters
tau_m       = PopParams.tau_m;       %membrane time constant (ms)
u_0         = PopParams.u_0;      %extenral drive (mV)
V_th        = PopParams.V_th;    %spike threshhold (mV)
V_reset     = PopParams.V_reset; %reset value (mV)

t_ref       = PopParams.t_ref;   %refractory period (ms)
delay_s     = PopParams.delay_s; %synaptic delay (ms)



%% Input: convert into function of t
% if isa(I_e, 'function_handle')
% elseif isequal(size(I_e),[1 1])
%     I_e = @(t) I_e;
% elseif length(I_e) == 2
%     I_e = @(t) transpose([I_e(1).*ones(1,EPopNum),     I_e(2).*ones(1,IPopNum)]);
% end

%% Variables

%Simulation Variables
V = zeros(PopNum,1);    %Membrane Potential
t_s = zeros(PopNum,1);  %synapse delay counter
t_r = zeros(PopNum,1); 	%refratory period counter
RI = 0;
u_ext = 0;

%Saved Variables
SimValues.t = nan(1,SaveTimeLength);
SimValues.V = nan(PopNum,SaveTimeLength);
SimValues.Input = nan(PopNum,SaveTimeLength);

spikes = nan(PopNum.*(SimTime+onsettime).*50,2); %assume mean rate 50Hz

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
    V0range = [0 V_th]; %make this neuron vector
    V(:,1) = V0range(1) + (1+p0spike).*diff(V0range).*rand(PopNum,1);
end
%% Time Loop
savecounter = 1;
timecounter = -onsettime-dt;
spikecounter = 0;
for tt=1:SimTimeLength
    %% Time Counter
    if SHOWPROGRESS && mod(tt,round(SimTimeLength./10))==0
        disp([num2str(round(100.*tt./SimTimeLength)),'% Done!']) %clearly, this needs improvement
    end
    %% Dynamics: update noise, V,s,w based on values in previous timestep
    
%     %V - Voltage Equation
%     dVdt =  (- g_L.*(V-E_L) ...                      %Leak
%              - g_w.*(V-E_w) ...                      %Adaptation
%              - g_e.*(V-E_e) - g_i.*(V-E_i) ...       %Synapses
%              + I_e(timecounter) + X_t)./C;           %External input
         
    dVdt =  (- V ...                                    %Leak
             + u_0 ...                                   %Constant Input
             + RI  ...                                   %Synapses
             + u_ext)./tau_m;                           %Time varying input
         
    %w - Adaptation Variable
    %dwdt = a_w.*(1-w) - b_w.*w; %Update to basic linear adaptation

    RI = 0; %reset the synaptic input to 0
    V   = V + dVdt.*dt;
    timecounter = round(timecounter+dt,4);  %Round to deal with computational error



    %% Spiking
    if any(V > V_th)
        %Find neurons that crossed threshold and record the spiketimes 
        spikeneurons = find(V > V_th);
        numspikers = length(spikeneurons);
        spikes(spikecounter+1:spikecounter+numspikers,:) = ...
            [timecounter.*ones(numspikers,1),spikeneurons];
        spikecounter = spikecounter+numspikers;
        
        %Start the synaptic delay counter
        t_s(spikeneurons) = delay_s;
        %Set spiking neurons refractory period 
        t_r(spikeneurons) = t_ref;
    end

    %%  Refractory period Countdowns
    if any(t_r > 0 | t_s > 0)
        refractoryneurons = t_r > 0;
        delayneurons = t_s > 0;
        
        %Hold voltage at rest
        V(refractoryneurons) = V_reset;
        %Count down the refractory periods
        t_r(refractoryneurons) = t_r(refractoryneurons) - dt;
        t_s(delayneurons) = t_s(delayneurons) - dt;
        
        %Pass along the spikes
        s = delayneurons & t_s<=0;
        RI = tau_m.*J_mat*s;
        
    end
        
    %% Add data to the output variables
    %Question: is accessing structure slower than doubles?
    if mod(timecounter,save_dt)==0 && timecounter>=0
         SimValues.t(savecounter)                 = timecounter;
         SimValues.V(:,savecounter)               = V;
         %SimValues.Input(:,savecounter)          = I_e(timecounter);
         
         savecounter = savecounter+1;
    end
    
    %%Idea: add a catch for silent network or excessive firing network?
end

%%
%Catch for no spiking in simulation error
spikes(spikecounter+1:end,:)=[];
if isempty(spikes); spikes = [nan nan]; end


%% Figure
if SHOWFIG
 
exneuron = randi(PopNum,1);
exspiketimes = spikes(spikes(:,2)==exneuron,1);
      
figure
    plot(spikes(:,1),spikes(:,2),'k.', 'Markersize' , 0.1)
    hold on
    plot([0 0],[0 PopNum],'r')
    xlabel('Time (ms)');ylabel('Neuron ID');title('Raster Plot');
    xlim([-onsettime SimTime]);ylim([0 PopNum+1]);
end
%% Output Structure

%Remove onset time
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
SimValues.WeightMat       = EE_mat+II_mat+EI_mat+IE_mat;


end