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
% 
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
PopParams.EPopNum = 8000;
PopParams.IPopNum = 2000;
PopParams.u_0 = 24;
PopParams.V_th =20;
PopParams.tau_m = 20;
PopParams.V_reset = 10;
PopParams.t_ref = 0.5;
PopParams.delay_s = 0.55;
PopParams.J = 0.2;
PopParams.J = 0.8;
PopParams.g = 5;
PopParams.Kee = 800;
PopParams.Kie = 800;
PopParams.Kei = 200;
PopParams.Kii = 200;

TimeParams.dt = 0.05;
TimeParams.SimTime = 100;
[SimValues] = Run_BrunelLIF(PopParams,TimeParams,'showprogress',true);