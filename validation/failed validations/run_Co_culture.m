%Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% May 2020. Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.
% Run a generic Co-culture simulation 


clear  
%close all 
%clc

% Difine flags 
flagCase = 'CoCulture';
kineticFlag = 'CoCulture';

%% set up 
reactor.V=1.2;            % reactor volume in L
reactor.pH(1)=1;          % Define whether pH is kept constant or not (yes)
reactor.pH(2)=6;          % pH value if it's controlled 
control = [1,1];          % REED controler and/or spikes 

%% prep structure parameters
parameters=get_parameters(reactor,flagCase,kineticFlag);
parameters.control = control;

% use most recent parameter values 
load('VeillionellaKinetics','par2names' ,'parameterReached')
load('parBootStrapBacillus','par2namesB' ,'parameterReportedB')

par2names = [par2names,par2namesB];
parameterReported = [parameterReached,parameterReportedB];

for i = 1:length(par2names)
    pos = strcmp(parameters.parAbb,par2names(i));
    parameters.parValues(pos) = parameterReported(i);
end
PRallPar(1,:) = par2names;
PRallPar(2,:) = num2cell(parameterReported);
disp('Os parámetros usado son:')
disp(PRallPar)
disp('')
% attitional tweeking
parameters.parValues(strcmp(parameters.parAbb,'lag_time')) = 0;


%% spike information BV5
spikeFlowrate = get_spike_data('Generic');
parameters.spikeFlowrate = spikeFlowrate;

%% initial conditions 
xInitial=get_initial_conditions(parameters,reactor,'Generic');

%% Reed data 
[Reed] = get_reed_data('Generic');
% give to structure parameters 
parameters.Reed = Reed; 
% tweek reed matrix to percentages 
ReedSwitch = 0;
parameters.ReedSwitch = ReedSwitch;
if ReedSwitch
parameters.Reed.RateMatrix(:,2) = 0.50;  %lactate
parameters.Reed.RateMatrix(:,4) = 0.25;  % proprionate 
parameters.Reed.RateMatrix(:,5) = 0.30;  % Acetate 
% parameters.Reed.RateMatrix(:,2) = percentLac;  %lactate
% parameters.Reed.RateMatrix(:,4) = percentPro;  % proprionate 
% parameters.Reed.RateMatrix(:,5) = percentAce;  % Acetate 
end


%% Manipulate variables
concFeed = 1000;
parameters.spikeFlowrate.compositionSpikes(2:end,1) = concFeed;


%% run ode15s
tsim = linspace(0,70,300);
tSpike =  spikeFlowrate.tSpikes;
tOde = union(tsim,tSpike); 
% optioins  
nneg = 1:length(xInitial);
options=odeset('RelTol',1e-4,'NonNegative',nneg,'MaxStep',0.1);
tic
[t,y] = ode15s(@f_mass_balances,tOde,xInitial,options,parameters);
toc


%% plot results 
title = 'Simulation of a generic Co-Culture experiment';
plot_CoCulture_simulation(t,y,parameters,title)










