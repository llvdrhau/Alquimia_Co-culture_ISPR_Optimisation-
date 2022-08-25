%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run multi objective optimisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% October 2021.Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.

% Run a generic Co-culture simulation with multObjectiveFunc.m 
% this second script changes the upper bounds of the spaces 
% additionaly the simulation is run with 600 points in stead of 300 to gfet
% the extraction in the reservoir more precise 

clear  
close all 
clc

% which algorithm should be run 
patternSwitch = 0; 
paretoSwitch = 1;
geneticSwitch = 0; 

% estimate glucose spike concentrations gluON = 1?
gluConc = 400; %gCOD/L
gluON = 1;

% choose save name 
saveName = 'MultiObj';
%% Simulation set up
for k = 1
% Difine flags 
preRun = 0;
flagCase = 'CoCulture';
kineticFlag = 'CoCulture';

%% set up 
reactor.V=1.2;            % reactor volume in L
reactor.pH(1)=1;          % Define whether pH is kept constant or not (yes)
reactor.pH(2)=7;          % pH value if it's controlled 
control = [1,1];          % REED controler and/or spikes 


%% parameters from bootstrap should be used!!! 
% load parameters 
parameters=get_parameters(reactor,flagCase,kineticFlag);
parameters.control = control;
%retrive bootstrap parameters and replace 
load('parBootStrapVellionella','par2namesV' ,'parameterReportedV')
parameterReached = parameterReportedV;
par2names = par2namesV;
load('parBootStrapBacillus','par2namesB' ,'parameterReportedB')

par2names = [par2names,par2namesB];
parameterReported = [parameterReached,parameterReportedB];

for i = 1:length(par2names)
    pos = strcmp(parameters.parAbb,par2names(i));
    parameters.parValues(pos) = parameterReported(i);
end

% attitional tweeking
parameters.parValues(strcmp(parameters.parAbb,'lag_time')) = 0;
parameters.parValues(strcmp(parameters.parAbb,'kdec_Xlac')) = 0.043; % could use the original parameter for decay from the bootstrap if you reestimate the amount of veil added 
% print to screen
PRallPar(1,:) = par2names;
PRallPar(2,:) = num2cell(parameterReported);
disp('Os parámetros usado son:')
disp(PRallPar)
disp('')

%% spike information BV5
spikeFlowrate = get_spike_data('Generic');
parameters.spikeFlowrate = spikeFlowrate;

%% initial conditions 
xInitial=get_initial_conditions(parameters,reactor,'Generic');

%% Reed data 
[Reed] = get_reed_data('Generic');
% give to structure parameters 
parameters.Reed = Reed; % extraction is more or less the mean of each experiment also cte 

ReedSwitch = 0; % 1 if you want to use percentage extraction (we don't)
parameters.ReedSwitch = ReedSwitch;



%% run ode
if preRun 
    for i=1
tsim = linspace(0,65,300);
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

%% transfer to Reservoir
reservoir = zeros(length(t),14); 
volReservoir = 1; %L 
for i = 2:length(t)
    states = y(i,:)';
    [~,~,transferReed]=f_mass_balances(t(i),states,parameters); 
    reservoir(i,:) = reservoir(i-1,:) + transferReed'*(t(i)-t(i-1))*states(1)/volReservoir;
end
compoundAbb=parameters.compoundAbb; 
lacR = reservoir(:,strcmp(compoundAbb,'Slac'));
proR = reservoir(:,strcmp(compoundAbb,'Spro'));
aceR = reservoir(:,strcmp(compoundAbb,'Sac'));
figure
plot(t,lacR,t,proR,t,aceR)

legend('Lac','Pro','Ace')
ylabel('concentration g/L') 
xlabel('Time h')
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run optimisation
%% decisionVar
% decisionVar = [glu conc (optional) /// TspikeVeil /// start reed1 /// spacing] 
spikeVeil = 4;
start = 12;
nSpace = (length(Reed.tStart)-2) * 2 + 2 ;
spaces = ones(1,nSpace)*3;
 
parameters.gluSwitch = gluON; 
if gluON
    decisionVar = [gluConc, spikeVeil, start, spaces];
else
    decisionVar = [spikeVeil, start, spaces];
end 


%% boundries 

lbSpace = ones(1,nSpace)*2.05;
ubSpace =  ones(1,nSpace)*15;
if gluON
    %lb = [50, 0, 5, lbSpace];
    lb = [50, 0, 0.1, 2.05, 0.1, 2.05, 0.1, 2.05, 0.1 ]; % change when reed can start, not sure if effective? 
    ub = [1000, 20, 25, ubSpace];
else
    lb = [0, 0.1, lbSpace];
    ub = [20, 25, ubSpace];
end


%% optimisation algorithms 
%% Multi objective solution (yield and productivity) with patern search algorithm
% original multiobjective function

tic
rng default % For reproducibility
costFunc = @(decisionVar)multiObjectiveFunc(decisionVar, parameters, xInitial);
nVar = length(lb);
%options = optimoptions('paretosearch','ParetoSetSize',200);
options = optimoptions('paretosearch','ParetoSetSize',250 ,'PlotFcn','psplotparetof'); % large ParetoSize
[decisionVarOut,fxVal]= paretosearch(costFunc,nVar,[],[],[],[],lb,ub,[],options);
toc
% save results
save(saveName,'decisionVarOut','fxVal','parameters','xInitial')





