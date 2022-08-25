%Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% May 2020. Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.
% Run a simulation of BV3 BV4, BV5

clear  
close all 
clc

% BV5
flagCase = 'CoCulture';
kineticFlag = 'CoCulture';

%% Reading experimental data
[data1,Names]=xlsread('expData.xlsx','BV4');

tExp=data1(:,strcmp(Names,'time'));

yExp1(:,1)=data1(:,strcmp(Names,'Glu'));
yExp1(:,2)=data1(:,strcmp(Names,'Lac'));
yExp1(:,3)=data1(:,strcmp(Names,'Pro'));
yExp1(:,4)=data1(:,strcmp(Names,'Ace'));
yExp1(:,5)=data1(:,strcmp(Names,'BM'));

%% set up 
reactor.V=1.295; %L
reactor.pH(1)=1;             % Define whether pH is kept constant or not
reactor.pH(2)=7;             % pH value if it's controlled

% REED controlers and/or spikes 
control = [1,1];

%prep structure parameters
parameters=get_parameters(reactor,flagCase,kineticFlag);
parameters.control = control;

%% use parameter values from bootstrap or first estimation 

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
%% spike information BV3
spikeFlowrate = get_spike_data('BV4');
parameters.spikeFlowrate = spikeFlowrate;

%% initial conditions 
xInitial=get_initial_conditions(parameters,reactor,yExp1);

%% Reed data 
[Reed] = get_reed_data('BV4');
% give to structure parameters 
parameters.Reed = Reed; 

%% Manipulate parameters 
% parameters.parValues(strcmp(parameters.parAbb,'km_su')) = 2.9298;
 parameters.parValues(strcmp(parameters.parAbb,'lag_time')) = 0;
% parameters.parValues(strcmp(parameters.parAbb,'Ysu')) = 0.0661;

 parameters.spikeFlowrate.compositionSpikes(1,11) = 1*1.34 ;     %Spike van Veillionella
 parameters.spikeFlowrate.compositionSpikes(2,1) = 500*1.067 ;     %1st spike 
% parameters.spikeFlowrate.compositionSpikes(4,1) = 1067*1 ;     %3d spike 
% parameters.spikeFlowrate.compositionSpikes(3,1) = 1067*1 ;     %4th spike


%% edit the extraction process with percent extraction 
ReedSwitch = 0;
parameters.ReedSwitch = ReedSwitch;

if ReedSwitch
parameters.Reed.RateMatrix(:,2) = 0.13;  %lactate
parameters.Reed.RateMatrix(:,4) = 0.25;  % proprionate 
parameters.Reed.RateMatrix(:,5) = 0.30;  % Acetate 
% parameters.Reed.RateMatrix(:,2) = percentLac;  %lactate
% parameters.Reed.RateMatrix(:,4) = percentPro;  % proprionate 
% parameters.Reed.RateMatrix(:,5) = percentAce;  % Acetate 
end


%% run ode15s
tsim = linspace(0,tExp(end),300);
tSpike =  spikeFlowrate.tSpikes;
tOde = union(tsim,tSpike); 
% optioins  
nneg = 1:length(xInitial);
options=odeset('RelTol',1e-4,'NonNegative',nneg,'MaxStep',0.1);
tic
[t,y] = ode45(@f_mass_balances,tOde,xInitial,options,parameters);
toc

%% R2  coefficient of determination
compoundNames = {'Sglu','Slac' 'Spro' 'Sac' 'BM' };
compoundAbb=parameters.compoundAbb;
glu = y(:,strcmp(compoundAbb,'Sglu'));
lac = y(:,strcmp(compoundAbb,'Slac'));
pro = y(:,strcmp(compoundAbb,'Spro'));
ace = y(:,strcmp(compoundAbb,'Sac'));
Xlac = y(:,strcmp(compoundAbb,'Xlac'));
Xbac = y(:,strcmp(compoundAbb,'Xsu'));

yR2 = [glu,lac,pro,ace,Xlac,Xbac];

[RMSENor, r2] = f_RMSE(tExp,yExp1,yR2,t,compoundNames);

%% plot results 
plot_CoCulture(tExp,yExp1,t,y,parameters,'BV4 Co-Culture',r2)


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
hold on 
tRV=data1(:,strcmp(Names,'timeRV'));
gluRV = data1(:,strcmp(Names,'GluRV'));
lacRV = data1(:,strcmp(Names,'LacRV'));
proRV = data1(:,strcmp(Names,'ProRV'));
acRV = data1(:,strcmp(Names,'AceRV'));
plot(tRV,lacRV,'bx',tRV,proRV,'rx',tRV,acRV,'yx')

legend('Lac','Pro','Ace')
ylabel('concentration g/L') 
xlabel('Time h')










