function [RSME,r2] = f_run_co_culture(kinetics,sheetData,parameterValues,parNames,giveTitleToGraf)
%Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% May 2020. Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.
% Run a simulation of BV5 or BV4, helps to determin ideal stoichiometry 

% clear  
% close all 
% clc

flagCase = 'CoCulture';
kineticFlag = 'CoCulture';

%% Reading experimental data
[data1,Names]=xlsread('expData.xlsx',sheetData);

tExp=data1(:,strcmp(Names,'time'));

yExp1(:,1)=data1(:,strcmp(Names,'Glu'));
yExp1(:,2)=data1(:,strcmp(Names,'Lac'));
yExp1(:,3)=data1(:,strcmp(Names,'Pro'));
yExp1(:,4)=data1(:,strcmp(Names,'Ace'));
yExp1(:,5)=data1(:,strcmp(Names,'BM'));

%% set up 
reactor.V=1.3; %L
reactor.pH(1)=1;             % Define whether pH is kept constant or not
reactor.pH(2)=7;             % pH value if it's controlled

% REED controler and/or spikes 
control = [1,1];

%prep structure parameters
parameters=get_parameters(reactor,flagCase,kineticFlag);
parameters.control = control;
parameters.kineticInhibition = kinetics;

%% spike information
spikeFlowrate = get_spike_data(sheetData);
parameters.spikeFlowrate = spikeFlowrate;

%% initial conditions 
xInitial=get_initial_conditions(parameters,reactor,yExp1);

%% Reed data 
[Reed] = get_reed_data(sheetData);
% give to structure parameters 
parameters.Reed = Reed; 

%% use parameter values from veilionella V3/4 and bootstrap B1/2
% parameters Veilionella
parameterReportedV =parameterValues; 
par2namesV = parNames;  
% parameters bacillus
load('parBootStrapBacillus','par2namesB' ,'parameterReportedB')
% Mech the two parameters togegther 
par2names = [par2namesV,par2namesB];
parameterReported = [parameterReportedV,parameterReportedB];

for i = 1:length(par2names)
    pos = strcmp(parameters.parAbb,par2names(i));
    parameters.parValues(pos) = parameterReported(i);
end

% attitional tweeking: make sure the lag time is zero 
parameters.parValues(strcmp(parameters.parAbb,'lag_time')) = 0;



%% edit the extraction process with mean removal rate  
meanExtractionRate = 1;
if meanExtractionRate
    ReedRateMatrix = Reed.RateMatrix;
    for i = 1:size(ReedRateMatrix,2)    % TODO check if this changes each row to it's corresponding mean rate
        parameters.Reed.RateMatrix(:,i) = mean(ReedRateMatrix(:,i));  %lactate
    end
end

%% run ode15s
tsim = linspace(0,tExp(end),300);
tSpike =  spikeFlowrate.tSpikes;
tOde = union(tsim,tSpike); 
% optioins  
nneg = 1:length(xInitial);
options=odeset('RelTol',1e-4,'NonNegative',nneg,'MaxStep',0.1);
tic
[t,y] = ode15s(@f_mass_balances,tOde,xInitial,options,parameters);
toc
%% plot results 
title = giveTitleToGraf;
plot_CoCulture(tExp,yExp1,t,y,parameters,title)

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


%% RMSE
compounds = {'Glu', 'Lac', 'Pro', 'Ace', 'BM'};

pos_glu = find((strcmp(parameters.compoundAbb,'Sglu')));
pos_lac = find((strcmp(parameters.compoundAbb,'Slac')));
pos_pro = find((strcmp(parameters.compoundAbb,'Spro')));
pos_ac = find((strcmp(parameters.compoundAbb,'Sac')));
pos_bm = find((strcmp(parameters.compoundAbb,'Xlac')));
pos_bm2 = find((strcmp(parameters.compoundAbb,'Xbac')));
%pos_ye = (strcmp(parameters.compoundAbb,'Sye'));
col = [pos_glu, pos_lac,pos_pro,pos_ac,pos_bm,pos_bm2];
ysim = y(:,col); 
[RSME, r2] = f_RMSE(tExp,yExp1,ysim,tsim,compounds);
end

