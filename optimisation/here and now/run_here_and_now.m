%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN HERE AND NOW script 
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% October 2021.Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.

% optimisation under uncertanty according to the methode of here_and_now
% uses objective function f_objcFunc_stochastic (max yield or productivity)
% optimisation algorithm: pattern search 


clear  
close all 
clc

%% Run a generic Co-culture simulation 

% estimate glucose spike concentrations gluON = 1?
gluConc = 200; %gCOD/L
gluON = 1;

% switch to see initial simulation set up 
preRun = 0;

% choose save name 
saveName = 'here_and_now_constrained2';
objective = 'productivity';
%% Simulation set up
for k = 1
% Difine flags 
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
parameters.parValues(strcmp(parameters.parAbb,'kdec_Xlac')) = 0.043; 
% could use the original parameter for decay from the bootstrap if you reestimate the amount of veil added 
% % print to screen
% PRallPar(1,:) = par2names;
% PRallPar(2,:) = num2cell(parameterReported);
% disp('Os parámetros usado son:')
% disp(PRallPar)
% disp('')

%% spike information BV5
spikeFlowrate = get_spike_data('Generic');
parameters.spikeFlowrate = spikeFlowrate;

%% initial conditions 
xInitial=get_initial_conditions(parameters,reactor,'Generic');

%% Reed data 
[Reed] = get_reed_data('Generic');
% give to structure parameters 
parameters.Reed = Reed; % extraction is more or less the mean of each experiment also cte 
ReedSwitch = 0; % not using percentage of extraction: turn switch off 


%% run ode45 to check simulation 
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
%% decisionVar initial values 
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


%% optimisation algorithms Patternsearch

% define disered options 
%options = optimoptions('patternsearch','Display','iter','PlotFcn',{'psplotbestf', 'psplotmeshsize'});
options = optimoptions('patternsearch','Display','iter','TolFun',1e-3,'TolX',1e-3); % skip plot func


switch objective
    case 'yield'
        % YIELD 
        decisionVar = [ 167.6723    4.5953   18.5137    2.0500    2.7480 2.0500    2.7480    3.1179    1.7875 ]; % so it doesn't go to other local minima!! (from pareto front)
        costFunc = @(decisionVar)f_objFunc_stochastic(decisionVar, parameters, xInitial,objective);
        tic
        [decisionVarOutYield,fValYield] = patternsearch(costFunc,decisionVar,[],[],[],[],lb,ub,[],options);
        toc
        %% get the distribution of yield
        objective = 'yield';
        [~,yieldDist] = f_objFunc_stochastic(decisionVarOutYield,parameters,xInitial,objective);
        %% save
        save(saveName,'decisionVarOutYield','yieldDist','parameters','xInitial')
    case 'productivity'
        % PRODUCTIVTY
        %decisionVar = [ 206    6.7109   10.2164    2.0500    0.1000    2.0500    0.1000    2.0500    0.1000 ]; % so it doesn't go to other local minima!! (from pareto front)
        decisionVar=  [206.0000    7.5000   11.0312    2.7500    0.2500    2.2500    0.1870    2.0625    0.1251]; % this should be close to the reult already
        costFunc = @(decisionVar)f_objFunc_stochastic(decisionVar, parameters, xInitial,objective);
        tic
        [decisionVarOutProd,fvalProd] = patternsearch(costFunc,decisionVar,[],[],[],[],lb,ub,[],options);
        toc
        %% get the distribution of yield and productivity
        % One problem how ever... is there possibilty of getting last calculated
        % distributions (yield and prod) off of patternsearch? because distribution won't be the same
        % due to random sampling during monte carlo simulation with underlying code.... maybe somthing with rng?
        objective ='productivity';
        [~,~,productivityDist] = f_objFunc_stochastic(decisionVarOutProd,parameters,xInitial,objective);
        %% save
        save(saveName,'decisionVarOutProd','productivityDist','parameters','xInitial')
end




