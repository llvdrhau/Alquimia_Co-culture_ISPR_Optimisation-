%Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% May 2020. Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.
% Run a generic Co-culture simulation with f_max_productivity2 


clear  
close all 
clc

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
parameters.Reed = Reed; % extraction is more or less the mean of each experiment also cte 

ReedSwitch = 0;
parameters.ReedSwitch = ReedSwitch;
if ReedSwitch
parameters.Reed.RateMatrix(:,2) = 0.30;  %lactate
parameters.Reed.RateMatrix(:,4) = 0.30;  % proprionate 
parameters.Reed.RateMatrix(:,5) = 0.30;  % Acetate  
end


%% run ode15s
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
spaces = ones(1,nSpace)*5;
gluConc = 400; %gCOD/L

gluON = 1; 
parameters.gluSwitch = gluON; 
if gluON
    decisionVar = [gluConc, spikeVeil, start, spaces];
else
    decisionVar = [spikeVeil, start, spaces];
end 


%% boundries 
lbSpace = ones(1,nSpace)*2.1;
ubSpace =  ones(1,nSpace)*20;
if gluON
    %lb = [50, 0, 5, lbSpace];
    lb = [50, 0, 0.1, lbSpace]; % change when reed can start, not sure if effective? 
    ub = [1000, 20, 20, ubSpace];
else
    lb = [0, 5, lbSpace];
    ub = [20, 20, ubSpace];
end


%% contour plots of f_max_yield to see a rough estimate of the optimum location 
for e =1
vx = linspace(lb(1),ub(1),7);
vy = linspace(lb(2),ub(2),7);
vz = linspace(lb(3),ub(3),7);
%list = [vx vy]

X = repmat(vx,length(vx),1);
Y = repmat(vy,length(vy),1);
Z = repmat(vz,length(vz),1);
 contureOn = false; 
if contureOn
    [X,Y] = ndgrid(vx,vy);
     [nx,ny]=size(X);
    
    for i=1:nx
        for j=1:ny
            parEst(1)=X(i,j);
            parEst(2) = Y(i,j);
            decisionVar = [parEst,3,2,2]
            [residual] = f_max_yield(decisionVar,parameters,xInitial)
            resNorm(i,j) = norm(residual);
        end
    end
    figure
    [c,h] = contour(X,Y,resNorm);
    clabel(c,h)
    
end
 
%% 3D contour plots of f_max_yield to see a rough estimate of the optimum location 

 contureOn3D = false; 
 if contureOn3D
     vx = linspace(lb(1),ub(1),7);
     vy = linspace(lb(2),ub(2),7);
     vz = linspace(lb(3),ub(3),7);
     
     [X,Y,Z] = meshgrid(vx,vy,vz);
     [nx,ny,nz]=size(X);
     
     for i=1:nx
         for j=1:ny
             for k = 1:nz
                 parEst(1)=X(i,j,k);
                 parEst(2) = Y(i,j,k);
                 parEst(3) = Z(i,j,k);
                 decisionVar = [parEst,2,2];
                 [residual] = f_max_yield(decisionVar,parameters,xInitial)
                 resNorm(i,j,k) = residual;
             end
         end
     end
  save('3D','resNorm')
 end
 
end


%% Multi objective solution (yield and productivity) with patern search algorithm
% % original multiobjective function
% tic
% rng default % For reproducibility
% costFunc = @(decisionVar)multiObjectiveFunc(decisionVar, parameters, xInitial);
% nVar = length(lb);
% %options = optimoptions('paretosearch','ParetoSetSize',200);
% options = optimoptions('paretosearch','ParetoSetSize',200 ,'PlotFcn','psplotparetof'); % large ParetoSize 
% [decisionVarOut,fxVal]= paretosearch(costFunc,nVar,[],[],[],[],lb,ub,[],options);
% toc
% 
%% Genetic search followed by patternsearch 
% tic
% costFunc = @(decisionVar)multiObjectiveFunc(decisionVar, parameters, xInitial);
% nVar = length(lb);
% options = optimoptions('gamultiobj','PlotFcn','gaplotpareto'); %  use this on parteo front to get more accurate  solutions 
% [decisionVarOut,fxVal]= gamultiobj(costFunc,nVar,[],[],[],[],lb,ub,[],options);
% toc

% Start patreto from Genetic 
% load genetic reults 
saveName = 'paretosearchGenetic';
load(saveName,'decisionVarOut','fxVal','parameters','xInitial')
decisionVarGenetic = decisionVarOut;
% run opt 
tic
costFunc = @(decisionVar)multiObjectiveFunc(decisionVar, parameters, xInitial);
nVar = length(lb);
options = optimoptions('paretosearch','ParetoSetSize',200,'InitialPoints',decisionVarGenetic); %  use this on parteo front to get more accurate  solutions 
[decisionVarOut,fxVal]= paretosearch(costFunc,nVar,[],[],[],[],lb,ub,[],options);
toc 

%% Save multiObj solutions  
% saveName = 'paretosearchStartingFromGenetic2';
saveName = 'paretosearchUpdated2';
save(saveName,'decisionVarOut','fxVal','parameters','xInitial')


%% plot 
% decisionVarOut
% title = 'Optimised generic extraction process';
% %[~,productivity,parametersOut,t,y] = f_max_productivity2(decisionVarOut,parameters,xInitial);
% [~,productivity,parametersOut,t,y] = f_max_yield(decisionVarOut,parameters,xInitial);
% plot_CoCulture_simulation(t,y,parametersOut,title)
% beep;
% pause(0.5)
% beep; 
% pause(0.5)
% beep; 

%disp(productivity)

%% Save intersting solutions 
%  saveName = 'patternYield4spikes65prod';
%  save(saveName,'decisionVarOut','parameters','xInitial')

% 






