% Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% October 2021.Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.

% Runs the saved solutions from the single objective optimisation and plots
% the results 

clear  
close all 
%clc 

%% Load saved result 
saveName = 'paretosearchOpt2_theOne';
load(saveName,'decisionVarOut','fxVal','parameters','xInitial')

decisionVarOutAll = decisionVarOut;
fxValAll = fxVal;

%sort the solutions based on yield
[A,I] = sort(fxValAll(:,2)); % sort by yield
allSorted = [fxValAll(I,:),decisionVarOutAll(I,:)];
yield = -allSorted(:,2);
prod = -allSorted(:,1);
%% extreem ends of the pareto front 
% 1st stragetgy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% max productivity 

x1 = allSorted(end,3:end); % solution 71 used in paper. 
% make plot 
[costFunc,parametersOut,t,y] = multiObjectiveFunc(x1,parameters,xInitial);
plot_Strategies(t,y,parametersOut)

disp('productivity')
disp(costFunc(1))
disp('Yield')
disp(costFunc(2))
disp('descision varibles')
disp(x1)

% 2nd stragetgy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% max yield 
 
x3 = allSorted(7,3:end); % decision variables 
% not the max yield becase this simulation is unnecesarlily streched out 
% make plot 
[costFunc,parametersOut,t,y,reservoirOg] = multiObjectiveFunc(x3,parameters,xInitial);
plot_Strategies(t,y,parametersOut)

disp('productivity')
disp(costFunc(1))
disp('Yield')
disp(costFunc(2))
disp('descision varibles')
disp(x3)

%% call the saved solution from SOO algorithm and compare 
saveName = 'SOO_3';
load(saveName,'decisionVarOutYield','parameters','xInitial')
decisionVarOut = decisionVarOutYield;
%% using pareto search cost function
[costFunc,parametersOut1,t1,y1,reservoirOG] = multiObjectiveFunc(decisionVarOut,parameters,xInitial);

%% plot
figure
titlePlot = 'pattern';
plot_CoCulture_simulation(t1,y1,parametersOut1,titlePlot)

disp('productivity')
disp(costFunc(1))
disp('Yield')
disp(costFunc(2))
disp('descision varibles')
disp(decisionVarOut)











