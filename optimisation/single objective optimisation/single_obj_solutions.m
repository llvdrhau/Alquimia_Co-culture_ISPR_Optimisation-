% Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% October 2021.Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.

% Runs the saved solutions from the single objective optimisation and plots
% the results 

clear  
%close all 
%clc 

%% call the saved solution
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











