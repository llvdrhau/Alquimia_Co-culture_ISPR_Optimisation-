%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualisation of Solution from the here and now script 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% October 2021.Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.

clear  
close all 
clc 
%% pattern Solutions  

% here and now interseting solutions 
% saveName = 'here_and_now_constrained';  % Average prod
% saveName = 'here_and_now_constrained2'; % Average prod (higher constraint)
% saveName = 'here_and_now_percentil_10'; % 90th percentile
% saveName = 'here_and_now_percentil_90'; % 10th percentile 

%% Soltion average productivity 
saveName = 'here_and_now_constrained2'; % objective variable: average distr
load(saveName,'decisionVarOutYield','decisionVarOutProd','parameters',...
    'xInitial','yieldDist','productivityDist')
decisionVarOut = decisionVarOutProd;
% get simulation outputs 
[costFunc,parametersOut,t1,y1] = multiObjectiveFunc(decisionVarOut,parameters,xInitial);

% plot
figure
histogram(productivityDist)
figure
plot_Strategies(t1,y1,parametersOut)

disp('AVERAGE solution')
disp('productivity')
disp(costFunc(1))
disp('Yield')
disp(costFunc(2))
disp('descision varibles')
disp(decisionVarOut)

figure
histogram(productivityDist)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10TH PERCENTILE SOLUTION
saveName = 'here_and_now_percentil_90_200MC'; %interseting solution
load(saveName,'decisionVarOutYield','decisionVarOutProd','parameters',...
    'xInitial','yieldDist','productivityDist')
decisionVarOut = decisionVarOutProd;

% using pareto search cost function
[costFunc,parametersOut,t1,y1,reservoirOG] = multiObjectiveFunc(decisionVarOut,parameters,xInitial);
% plot
figure
histogram(productivityDist)
plot_Strategies(t1,y1,parametersOut)

disp('10th PERCENTILE solution')
disp('productivity')
disp(costFunc(1))
disp(mean(productivityDist))
disp('Yield')
disp(costFunc(2))
%disp(mean(yieldDist))
disp('descision varibles')
disp(decisionVarOut)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 90TH PERCENTILE SOLUTION
saveName = 'here_and_now_percentil_10_200MC'; 
load(saveName,'decisionVarOutYield','decisionVarOutProd','parameters',...
    'xInitial','yieldDist','productivityDist')
decisionVarOut = decisionVarOutProd;

% using pareto search cost function
[costFunc,parametersOut,t1,y1,reservoirOG] = multiObjectiveFunc(decisionVarOut,parameters,xInitial);
% plot
figure
histogram(productivityDist)
plot_Strategies(t1,y1,parametersOut)

disp('90th PERCENTILE solution')
disp('productivity')
disp(costFunc(1))
disp('average solution from the distribution of solutions')
disp(mean(productivityDist))
disp('Yield')
disp(costFunc(2))
%disp(mean(yieldDist))
disp('descision varibles')
disp(decisionVarOut)
