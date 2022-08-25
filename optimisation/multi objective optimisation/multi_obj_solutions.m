%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualisation of Solution from the run_multi_obj_optimisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% October 2021.Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.

% Runs the saved solutions from the multi objective optimisation and plots
% the results 

clear 
close all 
clc

%% Load saved result 
saveName = 'paretosearchOpt2_theOne';
load(saveName,'decisionVarOut','fxVal','parameters','xInitial')

%% Visulise paretoFront
decisionVarOutAll = decisionVarOut;
fxValAll = fxVal;
close all


%sort the solutions based on yield
[A,I] = sort(fxValAll(:,2)); % sort by yield
allSorted = [fxValAll(I,:),decisionVarOutAll(I,:)];
yield = -allSorted(:,2);
prod = -allSorted(:,1);

% pareto front
figure
set(gcf,'color','w'); %background color
plot(prod(7:end),yield(7:end),'k*') % before solution 7: Unreqlistic solutions 

hold on 
% highlight intersting soolution 
plot(prod(end),yield(end),'gd', 'MarkerSize',15,'MarkerFaceColor','g') % productivity result
plot(prod(7),yield(7),'rd', 'MarkerSize',15,'MarkerFaceColor','r') % Yield result 
plot(prod(71),yield(71),'bd', 'MarkerSize',15,'MarkerFaceColor','b')% example used for trade-off
xlabel('Productivity (gCOD/L/h)','FontSize', 16)
ylabel('Yield (gCOD/gCOD)','FontSize', 16)
legend('','Maximizing the productivity','Maximizing the yield', 'Example of trade-off solution','FontSize', 16)

%% Extra plot 

figure
set(gcf,'color','w'); %background color

X = prod(7:end);
y = yield(7:end);

mdl = fitlm(X,y);
out = predict(mdl,X);


plot(X,out,'k','LineWidth',3)
%plot(prod(7:end),yield(7:end),'k*') % before solution 7: Unreqlistic solutions 

hold on 
% highlight intersting soolution 
plot(prod(end),yield(end),'gd', 'MarkerSize',15,'MarkerFaceColor','g') % productivity result
plot(prod(7),yield(7),'rd', 'MarkerSize',15,'MarkerFaceColor','r') % Yield result 
plot(prod(71),yield(71),'bd', 'MarkerSize',15,'MarkerFaceColor','b')% example used for trade-off
plot(0.9915,0.6382,'bd', 'MarkerSize',15,'MarkerFaceColor','y')% reultfrom the stochastic optimisation
xlabel('Productivity (gCOD/L/h)','FontSize', 16)
ylabel('Yield (gCOD/gCOD)','FontSize', 16)
legend('Aprox. Pareto Front','Maximizing the productivity','Maximizing the yield', 'Example of trade-off solution','Max Productivity (robust)' ,'FontSize', 16)
grid on



% 1st stragetgy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% max productivity 

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

%% 2nd strategy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% trade off solution 
x2 = allSorted(71,3:end); % solution 71 used in paper. 

% plot strategy
[costFunc,parametersOut,t,y] = multiObjectiveFunc(x2,parameters,xInitial);
plot_Strategies(t,y,parametersOut)

disp('productivity')
disp(costFunc(1))
disp('Yield')
disp(costFunc(2))
disp('descision varibles')
disp(x2)
%% 3d strategy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% max yield 
 
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


%% %%%%%%%%%%%%%%%%%% 
% This part was to check out the differential error from the transfer
% function (not really necesary any more)
manipulate = 0;
for l = 1
    if manipulate
        xman = allSorted(5,3:end);
        %xman =  decisionVarOutAll(yield == min(yield),:);
        % make plot
        titlePlot = 'max yield manipulated';
        
        % manipulate to upgrade the prod
         xman(2)= 2;
         %xman(4:end) = 4;
%         xman(6) = 4;
%         xman(end) = 4;  %4 is goed vb
        % x2(3) = 18.2;
        
        [costFunc,parametersOut,t,y,reservoirMan] = multiObjectiveFunc(xman,parameters,xInitial);
        figure
        plot_CoCulture_simulation(t,y,parametersOut,titlePlot)
        
        
        disp('productivity')
        disp(costFunc(1))
        disp('Yield')
        disp(costFunc(2))
        disp('descision varibles')
        disp(xman)
        
        % rename some variables
        yMan  = y;
        tMan = t;
        parMan = parametersOut;
        save('timeMan','tMan','parMan')
        
        %% lets see wat's going on
        compareReservoir = 0;
        if compareReservoir
            yOg = [y1(:,2:6),y1(:,11:12)].*y1(:,1); % convert to gCOD ipv gCOD/L
            SetOg = [yOg,reservoirOg];
            
            % IF conversion of time spteps is used (not in this case )
            % timeOgReservoir = 0:0.1:tOg(end);
            % reservoirOgConvert = interp1(timeOgReservoir,reservoirOg,tOg);
            % SetOg = [yOg,reservoirOgConvert];
            
            yManual = [yMan(:,2:6),yMan(:,11:12)].*yMan(:,1); % convert to gCOD ipv gCOD/L
            SetMan = [yManual,reservoirMan];
            
            % timeManReservoir = 0:0.1:tMan(end);
            % reservoirManConvert = interp1(timeManReservoir,reservoirMan,tMan);
            % SetMan = [yManual,reservoirManConvert];
            
            balanceOg = [sum(SetOg,2),tOg];
            balanceMan = [sum(SetMan,2),tMan];
            
            balanceCompare = [balanceOg,balanceMan];
            
            difference = balanceCompare(end,1)- balanceCompare(end,3)
            differenceReservoir = sum(reservoirMan(end,:)) - sum(reservoirOg(end,:))
            
            %close all
            
            plot(tOg,reservoirOg,tMan,reservoirMan )
        end
    end
end
