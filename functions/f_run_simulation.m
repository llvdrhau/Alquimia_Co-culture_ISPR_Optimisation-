function [t,y,parameters] = f_run_simulation(parameters,sheetData,control,xInitial)
% Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% October 2021.Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.
% Run a co_culture simulation 

parameters.control = control; 
%% spike information BV5
if control(2)
    spikeFlowrate = get_spike_data(sheetData);
    parameters.spikeFlowrate = spikeFlowrate;
end

%% Reed data 
if control(1)
    [Reed] = get_reed_data(sheetData);
    parameters.Reed = Reed; % give to structure parameters
end

%% edit the extraction process with mean removal rate  
meanExtractionRate = 1;
if meanExtractionRate && control(1)
    ReedRateMatrix = Reed.RateMatrix;
    for i = 1:size(ReedRateMatrix,2)    % TODO check if this changes each row to it's corresponding mean rate
        parameters.Reed.RateMatrix(:,i) = mean(ReedRateMatrix(:,i));  %lactate
    end
end

%% run ode15s
tsim = linspace(0,parameters.tEnd,100);
tSpike =  spikeFlowrate.tSpikes;
tOde = union(tsim,tSpike); 
% optioins  
nneg = 1:length(xInitial);
options=odeset('RelTol',1e-4,'NonNegative',nneg,'MaxStep',0.1);

[t,y] = ode15s(@f_mass_balances,tOde,xInitial,options,parameters);


end

