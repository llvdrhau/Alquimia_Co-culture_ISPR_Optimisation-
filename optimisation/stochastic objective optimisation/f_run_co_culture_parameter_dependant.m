function [t,y] = f_run_co_culture_parameter_dependant(tOde,xInitial,parameters)
% Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% October 2021.Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code..
%   Run a simulation of a co_culture experiment 

%% update parameter values 
%done in script f_montecarlo_here_and_now

%% update initial conditions 
%done in script f_montecarlo_here_and_now
%% run ode15s; 
% optioins  
nneg = 1:length(xInitial);
options=odeset('RelTol',1e-4,'NonNegative',nneg,'MaxStep',0.1);
[t,y] = ode45(@f_mass_balances,tOde,xInitial,options,parameters);




end

