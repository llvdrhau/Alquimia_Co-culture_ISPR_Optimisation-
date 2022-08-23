%% Kinetic model of anaerobic fermentation with optional methanogenesis
% Alberte Regueira. University of Santiago de Compostela. Spain
% November 2019. Please contact alberte.regueira@gmail.com if you
% intend to use this code. 

clear
close all
path =genpath( cd());
addpath(path)

reactor.V=1;            % L
reactor.Vgas=1;         % L
%reactor.TRH=1;          % d
reactor.pH(1)=1;        % Define whether pH is kept constant or not
reactor.pH(2)=7;        % pH value if it's controlled


%% Get parameters, initial conditions and feed flow(s)

%Get parameters and initial conditions
parameters=get_parameters(reactor);


%Get initial conditions
xInitial=get_initial_conditions(parameters,flagProteins,reactor);

%Get inflow 
%inputFeed = get_inflow(parameters,reactor);

%Simulation parameters
tEnd=100;                        % h
nIteration=min(5e2,tEnd*(6));    % one time node per 10 min
t=linspace(0,tEnd,nIteration);

%% Run simulation
tic
options = odeset('RelTol',1e-7,'AbsTol',1e-8);
[t,y]=ode15s(@f_mass_balances,t,xInitial,options,parameters);
toc


%% Plot and save results
plot_reactor_results(t,y, parameters);

FIN.parameters = parameters;
FIN.t = t;
FIN.y = y;
save batch_pH_7_xel.mat FIN
