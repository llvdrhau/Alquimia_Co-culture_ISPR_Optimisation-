function  [dxdt,pH,transferReed,Ivfa,growthRateLac]=f_mass_balances(t,states,parameters)
% Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% October 2021.Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.
% Alberte Regueira. University of Santiago de Compostela. Spain
% November 2019. Please contact alberte.regueira@gmail.com if you
% intend to use this code.
%%  Generic mass balance of a well stirred tank reactor
% Input: states (mass, composition) and feeding characteristics (flow, concentrations)
% Output: Derivates to integrate (mass balances) and pH value

M=states(1);            %Total mass (kg)
V = M; %assume rho = 1 kg/L
C=states(2:end);        %Concentration (kg/m3)

% retrive control options 
control = parameters.control; 


%% Feeding. To use when in a CSTR or fed-batch
% if isstruct(input)
%     %Inflow (1st component). A structure is used if the influent varies
%     %with time
%     mInflow=input.mInflow;          %m3/s
%     tFin=input.tFin;
%     %Define Fin and the inflow concentration
%     Fout=0;
%     Fin= interp1(tFin,mInflow(:,1),t);
%     inflowComposition= interp1(tFin,mInflow(:,2:end),t)';
% else
%     %If influent is constant
%     Fin=input(1);
%     inflowComposition = input(2:end);
%     Fout=Fin;
% end

flagLiquid = parameters.flagLiquid;
flagGas = parameters.flagGas;
% D = states(2:end)*0;
% D(flagLiquid) = 0;

if any(flagGas)
    [Qgas,Vgas]=f_gasPhase(states,parameters);
    D(flagGas) = Qgas/Vgas;
end

%% Add spikes
spikesActive = control(2);
if spikesActive 
    %inSpikes = C*0;
    [Fspikes,inSpikes]=f_spikes_flowrate(t, parameters);
else
    Fspikes = 0;
    inSpikes = C*0;
end

 %% cummulative inputs
Fin = Fspikes;
inflowComposition = (Fspikes*inSpikes)/(Fin+1e-9);  %the sum of every thing that flows in spikes
D=Fin/V;                    % dilution rate in the liquid phase (h-1)

%% get reaction rates
%if isfield(parameters,'ReedSwitch') && parameters.ReedSwitch == 1 
if parameters.control(1)
    [rt, pH,transferReed,Ivfa,growthRateLac]=f_rt(states,parameters,t);  % Get the vector of rates + transfer
    r=rt(2:end);                         % The first component is the total mass, the rest are the components
else
    [rt, pH]=f_rt(states,parameters,t);  % Get the vector of rates + transfer
    r=rt(2:end);  % The first component is the total mass, the rest are the components
end

%% Total mass balance in the liquid phase
dMdt=Fin+0;

% lag on reaction 
lag_time = parameters.parValues(strcmp(parameters.parAbb,'lag_time'));
factor = (tanh(10*(t-lag_time))+1)/2;
r = r*factor;

% Component balances
dCdt= D.*inflowComposition - D.*C + r;

dxdt=[dMdt;dCdt];
