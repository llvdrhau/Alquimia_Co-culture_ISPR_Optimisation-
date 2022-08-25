function [r,pH,transferReed,Ivfa,growthRateLac]=f_rt(states,parameters,t)
% Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% October 2021.Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.
%   Input: states, parameters and time
%   Output: reaction term, REED transfer, current inhibitio and growthrtate
%   Veillionella 


%states(states<1e-7) = 0;                   % Erase negative states
pHguess = parameters.pH.pHguess;
pH = f_pH(states,parameters,pHguess,t); 
parameters.pH.pHValue = pH;

% select the right kinetic model
switch parameters.kineticFlag
    case 'V_inhibition1'
        reac = f_reaction_V_inhibition1(states,parameters);
    case 'V_inhibition2'
        reac = f_reaction_V_inhibition2(states,parameters);
    case 'V_inhibition3'
        reac = f_reaction_V_inhibition3(states,parameters);
    case 'V_inhibition4'
        reac = f_reaction_V_inhibition4(states,parameters);
    case 'CoCulture'
        [reac,Ivfa,growthRateLac] = f_reaction_co_culture(states,parameters);
    otherwise  % for batch experiments
        reac = f_reaction(states,parameters);
end

% Select manner of REED extraction 
if parameters.control(1)
    if isfield(parameters,'ReedSwitch') && parameters.ReedSwitch == 1
        transferReed=f_REED_percent(states,parameters,t);
        r = reac - transferReed;
    else
        transferReed=f_REED(states,parameters,t);
        r = reac - transferReed;
    end
    
else
    r = reac;
end
end
