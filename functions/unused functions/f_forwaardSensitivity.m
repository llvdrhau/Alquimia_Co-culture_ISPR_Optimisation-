function [sensF, sensB] = f_forwaardSensitivity(parameters,par_name, par_value,xInitial,control,sheetData,epsilon)
%explanation: calculates the absolute sensitivity (forward and backward) of a parameter to
%all output parameters. 
%detailed explanation


flagCase = parameters.flagCase;
%typeExperiment = parameters.typeExperiment;

if nargin < 7
    epsilon = 1e-7;
end

originalParameters = parameters;          
[time,y1] = f_run_simulation(parameters,sheetData,control,xInitial);

            
%% Forward 

% change parameters with perturbation epsilon 
pos = find(strcmp(parameters.parAbb,par_name));
parameters.parValues(pos) = par_value *(1 + epsilon);

% !!! update Stochiometric matrix !!! 
stoTable = parameters.stoTable;
stoM = f_stoichiometric_matrix(stoTable,parameters);   %links the parameters                                                    %(e.g Yxs) 
parameters.stoM=stoM;                                   % with either actual value e.g. Yxs = 0,6  
                                                         % and caculates the stoichimetric coefficients
                                                          % like in the excel sheet stoichiometryModel1

% run a nieuw simulation 
[~,y2] = f_run_simulation(parameters,sheetData,control,xInitial);

% Forward sensitivity 
sensyF = (y2 - y1)/(epsilon*par_value);  %LET OP parValue is origineel parameter 

%% Backward

% update parameters
parameters = originalParameters;
parameters.parValues(pos) = par_value *(1 - epsilon);
% !!! update Stochiometric matrix !!! 
stoTable = parameters.stoTable;
stoM = f_stoichiometric_matrix(stoTable,parameters);   %links the parameters                                                    %(e.g Yxs) 
parameters.stoM=stoM;                                   % with either actual value e.g. Yxs = 0,6  
                                                         % and caculates the stoichimetric coefficients
                                                          % like in the excel sheet stoichiometryModel1
% run a nieuw simulation 
[~,y3] = f_run_simulation(parameters,sheetData,control,xInitial);


% backwards sensitivity
sensyB = (y1 - y3)/(epsilon*par_value);  %LET OP parValue is origineel parameter 

%%
% add time array to the sensitivity
sensF = [time sensyF];
sensB = [time sensyB];

% originele waarde param opnieuw toekennen aan base workspace!!!
parameters = originalParameters;
end

