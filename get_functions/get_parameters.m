function parameters=get_parameters(reactor,flagCase,kineticFlag)
%% Definition of model options and parameters
% Miguel Mauricio Iglesias. University of Santiago de Compostela. Spain
% November 2014. Please contact miguel.mauricio.iglesias@gmail.com if you
% intend to use this code.

%% Load the compounds used by each model

% Index guide
% M     Total Mass              1   kg m-3
% Ssu	Monosaccharides         2	kg COD m-3
% Slac  Total lactate           8   kg COD m-3
% Sye   Yeast extract
% Spro	Total propionate        7	kg COD m-3
% Sac	Total acetate           10	kg COD m-3
% Sh2	Hydrogen                11	kg COD m-3
% Stic	Inorganic carbon        13	kmole C m-3
% SIn	Inorganic nitrogen      14	kmole N m-3
% Xc	composites              16	kg COD m-3
% Xsu	Sugar degraders         20	kg COD m-3
% Xlac  Lactate degraders       25
% Scat	Cations                 30	kmole m-3
% San	Anions                  31	kmole m-3

% add flag case 

parameters.flagCase = flagCase;
parameters.kineticFlag = kineticFlag;

[~,~,compoundTable] = xlsread('dataFile.xlsx','Compounds');

compoundAbb = compoundTable(2:end,1);
compounds = compoundTable(2:end,2);
indexCompounds = compoundTable(2:end,3);
compoundUnits = compoundTable(2:end,4);
vectorMolarMass = compoundTable(2:end,5);
vectorFlagLiquid = compoundTable(2:end,7);
vectorFlagGas = compoundTable(2:end,8);

%Convert into numeric (double) variables
indexCompounds = [indexCompounds{:}]';
vectorFlagLiquid = [vectorFlagLiquid{:}]';
vectorFlagGas = [vectorFlagGas{:}]';
vectorMolarMass = [vectorMolarMass{:}]';

% Create the structure parameters. Compounds, compound indexes, liquid or
% gas flags and MW.
parameters.compounds = compounds;
parameters.indexCompounds = indexCompounds;
parameters.compoundAbb = compoundAbb;
parameters.compoundUnits = compoundUnits;
parameters.flagLiquid = logical(vectorFlagLiquid(2:end));
parameters.flagGas = logical(vectorFlagGas(2:end));
parameters.Mw = vectorMolarMass;

switch flagCase
    case 'CoCulture'
        % Get parameter table
        [~,~,parameterTable] = xlsread('dataFile.xlsx','parametersCoCulture');
        %Set the stoichiometric matrix
        [~,~,stoTable] = xlsread('dataFile.xlsx','StoichiometryCoCulture');
    case 'SpikeVeillonella'
         % Get parameter table
        [~,~,parameterTable] = xlsread('dataFile.xlsx','parametersCoCulture');
        %Set the stoichiometric matrix
        [~,~,stoTable] = xlsread('dataFile.xlsx','StoichiometryCoCulture');
        
    otherwise
        % Get parameter table
        [~,~,parameterTable] = xlsread('dataFile.xlsx','Parameters');
        %Set the stoichiometric matrix
        [~,~,stoTable] = xlsread('dataFile.xlsx','Stoichiometry');
end


%load parametersCalibratedBiofilm.mat
parAbb = parameterTable(2:end,1);       %Read the parameters abbreviations
parValues = parameterTable(2:end,3);    %Read the parameters values
parValues = [parValues{:}];             %Convert into a vector

parameters.parAbb = parAbb;
parameters.parValues = parValues;


%% Calculating fva.aa, fbu.aa, fpro.aa and fac.aa for the protein given

stoM = f_stoichiometric_matrix(stoTable,parameters);
parameters.stoM=stoM;
parameters.stoTable = stoTable;

%% Load the constant values and values needed for simulation
R = 8.314e3; %J/K/kmol
T = 273.15 + 25; %K
constants=[R T];
parameters.constants=constants;
nCompounds=length(indexCompounds);
parameters.simulation=nCompounds;
parameters.nStates=nCompounds;

% %Gas phase parameters
% Vgas = reactor.Vgas;     % m3
% Patm  =1.013e5;   % Pa
% %Ph2o = 1.2e5;                   %Pa

% gasPhase=[Vgas Patm];
% parameters.gasPhase = gasPhase;
%% Load the different characteristics of the model
%if any(strcmp(modelCharacteristics,'pH1'))
% acid base constants and flags for pH calculation
%% Load the different characteristics of the model

% acid base constants and flags for pH calculation
pHconstant=reactor.pH(1);           %Boolean variable. 1=Yes pH is constant and pH calculation is overridden by a given value of pH; 0=No, pH is determined at every t
pHValue=reactor.pH(2);              %If pH is not determined
parameters.pH.pHValue=pHValue;
parameters.pH.pHguess=pHValue;
parameters.pH.pHconstant=pHconstant;



