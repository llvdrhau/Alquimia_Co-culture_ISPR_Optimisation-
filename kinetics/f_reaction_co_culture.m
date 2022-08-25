function  [r,Ivfa,growthRateLac]=f_reaction_co_culture(states,parameters)

% kinetic equations and selection of inhibition model 

% Alberte Regueira & Lucas Van der Hauwaert. University of Santiago 
% de Compostela. Spain
% october 2020. Please contact alberte.regueira@gmail.com or lucas.vanderhauwaert@usc.es 
% if you intend to use this code.

r=0*states;                         % Initialize rates vector
parValues = parameters.parValues;
parAbb = parameters.parAbb;
compoundAbb=parameters.compoundAbb;
stoM=parameters.stoM;
pH = parameters.pH.pHValue;
H=10^(-pH);

%% Identify compounds
Sglu = states(strcmp(compoundAbb,'Sglu'));       % Glucose
Slac = states(strcmp(compoundAbb,'Slac'));       % Lactate
Sye = states(strcmp(compoundAbb,'Sye'));        % Yeast extract
Spro = states(strcmp(compoundAbb,'Spro'));       % Propionate
Sac = states(strcmp(compoundAbb,'Sac'));         % Acetate
Sh2 = states(strcmp(compoundAbb,'Sh2'));         % H2
Stic = states(strcmp(compoundAbb,'Stic'));       % TIC
SIn = states(strcmp(compoundAbb,'SIn'));         % TAN
Xc = states(strcmp(compoundAbb,'Xc'));           % Composites
Xsu = states(strcmp(compoundAbb,'Xsu'));         % Sugar degraders
Xlac = states(strcmp(compoundAbb,'Xlac'));       % Lactate degraders

%% Get parameter values
KS_IN = parValues(strcmp(parAbb,'KS_IN'));
km_su = parValues(strcmp(parAbb,'km_su'));
KS_su = parValues(strcmp(parAbb,'KS_su'));
km_lac = parValues(strcmp(parAbb,'km_lac'));
KS_lac = parValues(strcmp(parAbb,'KS_lac'));
km_ye = parValues(strcmp(parAbb,'km_ye'));
KS_ye = parValues(strcmp(parAbb,'KS_ye'));
pHUL_aa = parValues(strcmp(parAbb,'pHUL_aa'));
pHLL_aa = parValues(strcmp(parAbb,'pHLL_aa'));
KI_NH3 = parValues(strcmp(parAbb,'KI_NH3'));
kdec_Xsu = parValues(strcmp(parAbb,'kdec_Xsu'));
kdec_Xlac = parValues(strcmp(parAbb,'kdec_Xlac'));
Kac_in = parValues(strcmp(parAbb,'Kac_in'));
KI_VFA = parValues(strcmp(parAbb,'KI_VFA'));
Kac_ac = parValues(strcmp(parAbb,'Kac_ac'));
Kac_pro = parValues(strcmp(parAbb,'Kac_pro'));
Kac_bu = parValues(strcmp(parAbb,'Kac_bu'));
Kac_lac = parValues(strcmp(parAbb,'Kac_lac'));

%nieuw additions
Yyeb = parValues(strcmp(parAbb,'Yyeb'));
Yyev = parValues(strcmp(parAbb,'Yyev'));
KS_yeb = parValues(strcmp(parAbb,'KS_yeb'));
KS_yev = parValues(strcmp(parAbb,'KS_yev'));

% inhibition cte
KI_lac = parValues(strcmp(parAbb,'KI_lac'));
KI_pro = parValues(strcmp(parAbb,'KI_pro'));

%% Inhibition terms

% 1)Sum of inhibition
% TVFA = Sac + Spro + Slac;
% %TVFA =  Spro + Slac;
% Ivfa= 1/(1+(TVFA/KI_VFA));

% 2)separate contribution of inhibition
if isfield(parameters,'kineticInhibition')
    switch parameters.kineticInhibition
        case 'totalVFA'
            TVFA = Sac + Spro + Slac;
            %TVFA =  Spro + Slac;
            Ivfa= 1/(1+(TVFA/KI_VFA));
        case 'seperateProLac'
            Ivfa=1/(1+(Slac/KI_lac)) * 1/(1+(Spro/KI_pro)) ;
        case 'proOnly'
            Ivfa= 1/(1+(Spro/KI_pro)) ;
        case 'lacOnly'
            Ivfa=1/(1+(Slac/KI_lac)) ;
     end
else
    Ivfa=1/(1+(Slac/KI_lac)) * 1/(1+(Spro/KI_pro)) ;  % set this as default
end


%Yeast inhibition 
%Iye = 1/(1+(Slac/0.0015));
Iye = 1;

%% Process rates

% Glucose uptake (rho1)
rho(1,1) = km_su*(Sglu/(KS_su+Sglu))*Xsu;

% Yeast extract uptake Bacillus (rho2)
rho(2,1) = km_su*(Sye/(KS_yeb+Sye))*Xsu; %*Ivfa;

% Yeast extract uptake Veillionella (rho3)
rho(3,1) = km_lac*(Sye/(KS_yev+Sye))*Xlac*Ivfa*Iye;

% Lactate uptake (rho4)
rho(4,1) = km_lac*(Slac/(KS_lac+Slac))*Xlac*Ivfa;

% Sugar degraders decay (rho5)
rho(5,1) = kdec_Xsu*Xsu;

% Lactate degraders decay (rho6)
rho(6,1) = kdec_Xlac*Xlac;

rADM=stoM'*rho;
 
% %Rate assignment
r(strcmp(compoundAbb,'Sglu'))=rADM(1);           % Glucose
r(strcmp(compoundAbb,'Slac'))=rADM(2);           % Lactate
r(strcmp(compoundAbb,'Sye'))=rADM(3);            % Yeast extract
r(strcmp(compoundAbb,'Spro'))=rADM(4);           % Propionate
r(strcmp(compoundAbb,'Sac'))=rADM(5);            % Acetate
r(strcmp(compoundAbb,'Sh2'))=rADM(6);            % H2
r(strcmp(compoundAbb,'Stic'))=rADM(7);           % TIC
r(strcmp(compoundAbb,'SIn'))=rADM(8);            % TAN
r(strcmp(compoundAbb,'Xc'))=rADM(9);             % Composites
r(strcmp(compoundAbb,'Xsu'))=rADM(10);           % Sugar degraders
r(strcmp(compoundAbb,'Xlac'))=rADM(11);          % Lactate degraders


 growthRateLac  = rADM(11); % for analisis

end