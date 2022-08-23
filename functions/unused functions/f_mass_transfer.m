function transfer=f_mass_transfer(states,parameters)
%% Gas-liquid transfer between the liquid and the headspace

% Alberte Regueira. University of Santiago de Compostela. Spain
% November 2019. Please contact alberte.regueira@gmail.com if you
% intend to use this code.

parValues = parameters.parValues;
parAbb = parameters.parAbb;
compoundAbb = parameters.compoundAbb;
Mw = parameters.Mw;
pH = parameters.pH.pHValue;
H = 10^-pH;
transfer = states*0;

%% Identification of compounds and parameters
iTic = strcmp(compoundAbb,'Stic');
iH2 = strcmp(compoundAbb,'Sh2');
iCH4 = strcmp(compoundAbb,'Sch4');

Stic = states(iTic);      %TIC
Sh2 = states(iH2);        %H2
Sch4 = states(iCH4);      %CH4
Sgas_h2 = states(strcmp(compoundAbb,'Sgas_h2'));    %H2 in gas headspace
Sgas_co2 = states(strcmp(compoundAbb,'Sgas_co2'));  %CO2 in gas headspace
Sgas_ch4 = states(strcmp(compoundAbb,'Sgas_ch4'));  %CH4 in gas headspace

Kac_tic = parValues(strcmp(parAbb,'Kac_tic'));
kLa = parValues(strcmp(parAbb,'kLa'));
%kLa= kLa+5000*(7-pH);
H_co2 = parValues(strcmp(parAbb,'H_co2'));
H_ch4 = parValues(strcmp(parAbb,'H_ch4'));
H_h2 = parValues(strcmp(parAbb,'H_h2'));

R = parameters.constants(1);    % Pa·L/mol·K
T = parameters.constants(2);    % K

Vgas = parameters.gasPhase(1);  % L

%% Stripping of CO2, H2 and CH4

% Partial pressure (Pa)
Ph2 = Sgas_h2*R*T/Mw(iH2);      %Pa
Pch4 = Sgas_ch4*R*T/Mw(iCH4);   %Pa
Pco2 = Sgas_co2*R*T/Mw(iTic);   %Pa

Sco2 = Stic*(H/(H+Kac_tic));

% Transfer
transfer(iH2) = kLa*(Mw(iH2)*H_h2*Ph2-Sh2);         % gCOD/h/L
transfer(iCH4) = kLa*(Mw(iCH4)*H_ch4*Pch4-Sch4);    % gCOD/h/L
transfer(iTic) = kLa*(Mw(iTic)*H_co2*Pco2-Sco2);    % gCOD/h/L

M=states(1);            % Total mass (kg)
rho=1;                  % kg/L    
V=M/rho;                % L

transfer(strcmp(compoundAbb,'Sgas_h2')) = -transfer(iH2)*V/Vgas;       % gCOD/h/L
transfer(strcmp(compoundAbb,'Sgas_ch4')) = -transfer(iCH4)*V/Vgas;     % gCOD/h/L
transfer(strcmp(compoundAbb,'Sgas_co2')) = -transfer(iTic)*V/Vgas;     % gCOD/h/L
