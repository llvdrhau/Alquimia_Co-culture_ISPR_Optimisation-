function [par2names,parameterReached] = f_run_parameter_estimation_BV(loadFile,parNames)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


%Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% May 2020. Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.
% Runs a parameter estimation on the cocultures 

%load('testParEstimationBV5')
load(loadFile)

% Initial guessing values
km_su=2;           % 1  Glucose -> Lac
km_ye=1;           % 2  YE -> ace + pro 
km_lac=10;          % 3  3Lac -> 2Pro + Ac
KS_su=0.1;         % 4  
KS_ye=0.1;         % 5
KS_lac=0.27;          % 6
Ysu = 0.2;         % 7
Yye = 0.6;         % 8
Ylac = 0.022;      % 9
KI_VFA = 5;       % 10
lag_time = 2;      % 11
K_YE = 0.68;        % 12
kdec_Xlac = 0.03;  % 13
KI_pro = 6;        % 14
KI_lac = 10;      %15


allNames = {'km_su' 'km_ye' 'km_lac' 'KS_su' 'KS_ye' 'KS_lac','Ysu','Yye','Ylac','KI_VFA' 'lag_time','K_YE', 'kdec_Xlac', 'KI_pro', 'KI_lac'};
idx = zeros(1,length(parNames));
for i = 1:length(parNames)
    idx(i) = find(strcmp(allNames,parNames(i)),2);
end
%idx = [14 15];



par2estimate = [km_su km_ye km_lac KS_su KS_ye KS_lac Ysu Yye Ylac KI_VFA lag_time K_YE kdec_Xlac KI_pro KI_lac]; % provide initial guess
par2names = {'km_su' 'km_ye' 'km_lac' 'KS_su' 'KS_ye' 'KS_lac','Ysu','Yye','Ylac','KI_VFA' 'lag_time','K_YE', 'kdec_Xlac', 'KI_pro', 'KI_lac'};
par2estimate = par2estimate(idx);
par2names = par2names(idx);

options = optimset('display', 'iter','tolfun',1.0e-06, 'tolx',1.0e-5, 'maxfunevals', 2000);
lb=par2estimate*0;
ub=par2estimate*10;


%% Parameter first estimation
yExp = yExp1(:,1:end-1); % don't include biomass for estimation because of the mix?

tic
caso = 'normal';
%[parameterReached,normResidualresidual, residual, ~, ~,~, jacobian] = lsqnonlin(@f_residual_BV5,par2estimate,lb,ub,options,tExp,yExp,xInitial,parameters,par2names,caso);
[parameterReached] = fmincon(@f_residual_BV5,par2estimate,[],[],[],[],lb,ub,[],[],tExp,yExp,xInitial,parameters,par2names,caso);
toc

%%
PRINT(1,:) = par2names;
PRINT(2,:) = num2cell(parameterReached);
disp('Os parámetros estimados son:')
disp(PRINT)
disp('')

%% plot results 

[~,t,y] = f_residual_BV5(parameterReached,tExp,yExp,xInitial,parameters,par2names,caso);
title = 'BV5 Co-Culture inhibitions estimated';
plot_CoCulture(tExp,yExp1,t,y,parameters,title)



end

