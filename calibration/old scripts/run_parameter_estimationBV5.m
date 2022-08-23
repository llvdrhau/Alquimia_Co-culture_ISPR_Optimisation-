%Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% May 2020. Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.
% Run a parameter estimation on BV5

clear 
close all 

%load('testParEstimationBV5')
load('parEstimationBV5')

% manipulate spikes 2X
parameters.spikeFlowrate.compositionSpikes(2,11) = 1.28*1.34 ;     %Spike van Veillionella
parameters.spikeFlowrate.compositionSpikes(5,1) = 1067*1 ;     %last spike 
parameters.spikeFlowrate.compositionSpikes(4,1) = 1067*1 ;     %3th spike 
parameters.spikeFlowrate.compositionSpikes(3,1) = 1067*1 ;     %2nd spike 
% parameters.spikeFlowrate.compositionSpikes(1,1) = 0 ;     %1st spike 
%xInitial(2) = yExp1(1,1)+5;    

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
kdec_Xlac = 0.051;  % 13
KI_pro = 9;        % 14
KI_lac = 15;      %15
spike = 1.2*1.38;        % 16


idx = [13 16];



par2estimate = [km_su km_ye km_lac KS_su KS_ye KS_lac Ysu Yye Ylac KI_VFA lag_time K_YE kdec_Xlac KI_pro KI_lac spike]; % provide initial guess
par2names = {'km_su' 'km_ye' 'km_lac' 'KS_su' 'KS_ye' 'KS_lac','Ysu','Yye','Ylac','KI_VFA' 'lag_time','K_YE', 'kdec_Xlac', 'KI_pro', 'KI_lac' 'spike'};
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



