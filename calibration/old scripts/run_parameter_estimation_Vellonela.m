%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter estimation using non-linear least squares
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alberte Regueira. University of Santiago de Compostela. Spain
% November 2019. Please contact alberte.regueira@gmail.com if you
% intend to use this code.

clear;  
% close all;  
graficas = 1;
flagCase = 'batchVellonela';
kineticFlag = 'V_inhibition1';
%% Reading experimental data
[data1,Names]=xlsread('VellonelaBatch1.xlsx','MATLAB');
[data2,Names]=xlsread('VellonelaBatch2.xlsx','MATLAB');

tExp=data1(:,strcmp(Names,'time'));
yExp1(:,1)=data1(:,strcmp(Names,'Lac'));
yExp1(:,2)=data1(:,strcmp(Names,'Pro'));
yExp1(:,3)=data1(:,strcmp(Names,'Ace'));
yExp1(:,4)=data1(:,strcmp(Names,'BM'));

% tExp2=data2(:,strcmp(Names,'t'));
yExp2(:,1)=data2(:,strcmp(Names,'Lac'));
yExp2(:,2)=data2(:,strcmp(Names,'Pro'));
yExp2(:,3)=data2(:,strcmp(Names,'Ace'));
yExp2(:,4)=data2(:,strcmp(Names,'BM'));

yExp = [yExp1 yExp2];

%% Define model parameters

reactor.V=1; %L
reactor.pH(1)=1;             % Define whether pH is kept constant or not
reactor.pH(2)=7;             % pH value if it's controlled

% controlers and/or spikes 
control = [0,0];

parameters=get_parameters(reactor,flagCase,kineticFlag);
parameters.flagCase = flagCase;
parameters.control = control;
parameters.kineticFlag = kineticFlag;

xInitial(1,:)=get_initial_conditions(parameters,reactor,yExp1);
xInitial(2,:)=get_initial_conditions(parameters,reactor,yExp2);

%% Define estimation parameters

% Set parameters manually
% GLUCOSE UPTAKE SHOULD BE ZERO
parameters.parValues(strcmp(parameters.parAbb,'km_su')) = 0;
% lactate uptake rate 
parameters.parValues(strcmp(parameters.parAbb,'km_lac')) = 6;
parameters.parValues(strcmp(parameters.parAbb,'KS_lac')) = 0.2735; %0.5;
parameters.parValues(strcmp(parameters.parAbb,'kdec_Xlac')) = 0;%0.05;
%yeast uptake rate 
parameters.parValues(strcmp(parameters.parAbb,'Yye')) = 0.02;
parameters.parValues(strcmp(parameters.parAbb,'km_ye')) = 1;
parameters.parValues(strcmp(parameters.parAbb,'KS_ye')) = 15; %0.5
parameters.parValues(strcmp(parameters.parAbb,'K_YE')) = 1;
% NO Xglu DECAY
parameters.parValues(strcmp(parameters.parAbb,'kdec_Xsu')) = 0;
% VFA inhibition 
parameters.parValues(strcmp(parameters.parAbb,'KI_VFA')) = 1000;
parameters.parValues(strcmp(parameters.parAbb,'lag_time')) = 0;
% Initial guessing values
km_su=2;           % 1  Glucose -> Lac
km_ye=1;           % 2  YE -> ace + pro 
km_lac=3;          % 3  3Lac -> 2Pro + Ac
KS_su=0.1;         % 4  
KS_ye=0.1;         % 5
KS_lac=1;          % 6
Ysu = 0.2;         % 7
Yye = 0.6;         % 8
Ylac = 0.1;        % 9
KI_VFA = 26;       % 10
lag_time = 2;      % 11
K_YE = 0.68;       % 12
kdec_Xlac = 0.06;  % 13

%idx = [2 5 8 9 10 11 12 13]; % Parameters to be determined from the list above
%idx = [2 3 5 8 9 10 11 12];
idx = [3 6 9 13];


par2estimate = [km_su km_ye km_lac KS_su KS_ye KS_lac Ysu Yye Ylac KI_VFA lag_time K_YE kdec_Xlac]; % provide initial guess
par2names = {'km_su' 'km_ye' 'km_lac' 'KS_su' 'KS_ye' 'KS_lac','Ysu','Yye','Ylac','KI_VFA' 'lag_time','K_YE', 'kdec_Xlac'};
par2estimate = par2estimate(idx);
par2names = par2names(idx);

options = optimset('display', 'iter','tolfun',1.0e-06, 'tolx',1.0e-5, 'maxfunevals', 2000);
lb=par2estimate*0;
% lb(strcmp(par2names,'lag_time')) = 3;
ub=par2estimate*10;
% ub(end)=1;

%% Parameter first estimation
tic
caso = 'normal';
[parameterReached,normResidualresidual, residual, ~, ~,~, jacobian] = lsqnonlin(@f_residual_V_Batch,par2estimate,lb,ub,options,tExp,yExp,xInitial,parameters,par2names,caso);
toc
PR(1,:) = par2names;
PR(2,:) = num2cell(parameterReached);
disp('Os parámetros estimados son:')
disp(PR)
disp('')
%% Solution of the model with the first parameters estimation
caso = 'normal';

% maybe use the residual to get simulataion data 
[t1,y1] = f_solution(xInitial(1,:),parameters,tExp,parameterReached,par2names,caso);
[t2,y2] = f_solution(xInitial(2,:),parameters,tExp,parameterReached,par2names,caso);

y = {y1 y2};
t = {t1 t2};
yExp = {yExp1 yExp2};
tExpGraf = {tExp tExp};

%plot_reactor_results(t,y,parameters)

compounds = {'Lac', 'Pro', 'Ace', 'BM'};

pos_lac = find((strcmp(parameters.compoundAbb,'Slac')));
pos_pro = find((strcmp(parameters.compoundAbb,'Spro')));
pos_ac = find((strcmp(parameters.compoundAbb,'Sac')));
pos_bm = find((strcmp(parameters.compoundAbb,'Xlac')));
%pos_ye = (strcmp(parameters.compoundAbb,'Sye'));
col = [pos_lac,pos_pro,pos_ac,pos_bm];
ysim = y1(:,col); 
[RMSE,r2] = f_RMSE(tExp,yExp1,ysim,t1,compounds); 

%% plot graphs 
if graficas
    title = 'Batch Veillionella';
    plot_V_batch(tExpGraf,yExp,t,y,parameters,title,r2)
end



