function [par2names,parameterReached] = f_run_veillionella_spike_estimation(kinetics)
%UNTITLED17 Summary of this function goes here
%   Detailed explanation goes here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter estimation using non-linear least squares
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alberte Regueira. University of Santiago de Compostela. Spain
% November 2019. Please contact alberte.regueira@gmail.com if you
% intend to use this code.
 
flagCase = 'SpikeVeillonella';
kineticFlag = 'CoCulture';
%% Reading experimental data
[data1,~]=xlsread('expData.xlsx','V3');
[data2,Names]=xlsread('expData.xlsx','V4');

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

% Reed controler (pos1) and/or spikes (pos2)
control = [0,1];

%% prep structure parameters
parameters=get_parameters(reactor,flagCase,kineticFlag);
parameters.flagCase = flagCase;
parameters.control = control;
parameters.kineticFlag = kineticFlag;
parameters.kineticInhibition = kinetics;

%% spike information exp3
spikeFlowrate = get_spike_data('V34');
parameters.spikeFlowrate = spikeFlowrate;

%% initial conditions 
xInitial(1,:)=get_initial_conditions(parameters,reactor,yExp1);
xInitial(2,:)=get_initial_conditions(parameters,reactor,yExp2);

%% Define estimation parameters

% Set parameters manually
% GLUCOSE UPTAKE SHOULD BE ZERO
parameters.parValues(strcmp(parameters.parAbb,'km_su')) = 0;
% lactate uptake rate  = estimated from batch
parameters.parValues(strcmp(parameters.parAbb,'Ylac')) =  0.0371;    %0.0389;
parameters.parValues(strcmp(parameters.parAbb,'km_lac')) = 10.5179;     %7.635;
parameters.parValues(strcmp(parameters.parAbb,'KS_lac')) = 0.1; %0.2735; %0.005; %0.2735;   %0.5;
parameters.parValues(strcmp(parameters.parAbb,'kdec_Xlac')) = 0.0729; %0.05;

%yeast uptake rate Bacillus => their is none 
parameters.parValues(strcmp(parameters.parAbb,'Yyeb')) = 0;
parameters.parValues(strcmp(parameters.parAbb,'KS_yeb')) = 0.5; %0.5

%yeast uptake rate veillionella 
parameters.parValues(strcmp(parameters.parAbb,'Yyev')) = 0.08; %0.02;
parameters.parValues(strcmp(parameters.parAbb,'KS_yev')) = 0.5;  %15; %0.5
parameters.parValues(strcmp(parameters.parAbb,'K_YE')) = 1;
%parameters.parValues(strcmp(parameters.parAbb,'km_ye')) = 1;  %not
%necesary is the same as km_lac

% NO Xglu DECAY
parameters.parValues(strcmp(parameters.parAbb,'kdec_Xsu')) = 0;
% VFA inhibition 
parameters.parValues(strcmp(parameters.parAbb,'KI_VFA')) = 1000;
parameters.parValues(strcmp(parameters.parAbb,'lag_time')) = 0;

% VFA inhibition 
parameters.parValues(strcmp(parameters.parAbb,'KI_pro')) = 1000;
parameters.parValues(strcmp(parameters.parAbb,'KI_lac')) = 1000;


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
KI_pro = 8;        % 14
KI_lac = 10;      %15

switch kinetics
    case 'totalVFA'
        idx = [3 9 10 13];
        %idx = [10];
    case 'seperateProLac'
        idx = [3 9 14 15 13];
        %idx = [14 15];
    case 'proOnly'
        idx = [3 9 14 13];
        %idx = [14];
    case 'lacOnly'
        idx = [3 9 15 13];
        %idx = [15];
end
    

par2estimate = [km_su km_ye km_lac KS_su KS_ye KS_lac Ysu Yye Ylac KI_VFA lag_time K_YE kdec_Xlac KI_pro KI_lac]; % provide initial guess
par2names = {'km_su' 'km_ye' 'km_lac' 'KS_su' 'KS_ye' 'KS_lac','Ysu','Yye','Ylac','KI_VFA' 'lag_time','K_YE', 'kdec_Xlac', 'KI_pro', 'KI_lac'};
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
[parameterReached,normResidualresidual, residual, ~, ~,~, jacobian] = lsqnonlin(@f_residual_V_spike,par2estimate,lb,ub,options,tExp,yExp,xInitial,parameters,par2names,caso);
toc
PR(1,:) = par2names;
PR(2,:) = num2cell(parameterReached);
disp(['Estimated parameters from:',kinetics])
disp(PR)
disp('')
%% Solution of the model with the first parameters estimation
caso = 'normal';
% check these inhibitions 
% parameterReached(3) = 8.026;
% parameterReached(4) = 15.43;

% maybe use the residual to get simulataion data 
[t1,y1] = f_solution(xInitial(1,:),parameters,tExp,parameterReached,par2names,caso);
[t2,y2] = f_solution(xInitial(2,:),parameters,tExp,parameterReached,par2names,caso);

%% RMSE
compounds = {'Lac', 'Pro', 'Ace', 'BM'};

pos_lac = find((strcmp(parameters.compoundAbb,'Slac')));
pos_pro = find((strcmp(parameters.compoundAbb,'Spro')));
pos_ac = find((strcmp(parameters.compoundAbb,'Sac')));
pos_bm = find((strcmp(parameters.compoundAbb,'Xlac')));
%pos_ye = (strcmp(parameters.compoundAbb,'Sye'));
col = [pos_lac,pos_pro,pos_ac,pos_bm];
ysim = y1(:,col); 
[RMSE, r2] = f_RMSE(tExp,yExp1,ysim,t1,compounds); 

%% plot graphs
y = {y1 y2};
t = {t1 t2};
out = [t1 y1];

yExpGraf = {yExp1 yExp2};
tExpGraf = {tExp tExp};

    titleGraf = ['Inhibition:',kinetics];
    plot_V_batch(tExpGraf,yExpGraf,t,y,parameters,titleGraf, r2)

 
end

