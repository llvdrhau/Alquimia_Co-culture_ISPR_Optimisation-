%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter estimation using non-linear least squares
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alberte Regueira. University of Santiago de Compostela. Spain
% November 2019. Please contact alberte.regueira@gmail.com if you
% intend to use this code.

clear;  
% close all; 

% options 
graficas = 1;
flagCase = 'batchBacilus';
kineticFlag = 'GiveAName';
nr = 3;
%% Reading experimental data
[data1,Names]=xlsread('Biochem 5 Reactor 1 BMH.xlsx','MATLAB');
[data2]=xlsread('Biochem 5 Reactor 2 BMH.xlsx','MATLAB');
[data3]=xlsread('Biochem 4 R3 BMH.xlsx','MATLAB','A2:D11');

% time cell 
tExp1=data1(:,strcmp(Names,'t'));
tExp2=data2(:,strcmp(Names,'t'));
tExp3=data3(:,strcmp(Names,'t'));

tExp= {tExp1,tExp2,tExp3};

% data cell
yExp1(:,1)=data1(:,strcmp(Names,'Glucose'));
yExp1(:,2)=data1(:,strcmp(Names,'Lactate'));
yExp1(:,3)=data1(:,strcmp(Names,'BM'));

% tExp2=data2(:,strcmp(Names,'t'));
yExp2(:,1)=data2(:,strcmp(Names,'Glucose'));
yExp2(:,2)=data2(:,strcmp(Names,'Lactate'));
yExp2(:,3)=data2(:,strcmp(Names,'BM'));

yExp3(:,1)=data3(:,strcmp(Names,'Glucose'));
yExp3(:,2)=data3(:,strcmp(Names,'Lactate'));
yExp3(:,3)=data3(:,strcmp(Names,'BM'));

yExp = {yExp1, yExp2, yExp3};

%% Define model parameters

reactor.V=1; %L
reactor.pH(1)=1;             % Define whether pH is kept constant or not
reactor.pH(2)=7;             % pH value if it's controlled

% parameters for each experiment 
parameters1 = get_parameters(reactor,flagCase);
parameters2 = get_parameters(reactor,flagCase);
parameters3 = get_parameters(reactor,flagCase);

%give the kineticFlag to the parameters
parameters1.kineticFlag = kineticFlag;
parameters2.kineticFlag = kineticFlag;
parameters3.kineticFlag = kineticFlag;

% are spike aded or feeding? first element feed second Spkies
control1 = [0, 0]; %Spikes in experiment 3  
control2 = [0, 0]; %Spikes in experiment 3  
control3 = [0, 1]; %Spikes in experiment 3  

parameters1.control = control1;
parameters2.control = control2;
parameters3.control = control3;

% spike information exp3

spikeFlowrate.tSpikes = [14, 20]; %h
spikeFlowrate.volumeSpikes = [0.035, 0.035]; % in L 
spikeFlowrate.compositionSpikes = zeros(14,2);
spikeFlowrate.compositionSpikes(2,1) = 400*1.07 ; % in gCOD/L, only glucose spikes   
spikeFlowrate.compositionSpikes(2,2) = 400*1.07 ; % in gCOD/L, only glucose spikes % composition at each time
parameters3.spikeFlowrate = spikeFlowrate;

%initial conditions 
xInitial(1,:)=get_initial_conditions(parameters1,reactor,yExp1);
xInitial(2,:)=get_initial_conditions(parameters2,reactor,yExp2);
xInitial(3,:)=get_initial_conditions(parameters3,reactor,yExp3);


% give to over all structure parameters 
parameters.parameters1 = parameters1;
parameters.parameters2 = parameters2;
parameters.parameters3 = parameters3;
%% Define estimation parameters

% Set parameters manually
names = fieldnames(parameters);
for i = 1:length(names)
    par = parameters.(names{i});
    
    par.parValues(strcmp(par.parAbb,'km_su')) = 0;
    par.parValues(strcmp(par.parAbb,'km_ye')) = 0.1;
    par.parValues(strcmp(par.parAbb,'km_lac')) = 0;
    par.parValues(strcmp(par.parAbb,'KS_su')) = 0.1;
    par.parValues(strcmp(par.parAbb,'KS_ye')) = 0.1;
    par.parValues(strcmp(par.parAbb,'KS_lac')) = 0.5;
    par.parValues(strcmp(par.parAbb,'kdec_Xlac')) = 0;
    par.parValues(strcmp(par.parAbb,'kdec_Xsu')) = 0;
    par.parValues(strcmp(par.parAbb,'KI_VFA')) = 1e10;
    
    parametersAll.(names{i})= par;
end


% Initial guessing values
km_su=4.276;           % 1  Glucose -> Lac
km_ye=0.1;           % 2  YE -> Lac
km_lac=1;          % 3  3Lac -> 2Pro + Ac
KS_su=0.1;           % 4  
KS_ye=0.1;           % 5
KS_lac=1;          % 6
Ysu = 0.0714;         % 7
Yye = 0.1;         % 8
Ylac = 0.1;        % 9
KI_VFA = 5;        % 10
lag_time = 2.2013;      % 11
K_YE = 0.4184;        % 12    Is the % of yeast eventualy consumed!

idx = [1 7 11 12]; % Parameters to be determined from the list above
%idx = [1 7];
par2estimate = [km_su km_ye km_lac KS_su KS_ye KS_lac Ysu Yye Ylac KI_VFA lag_time K_YE]; % provide initial guess
par2names = {'km_su' 'km_ye' 'km_lac' 'KS_su' 'KS_ye' 'KS_lac','Ysu','Yye','Ylac','KI_VFA' 'lag_time','K_YE'};
par2estimate = par2estimate(idx);
par2names = par2names(idx);


%% Parameter first estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nr = 2 ;
if nr == 2
%just the two batch experiments comment for three experiments 
tExp = {tExp1,tExp2};
yExp = {yExp1,yExp2};

parameters2batch.parameters1 = parametersAll.parameters1;
parameters2batch.parameters2 = parametersAll.parameters2;
parametersAll = parameters2batch;
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimation 

% lower and upper bound 
lb=par2estimate*0;    % lb(strcmp(par2names,'lag_time')) = 3;
ub=par2estimate*10;    % ub(end)=1;
caso = 'normal';

%lsqnonlin
%options 
options = optimset('display', 'iter','tolfun',1.0e-06, 'tolx',1.0e-5, 'maxfunevals', 2000);
tic
[parameterReached,normResidualresidual, residual, ~, ~,~, jacobian] = lsqnonlin(@f_residual,par2estimate,lb,ub,options,tExp,yExp,xInitial,parametersAll,par2names,caso,'lsq');
toc


% Estimation fconmin 
%optionsFcon=optimoptions('lsqnonlin','Display','iter-detailed');
% tic
% [parameterReached,normResidualresidual, residual, ~, ~,~, jacobian] = fmincon(@f_residual,par2estimate,[],[],[],[],lb,ub,[],[],tExp,yExp,xInitial,parametersAll,par2names,caso,'fconmin');
% toc

%% display on screen 
PR(1,:) = par2names;
PR(2,:) = num2cell(parameterReached);
disp('Os parámetros estimados son:')
disp(PR)
disp('')
%% Solution of the model with the first parameters estimation
caso = 'normal';

[t1,y1] = f_solution(xInitial(1,:),parameters.parameters1,tExp1,parameterReached,par2names,caso);
[t2,y2] = f_solution(xInitial(2,:),parameters.parameters2,tExp2,parameterReached,par2names,caso);
[t3,y3] = f_solution(xInitial(3,:),parameters.parameters3,tExp3,parameterReached,par2names,caso);

y = {y1,y2,y3}; 
t = {t1,t2,t3};

if nr == 2
    y = {y1,y2};
    t = {t1,t2};
end

if graficas
    f_plot_4grafs(tExp,yExp,t,y,parameters.parameters1)
end

exp3 =[t3,y3];
%plot_reactor_results(t,y,parameters)


%% %TO BE CONTINUED FROM HERE ON (uncomment)
goOn = 0;
if goOn == 1 
for i=1 

%% Sensitivity of the model outcome with the parameters
[dydp,dyNDdp] = f_sensitivity(y,idx,parameterReached,xInitial,parameters,tExp,par2names);
tijd = 1:size(dydp,1);
hold off
hold on
for i = 1:size(dydp,2)
plot(tijd, dydp(:,i))
legend(par2names)
end
hold off

%% Identifiability analysis
[col,dmsqr,sNorm,Combi] = f_ident_analysis(dydp,dyNDdp,idx,par2names);
disp('A sensibilidade (dmsqr) dos parámetros é:')
disp(dmsqr)
disp('')

%% Bootstrap method
[parameterReported, parameterMatrix, parameterReported_dev,CImatrix] = f_bootstrap(yExp,par2estimate,lb,ub,options,tExp,xInitial,parameters,par2names,parameterReached);

if graficas
    % Histogram of the parameters space estimated with the Boostrap method
    for i=1:size(parameterMatrix,2)
        figure
        hist(parameterMatrix(:,i))
        title(par2names{i})
    end
    
    % Correlation coefficient between parameters
    for i=1:size(parameterMatrix,2)
        coefs = corrcoef(parameterMatrix(:,2),parameterMatrix(:,i));
        m_coefs(1,i) = coefs(1,2);
    end
end
%% Monte Carlo method

[Matrix1, Matrix2] = f_montecarlo(par2estimate,CImatrix,xInitial,parameters,tExp,par2names,t);

if graficas
    f_plot_montecarlo(Matrix1,Matrix2,tExp,yExp)
end

%% RMSE determination
[ RMSE, RMSE_ind] = f_RMSE(tExp,yExp,Matrix1,Matrix2);


%% Organize results and saving
resultados_finais(1,:) = par2names;
resultados_finais(2,:) = num2cell(parameterReported);
resultados_finais(3:4,:) = CImatrix(2:3,:);
resultados_finais{1,length(par2names)+1} = 'RMSE';
resultados_finais(2:3,length(par2names)+1) = num2cell(RMSE);

FIN.RMSE = RMSE;
FIN.RMSE_ind = RMSE_ind;
FIN.Matrix = Matrix;
FIN.parameterReported = parameterReported;
FIN.parameterMatrix = parameterMatrix;
FIN.parameterReported_dev = parameterReported_dev;
FIN.CImatrix = CImatrix;
FIN.dmsqr = dmsqr;
FIN.parameterReached = parameterReached;
FIN.resultados_finais = resultados_finais;

save resultados.mat
% print(gcf, 'X_pH7.jpg', '-djpeg', '-r300' );
end
end