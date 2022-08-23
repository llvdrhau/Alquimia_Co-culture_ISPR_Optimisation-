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
kineticFlag = 'V_inhibition4';
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

parameters=get_parameters(reactor,flagCase);
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
parameters.parValues(strcmp(parameters.parAbb,'km_lac')) = 1;
parameters.parValues(strcmp(parameters.parAbb,'KS_lac')) = 0.5;
parameters.parValues(strcmp(parameters.parAbb,'kdec_Xlac')) = 0;%0.05;
%yeast uptake rate 
parameters.parValues(strcmp(parameters.parAbb,'km_ye')) = 0;
parameters.parValues(strcmp(parameters.parAbb,'KS_ye')) = 0.1;
% NO Xglu DECAY
parameters.parValues(strcmp(parameters.parAbb,'kdec_Xsu')) = 0;

% VFA inhibition 
parameters.parValues(strcmp(parameters.parAbb,'KI_pro')) = 8;
parameters.parValues(strcmp(parameters.parAbb,'KI_lac')) = 10;

parameters.parValues(strcmp(parameters.parAbb,'lag_time')) = 0;
% Initial guessing values

% Initial guessing values
km_su=2;           % 1  Glucose -> Lac
km_ye=3;           % 2  YE -> ace + pro 
km_lac=1;          % 3  3Lac -> 2Pro + Ac
KS_su=0.1;           % 4  
KS_ye=0.1;           % 5
KS_lac=1;          % 6
Ysu = 0.2;         % 7
Yye = 0.1;         % 8
Ylac = 0.1;        % 9
KI_VFA = 0;        %10

lag_time = 2;      % 11
K_YE = 0.68;        % 12
kdec_Xlac = 0.03;  % 13
KI_pro = 8;        % 14
KI_lac = 10;      %15


idx = [3 9 2 8 13 15]; % Parameters to be determined from the list above

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
if graficas
    title = 'Substrate inhibition Lactate';
    plot_V_batch(tExpGraf,yExp,t,y,parameters,title)
end
%plot_reactor_results(t,y,parameters)



%% Extended analysis
for i=1
% 
% %% Sensitivity of the model outcome with the parameters
% [dydp,dyNDdp] = f_sensitivity(y,idx,parameterReached,xInitial,parameters,tExp,par2names);
% tijd = 1:size(dydp,1);
% hold off
% hold on
% for i = 1:size(dydp,2)
% plot(tijd, dydp(:,i))
% legend(par2names)
% end
% hold off
% 
% %% Identifiability analysis
% [col,dmsqr,sNorm,Combi] = f_ident_analysis(dydp,dyNDdp,idx,par2names);
% disp('A sensibilidade (dmsqr) dos parámetros é:')
% disp(dmsqr)
% disp('')
% 
% %% Bootstrap method
% [parameterReported, parameterMatrix, parameterReported_dev,CImatrix] = f_bootstrap(yExp,par2estimate,lb,ub,options,tExp,xInitial,parameters,par2names,parameterReached);
% 
% if graficas
%     % Histogram of the parameters space estimated with the Boostrap method
%     for i=1:size(parameterMatrix,2)
%         figure
%         hist(parameterMatrix(:,i))
%         title(par2names{i})
%     end
%     
%     % Correlation coefficient between parameters
%     for i=1:size(parameterMatrix,2)
%         coefs = corrcoef(parameterMatrix(:,2),parameterMatrix(:,i));
%         m_coefs(1,i) = coefs(1,2);
%     end
% end
% %% Monte Carlo method
% 
% [Matrix1, Matrix2] = f_montecarlo(par2estimate,CImatrix,xInitial,parameters,tExp,par2names,t);
% 
% if graficas
%     f_plot_montecarlo(Matrix1,Matrix2,tExp,yExp)
% end
% 
% %% RMSE determination
% [ RMSE, RMSE_ind] = f_RMSE(tExp,yExp,Matrix1,Matrix2);
% 
% 
% %% Organize results and saving
% resultados_finais(1,:) = par2names;
% resultados_finais(2,:) = num2cell(parameterReported);
% resultados_finais(3:4,:) = CImatrix(2:3,:);
% resultados_finais{1,length(par2names)+1} = 'RMSE';
% resultados_finais(2:3,length(par2names)+1) = num2cell(RMSE);
% 
% FIN.RMSE = RMSE;
% FIN.RMSE_ind = RMSE_ind;
% FIN.Matrix = Matrix;
% FIN.parameterReported = parameterReported;
% FIN.parameterMatrix = parameterMatrix;
% FIN.parameterReported_dev = parameterReported_dev;
% FIN.CImatrix = CImatrix;
% FIN.dmsqr = dmsqr;
% FIN.parameterReached = parameterReached;
% FIN.resultados_finais = resultados_finais;
% 
% save resultados.mat
% % print(gcf, 'X_pH7.jpg', '-djpeg', '-r300' );
% 
end