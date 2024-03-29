%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter estimation using non-linear least squares Bacillus data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lucas Van der Hauwaert & Alberte Regueira.  
% University of Santiago de Compostela. Spain
% October 2021. Please contact alberte.regueira@gmail.com or 
% lucas.vanderhauwaert@usc.es if you
% intend to use this code.

clear;  
close all;  

% switchs to run the bootstrap or motecarlo 
% if 0 saved data is used to plot the graphs and distributions
graficas = 1;
runNewBoot = 0;
newMonteCarlo = 0; 

flagCase = 'batchBacilus';
kineticFlag = 'anything';

%% Reading experimental data
[data1,~]=xlsread('Biochem 5 Reactor 1 BMH.xlsx','MATLAB');
[data2,Names]=xlsread('Biochem 5 Reactor 2 BMH.xlsx','MATLAB');

tExp=data1(:,strcmp(Names,'t'));
yExp1(:,1)=data1(:,strcmp(Names,'Glucose'));
yExp1(:,2)=data1(:,strcmp(Names,'Lactate'));
yExp1(:,3)=data1(:,strcmp(Names,'BM'));

% tExp2=data2(:,strcmp(Names,'t'));
yExp2(:,1)=data2(:,strcmp(Names,'Glucose'));
yExp2(:,2)=data2(:,strcmp(Names,'Lactate'));
yExp2(:,3)=data2(:,strcmp(Names,'BM'));

yExp = [yExp1 yExp2];

%% Define model parameters

reactor.V=1; %L
reactor.pH(1)=1;             % Define whether pH is kept constant or not
reactor.pH(2)=7;             % pH value if it's controlled

parameters=get_parameters(reactor,flagCase,kineticFlag); % initial parameters
parameters.control = [0 0];
parameters.kineticFlag = kineticFlag;


xInitial(1,:)=get_initial_conditions(parameters,reactor,yExp1);
xInitial(2,:)=get_initial_conditions(parameters,reactor,yExp2);

%% Define estimation parameters

% Set parameters manually
parameters.parValues(strcmp(parameters.parAbb,'km_su')) = 3;
parameters.parValues(strcmp(parameters.parAbb,'km_ye')) = 0;
parameters.parValues(strcmp(parameters.parAbb,'km_lac')) = 0;
parameters.parValues(strcmp(parameters.parAbb,'KS_su')) = 0.1;
parameters.parValues(strcmp(parameters.parAbb,'KS_ye')) = 0.1;
parameters.parValues(strcmp(parameters.parAbb,'KS_lac')) = 0.5;
parameters.parValues(strcmp(parameters.parAbb,'kdec_Xlac')) = 0;
parameters.parValues(strcmp(parameters.parAbb,'kdec_Xsu')) = 0;
parameters.parValues(strcmp(parameters.parAbb,'KI_VFA')) = 1e10;
parameters.parValues(strcmp(parameters.parAbb,'lag_time')) = 0;
parameters.parValues(strcmp(parameters.parAbb,'K_YE')) = 0;

% Initial guessing values
km_su=3;           % 1  Glucose -> Lac
km_ye=3;           % 2  YE -> Lac
km_lac=1;          % 3  3Lac -> 2Pro + Ac
KS_su=0.1;           % 4  
KS_ye=0.1;           % 5
KS_lac=1;          % 6
Ysu = 0.2;         % 7
Yye = 0.1;         % 8
Ylac = 0.1;        % 9
KI_VFA = 5;        % 10
lag_time = 10;      % 11
K_YE = 0.8;        % 12
K_YE = 0.4;

% Parameters to be determined from the list above
idx = [1 7 12];
%idx = [1 7];
%idx = [1 7 11 12];

par2estimate = [km_su km_ye km_lac KS_su KS_ye KS_lac Ysu Yye Ylac KI_VFA lag_time K_YE]; % provide initial guess
par2names = {'km_su' 'km_ye' 'km_lac' 'KS_su' 'KS_ye' 'KS_lac','Ysu','Yye','Ylac','KI_VFA' 'lag_time','K_YE'};
par2estimate = par2estimate(idx);
par2names = par2names(idx);

% options for lsqnonlin
options = optimset('display', 'iter','tolfun',1.0e-06, 'tolx',1.0e-5, 'maxfunevals', 2000);
% upper and lower bounds 
lb=par2estimate*0;
ub=par2estimate*10;

%% Parameter first estimation
tic
caso = 'normal';
[parameterReached,normResidualresidual, residual, ~, ~,~, jacobian] = lsqnonlin(@f_residual_bacillus,par2estimate,lb,ub,options,tExp,yExp,xInitial,parameters,par2names,caso);
toc
PR(1,:) = par2names;
PR(2,:) = num2cell(parameterReached);
disp('Os parámetros estimados son:')
disp(PR)
disp('')
%% Solution of the model with the first parameters estimation
caso = 'normal';
%[t,~,y1,y2] = f_solution_bacillus(xInitial,parameters,tExp,parameterReached,par2names,caso);
[t,~,y1,y2] = f_solution_bacillus(xInitial,parameters,tExp,parameterReached,par2names,caso);

y = {y1 y2};
tGraf = {t t};
yExpGraf = {yExp1 yExp2};
tExpGraf = {tExp tExp};
if graficas
    f_plot_bacillus(tExpGraf,yExpGraf,tGraf,y,parameters)
end


%% RMSE

% select compounds 
compounds = {'Glu', 'Lac', 'BM'}; % so names of compounds are displayed
% positions (columns) of data in y 

pos_glu = find((strcmp(parameters.compoundAbb,'Sglu')));
pos_lac = find((strcmp(parameters.compoundAbb,'Slac')));
pos_pro = find((strcmp(parameters.compoundAbb,'Spro')));
pos_ac = find((strcmp(parameters.compoundAbb,'Sac')));
pos_bm = find((strcmp(parameters.compoundAbb,'Xsu')));

col = [pos_glu,pos_lac,pos_bm]; % define which elements position 
ysim = y1(:,col);  % values of simulation from choosen compounds

% seperate datasets 

% %dataset 1
% [RMSE, r2] = f_RMSE(tExp,yExp1,ysim,t,compounds); % calculates RSME and R^2
% %dataset 2
% ysim = y2(:,col); 
% [RMSE2, r3] = f_RMSE(tExp,yExp2,ysim,t,compounds); 

% collective datasets
ysim = {y1(:,col) y2(:,col)};
tAll = {t t};
yExpAll = {yExp1 yExp2};
tExpAll = {tExp tExp};

[RMSE2, rSquared] = f_RMSE(tExpAll,yExpAll,ysim,tAll,compounds); 
%% Bootstrap method

[parameterReported, parameterMatrix, parameterReported_dev,CImatrix] = f_bootstrap(yExp,par2estimate,lb,ub,options,tExp,xInitial,parameters,par2names,parameterReached,runNewBoot);
path =cd();
saveName = '\saved data\bootstrapResultsBacillus2';
savePath = append(path,saveName);

save(savePath,...
   'parameterReported','parameterReported', 'parameterMatrix', 'parameterReported_dev','CImatrix')

% save('C:\BioChemCode_official\calibration\saved data\bootstrapResultsBacillus2',...
%     'parameterReported','parameterReported', 'parameterMatrix', 'parameterReported_dev','CImatrix')
%% Histogram of the parameters space estimated with the Boostrap method
%load('bootstrapResultsBacillus2')

nPar = size(parameterMatrix,2);
parNamesLx = {'km_{su}' 'km_{ye}' 'km_{lac}' 'KS_{su}' 'KS_{ye}' 'KS_{lac}','Y_{su}','Y_{ye}','Y_{lac}','KI_{VFA}' 'lag_time','K_{YE}', 'kdec_{Xlac}', 'KI_{pro}', 'KI_{lac}'};
parUnits = {'(h^{-1})','(gCOD_{X}/gCOD_{glu})', '(-)'};
parLatex = parNamesLx(idx);
if graficas
    for i=1:nPar
        %figure
        subplot(1,3,i)
        set(gcf,'color','w'); %background color
        histogram(parameterMatrix(:,i))
        txtX = [parLatex{i} parUnits{i}];
        xlabel(txtX)
        ylabel('Frequency')
    end
end

%% Monte Carlo method
meanPar = parameterReported; 
sdPar = cell2mat(parameterReported_dev(2,:));

if newMonteCarlo
    tic
    [Matrix1, Matrix2] = f_montecarlo(par2estimate,parameterMatrix,xInitial,parameters,tExp,par2names,meanPar,sdPar);
    toc

    saveName = '\saved data\MonteCarloBacillus';
    savePath = append(path,saveName);
    save(savePath,...
        'Matrix1', 'Matrix2','tExp','yExp')

end

  
%% montecarlo plots
if graficas
    load('MonteCarloBacillus')
    f_plot_montecarlo(Matrix1,Matrix2,tExp,yExp)
end

%% rewrite parameters to a save file 

parBootStrapBacillus = parameters.parValues;
for i = 1:length(par2names)
    pos = strcmp(parameters.parAbb,par2names(i));
    parBootStrapBacillus(pos) = parameterReported(i);
end
par2namesB = par2names;
parameterReportedB = parameterReported; 

saveName = '\saved data\parBootStrapBacillus';
savePath = append(path,saveName);
save(savePath, ...
    'par2namesB' ,'parameterReportedB', 'parBootStrapBacillus')

%% check RMSE of bootstraped parameters 

[~,~,y1,y2] = f_solution_bacillus(xInitial,parameters,tExp,meanPar,par2names,caso); % meanPar are the bootstraped parameters, par2names are the names of the bootstraped parameters  

ysim = {y1(:,col) y2(:,col)};
tAll = {t t};
yExpAll = {yExp1 yExp2};
tExpAll = {tExp tExp};
[RMSE2, rSquared] = f_RMSE(tExpAll,yExpAll,ysim,tAll,compounds); 

