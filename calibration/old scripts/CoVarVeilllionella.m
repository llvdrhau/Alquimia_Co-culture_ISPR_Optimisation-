%% co-Variance of veillionella spike experiments
clear
close all
clc

%% Define experiment 
flagCase = 'SpikeVeillonella';
kinetics = 'seperateProLac';
kineticFlag = 'CoCulture' ; 
sheetData = 'V34';
control = [0 1];

% REACTOR 
reactor.V=1; %L
reactor.pH(1)=1;             % Define whether pH is kept constant or not
reactor.pH(2)=7;             % pH value if it's controlled

%% prep parameter structure 
% parameters=get_parameters(reactor,flagCase,kineticFlag)

load('parametersCorrect','parCorrect') % parameters of estimated veillionella spikes experiments 
parameters = parCorrect; 
parameters.flagCase = flagCase;
parameters.kineticInhibition = kinetics;


%% Initial conditions 
[dataExp1,tExp] = get_exp_data('V3');
[dataExp2,~] = get_exp_data('V4');

tEnd = tExp(end);
parameters.tEnd = tEnd; 
xInitial=get_initial_conditions(parameters,reactor,dataExp1);

%% prep the parameters of interest to simulate 
load('bootstrapResultsVeillionella.mat','parameterReported_dev')
parName = parameterReported_dev(1,:);
parName{end+1} = 'Yyev';
parName{end+1} = 'KS_lac';
parName{end+1} = 'KS_yev';

parVal = zeros(1,length(parName));
for i = 1:length(parName)
    pos = strcmp(parameters.parAbb,parName(i));
    parVal(i) = parameters.parValues(pos);
end

% attitional tweeking: make sure the lag time is zero 
parameters.parValues(strcmp(parameters.parAbb,'lag_time')) = 0;

%% get Covariance matrix 
epsilon = 1e-7;
sigma = ones(1,length(parameters.compoundAbb));  % Scaling of outputs (Q matirx for FIM)
% more or less max values 
sigma(3) = 14; %lactate 
sigma(5) = 20; %pro
sigma(6) = 5.5; %ace
sigma(12) = 0.45; %lac BM  

FIM = f_fischerInformationMatrix(parameters,parVal,parName,xInitial,control,sigma,sheetData);
C = pinv(FIM);  

%C = (FIM2)^(-1);   
%C = inv(FIM2);            % warning close to singular!!! maybe use pinv? 

corr = zeros(length(C));
sigRow = diag(C);
sigCol = diag(C);

for i= 1:length(C)      %row
    for j = 1:length(C) %col
       corr(i,j)= C(i,j)/(sqrt(sigRow(i)*sigCol(j)));      
    end 
end


%% print to command window 
% preallocate 
n = length(parName); 
PR = cell(n+1,n+1);
PR(1,2:end) = parName;
PR(2:end,1) = parName; 
PR(2:end,2:end) = num2cell(corr); 

disp(PR)



%% get the forward sensitivity of the simulation outputs  
allsens = [];
p = length(parName); 
for i=1:p
    sens = f_forwaardSensitivity(parameters,parName{i}, parVal(i), xInitial,control,sheetData);
    allsens = cat(3, allsens, sens);                  
end

%% plot sensitivities
abb = {'Slac','Spro','Sac','Xlac'};
s = length(abb);
pos2 = zeros(1,s);
for i= 1:s
    pos2(i) = find(strcmp(parameters.compoundAbb,abb(i)));
end


states = parameters.compounds(pos2);
t = allsens(:,1,1);
poslegend = zeros(1,length(parVal));

for k = 1:s
    y = squeeze(allsens(:,pos2(k)+1,:));  % +1 because first pos is time 
    ycheck = y.*parVal;
    subplot(2,2,k)
    hold on
    for i = 1:length(parVal)
        yState = y(:,i)*parVal(i);
        poslegend(i) = 1;
        plot(t,yState)
        
    end
    legend(parName(logical(poslegend)))
    xlabel('Time in h')
    ylabel('Relative sensitivity')
    title(char(states{k}))
    poslegend = zeros(1,length(parVal));
    hold off
end




