%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VALIDATE experiment BV5 with bootstrap parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% October 2021.Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.

% Run a simulation of experiment BV5: co-culture with extraction cycles 

clear  
close all 
clc

% BV5
flagCase = 'CoCulture';
kineticFlag = 'CoCulture';
kinetics = 'seperateProLac';
%% Reading experimental data
[data1,Names]=xlsread('expData.xlsx','BV5');

tExp=data1(:,strcmp(Names,'time'));

yExp1(:,1)=data1(:,strcmp(Names,'Glu'));
yExp1(:,2)=data1(:,strcmp(Names,'Lac'));
yExp1(:,3)=data1(:,strcmp(Names,'Pro'));
yExp1(:,4)=data1(:,strcmp(Names,'Ace'));
yExp1(:,5)=data1(:,strcmp(Names,'BM'));

%% set up 
reactor.V=1.3; %L   starting at 1,3 L because we'll assume that the 100 ml Bacillus is in the starting conditions 
reactor.pH(1)=1;             % Define whether pH is kept constant or not
reactor.pH(2)=7;             % pH value if it's controlled

% REED controler and/or spikes 
control = [1,1];

%prep structure parameters
parameters=get_parameters(reactor,flagCase,kineticFlag);
parameters.control = control;
parameters.kineticInhibition = kinetics;

%% spike information BV5
spikeFlowrate = get_spike_data('BV5');
parameters.spikeFlowrate = spikeFlowrate;

%% initial conditions 
xInitial=get_initial_conditions(parameters,reactor,yExp1);

%% Reed data 
[Reed] = get_reed_data('BV5');
% give to structure parameters 
parameters.Reed = Reed; 

%% Manipulate parameters: original tweeking 

parameters.spikeFlowrate.compositionSpikes(2,11) = 1.65 ;     %Spike van Veillionella
parameters.spikeFlowrate.compositionSpikes(5,1) = 1067*1 ;     %last spike 
parameters.spikeFlowrate.compositionSpikes(4,1) = 1067*1 ;     %3th spike 
parameters.spikeFlowrate.compositionSpikes(3,1) = 1067*1 ;     %2nd spike 
parameters.spikeFlowrate.compositionSpikes(1,1) = 1067*1/2 ;     %1st spike 
%xInitial(2) = yExp1(1,1)+5;    

%% use parameter values from bootstrap 

load('parBootStrapVellionella','par2namesV' ,'parameterReportedV')
parameterReached = parameterReportedV;
par2names = par2namesV;
load('parBootStrapBacillus','par2namesB' ,'parameterReportedB')

par2names = [par2names,par2namesB];
parameterReported = [parameterReached,parameterReportedB];

for i = 1:length(par2names)
    pos = strcmp(parameters.parAbb,par2names(i));
    parameters.parValues(pos) = parameterReported(i);
end

% attitional tweeking
parameters.parValues(strcmp(parameters.parAbb,'lag_time')) = 0;
parameters.parValues(strcmp(parameters.parAbb,'kdec_Xlac')) = 0.043; % could use the original parameter for decay from the bootstrap if you reestimate the amount of veil added 
% print to screen
PRallPar(1,:) = par2names;
PRallPar(2,:) = num2cell(parameterReported);
disp('Os parámetros usado son:')
disp(PRallPar)
disp('')

%% Reed start and stop  
for i = 1:length(parameters.Reed.tStart)
    posSt(i) = find(tExp <= parameters.Reed.tStart(i),1,'last');
    posEnd(i) = find(tExp <= parameters.Reed.tEnd(i),1,'last');
end


%% percentage of extraction 
percentPro = 1-yExp1(posEnd,3)./yExp1(posSt,3);
percentAce = 1-yExp1(posEnd,4)./yExp1(posSt,4);
percentLac = 1-yExp1(posEnd,2)./yExp1(posSt,2);

startConcentrations = [yExp1(posSt,3),yExp1(posSt,4),yExp1(posSt,2)];
PR2(1,:) = {'pro','ace', 'Lac'};
PR2(2:5,:) = num2cell(startConcentrations);
disp('Begin concentration of extraction in the broth')
disp(PR2) 
disp('')


percentExtracted = [percentPro,percentAce,percentLac];
PR2(1,:) = {'pro','ace', 'Lac'};
PR2(2:5,:) = num2cell(percentExtracted);
disp('Percent extracted from the broth reactor')
disp(PR2) 
disp('')

meanVFA= [mean(percentPro),mean(percentAce), mean(percentLac)];
PR(1,:) = {'pro','ace', 'Lac'};
PR(2,:) = num2cell(meanVFA);
disp('Mean percentage of extraction')
disp(PR)
disp('')

originalExtractionRate = parameters.Reed.RateMatrix; 
OgRates = [originalExtractionRate(:,4),originalExtractionRate(:,5),originalExtractionRate(:,2)];
PR2(1,:) = {'pro','ace', 'Lac'};
PR2(2:5,:) = num2cell(OgRates);
disp('Extraction rates (g/L/h) determined from data')
disp(PR2) 
disp('')


rLac = mean(parameters.Reed.RateMatrix(:,2));  %lactate
rPro = mean(parameters.Reed.RateMatrix(:,4));   % proprionate 
rAce = mean(parameters.Reed.RateMatrix(:,5));   % Acetate 

meanExtractionRate = [rPro rAce rLac];
PR5(1,:) = {'pro','ace', 'Lac'};
PR5(2,:) = num2cell(meanExtractionRate);
disp('Mean extraction rates (g/L/h) determined from data')
disp(PR5) 
disp('')
%% edit the extraction process with percent extraction ReedSwitch = 1
% not really used in the end 
ReedSwitch = 0;
parameters.ReedSwitch = ReedSwitch;

if ReedSwitch
    parameters.Reed.RateMatrix(:,2) = 0.30; %0.50;  %lactate
    parameters.Reed.RateMatrix(:,4) = 0.30; %0.25;  % proprionate
    parameters.Reed.RateMatrix(:,5) = 0.30;   %0.30;  % Acetate

end

%% edit the extraction process with mean removal rate  
meanExtractionRate = 1;
if meanExtractionRate
    parameters.Reed.RateMatrix(:,2) = mean(parameters.Reed.RateMatrix(:,2));  %lactate
    parameters.Reed.RateMatrix(:,4) = mean(parameters.Reed.RateMatrix(:,4));   % proprionate
    parameters.Reed.RateMatrix(:,5) = mean(parameters.Reed.RateMatrix(:,5));   % Acetate
end


%% run ode15s
tsim = linspace(0,tExp(end),300);
tSpike =  spikeFlowrate.tSpikes;
tOde = union(tsim,tSpike); 
% optioins  
nneg = 1:length(xInitial);
options=odeset('RelTol',1e-4,'NonNegative',nneg,'MaxStep',0.1);
tic
[t,y] = ode15s(@f_mass_balances,tOde,xInitial,options,parameters);
toc
%% R2  coefficient of determination
compoundNames = {'Sglu','Slac' 'Spro' 'Sac' 'BM' };
compoundAbb=parameters.compoundAbb;
glu = y(:,strcmp(compoundAbb,'Sglu'));
lac = y(:,strcmp(compoundAbb,'Slac'));
pro = y(:,strcmp(compoundAbb,'Spro'));
ace = y(:,strcmp(compoundAbb,'Sac'));
Xlac = y(:,strcmp(compoundAbb,'Xlac'));
Xbac = y(:,strcmp(compoundAbb,'Xsu'));

yR2 = [glu,lac,pro,ace,Xlac,Xbac];

[RMSENor, r2] = f_RMSE(tExp,yExp1,yR2,t,compoundNames);

%% plot results 
title = 'BV5 Co-Culture';
plot_CoCulture(tExp,yExp1,t,y,parameters,title,r2)

%% transfer to Reservoir
reservoir = zeros(length(t),14); 
volReservoir = 1; %L 
for i = 2:length(t)
    states = y(i,:)';
    [~,~,transferReed]=f_mass_balances(t(i),states,parameters); 
    reservoir(i,:) = reservoir(i-1,:) + transferReed'*(t(i)-t(i-1))*states(1)/volReservoir;
end
 
lacR = reservoir(:,strcmp(compoundAbb,'Slac'));
proR = reservoir(:,strcmp(compoundAbb,'Spro'));
aceR = reservoir(:,strcmp(compoundAbb,'Sac'));

figure
plot(t,lacR,t,proR,t,aceR)
hold on 
tRV=data1(1:8,strcmp(Names,'timeRV'));
gluRV = data1(1:8,strcmp(Names,'GluRV'));
lacRV = data1(1:8,strcmp(Names,'LacRV'));
proRV = data1(1:8,strcmp(Names,'ProRV'));
aceRV = data1(1:8,strcmp(Names,'AceRV'));
% plot reservoir
plot(tRV,lacRV,'bx',tRV,proRV,'rx',tRV,aceRV,'yx')


legend('Lac','Pro','Ace')
ylabel('concentration g/L') 
xlabel('Time h')

%% RMSE and regresion coeficients of reservoir data 
compoundNames = {'Lac','Pro','Ace'};
yExpReservoir = [lacRV,proRV,aceRV]; 
yRes = [lacR,proR,aceR];
[ nn, r2Res] = f_RMSE(tRV,yExpReservoir,yRes,t,compoundNames);

%% save outputs and figures 
%save('ReservoirData','reservoir','t','tRV','lacRV','tRV','proRV','tRV','aceRV') % to give function so it can be plotted 
path =cd();
saveName = '\ReservoirData';
savePath = append(path,saveName);
save(savePath...
,'reservoir','t','tRV','lacRV','tRV','proRV','tRV','aceRV') % to give function so it can be plotted 
%% alternative plot/ the fancy one 
close all 
plot_CoCulture_reservoir(tExp,yExp1,t,y,parameters,r2,'ReservoirData')
%% calculate yield and productivity 
pos_glu = strcmp(parameters.compoundAbb,'Sglu');
pos_pro = strcmp(parameters.compoundAbb,'Spro');

gluAdded  = sum(parameters.spikeFlowrate.volumeSpikes.* ... % in gCOD
    parameters.spikeFlowrate.compositionSpikes(:,1)) ...  % added by spike
    + y(end,pos_glu)*y(end,1) + ...     % glucose left over
    + y(1,pos_glu)*y(1,1);    %present in the begining of the batch

propMade = reservoir(end,pos_pro);

%decisionVar
productivity = reservoir(end,pos_pro)/t(end);
yieldActual = propMade/gluAdded ;

disp('productivity')
disp(productivity)

disp('yield')
disp(yieldActual)

%% extra snipits of code to check mass balances  
for l=1
%% extraction balance
% %in broth 
% proEx = yExp1(posSt,3)-yExp1(posEnd,3)*1.3;
% aceEx = yExp1(posSt,4)-yExp1(posEnd,4)*1.3;
% lacEx = yExp1(posSt,2)-yExp1(posEnd,2)*1.3;
% Ex = [proEx,aceEx,lacEx];
% % in reservoir 
% Reservoir = [lacRV,proRV,acRV];
% pos = find(~isnan(lacRV),1,'last');
% yRv = Reservoir(1:pos,:);
% 
% counter = 0; 
% for i = 2:2:size(yRv,1)
%     counter = counter+1;
%     proRvEX(counter,1) = yRv(i-1,2)-yRv(i,2)*1;
%     aceRvEX(counter,1) = yRv(i-1,1)-yRv(i,1)*1;
%     lacRvEX(counter,1) = yRv(i-1,3)-yRv(i,3)*1;
% end
% RvEX = [proRvEX aceRvEX lacRvEX];
% RvEX = abs(RvEX); 
% 
% PRres(1,:) = {'pro','ace', 'Lac'};
% PRres(2:5,:) = num2cell(Ex-RvEX);
% disp('difference between extraction from broth and transfer to reservoir in grams') % doesn't make sens for lactate because it's being consumed aswell 
% disp(PRres)
% disp('')

%% save parameters for Parameter estimation 
%save('parEstimationBV5') 
%reservoirBV5 = reservoir; 
%save('BV5Sim','y','t','reservoirBV5')


end
