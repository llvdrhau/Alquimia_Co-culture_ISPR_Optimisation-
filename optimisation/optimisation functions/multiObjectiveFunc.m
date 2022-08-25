function [costFunc,parametersOut,t,y,reservoir] = multiObjectiveFunc(decisionVar,parameters,xInitial)
%% Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% October 2021.Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.
%    This script calculates the yield and productivity of proprionate as 
%    objective variables 

%% Select correct data 

if isfield(parameters,'gluSwitch') && parameters.gluSwitch
    % SPIKE Concentration
    parameters.spikeFlowrate.compositionSpikes(2:end,1) = decisionVar(1);
    % SPIKE VEILLIONELLA
    parameters.spikeFlowrate.tSpikes(1) = decisionVar(2);
    % FIRST REED
    startPos = decisionVar(3);   
    % SPACING 
    spacing =  decisionVar(4:end);
else
    % SPIKE VEILLIONELLA
    if isfield(parameters,'exp') && strcmp(parameters.exp,'BV5')
        %the second input is the spike for veilionella
        parameters.spikeFlowrate.tSpikes(2) = decisionVar(1);
    else
        parameters.spikeFlowrate.tSpikes(1) = decisionVar(1);
    end
    % FIRST REED
    startPos = decisionVar(2);
    % SPACING 
    spacing =  decisionVar(3:end);
end

%% Re-position the reed/spkie data 
parameters.Reed.tStart(1) = startPos;
parameters.Reed.tEnd(1) = startPos+2;    % stop 2 hours later

startPosNew = startPos;
posSpike = 1; 
posReed = 1; 
% updating the folwing spike/reed times 
for i = 1:length(spacing)
    if rem(i,2) % if odd update spike time
        if isfield(parameters,'exp') && strcmp(parameters.exp,'BV5') % in case of bv5
            parameters.spikeFlowrate.tSpikes(posSpike+2) = startPosNew + spacing(i);
            posSpike = posSpike +1; 
        else % for any other case (the generic one) 
            parameters.spikeFlowrate.tSpikes(posSpike+1) = startPosNew + spacing(i);
             posSpike = posSpike +1;
        end   
    else % if i even update REED time 
        parameters.Reed.tStart(posReed+1) = startPosNew + spacing(i);
        parameters.Reed.tEnd(posReed+1) = startPosNew + spacing(i) + 2; % stop 2 hours later
        posReed = posReed +1; 
    end
    startPosNew = startPosNew +spacing(i);
end

%% End of sim when last reed has stoped 
tEndSim =parameters.Reed.tEnd(end);
   
%% Run ODE
nneg = 1:length(xInitial);
options=odeset('RelTol',1e-4,'NonNegative',nneg,'MaxStep',0.1);
tSim1 = linspace(0,tEndSim,1000); %300 steps original 

%for steps to include the spike and REED times steps
tSim2 = union(tSim1, parameters.Reed.tStart);
tSim3 = union(tSim2, parameters.Reed.tEnd);
tSim = union(tSim3,  parameters.spikeFlowrate.tSpikes);

% ode 
[t,y] = ode45(@f_mass_balances,tSim,xInitial,options,parameters);

% reservoir data
pos_glu = strcmp(parameters.compoundAbb,'Sglu');
pos_pro = strcmp(parameters.compoundAbb,'Spro');
%pos_lac = strcmp(parameters.compoundAbb,'Slac');



%% Reservoir compounds from extraction
reservoir = zeros(length(t),14); 
helpReservoir = zeros(1,14); 
volReservoir = 1; %L 

for i = 2:length(t)    
    states = y(i,:)';
    [~,~,transferReed]=f_mass_balances(t(i),states,parameters); 
    helpReservoir = helpReservoir + transferReed'*(t(i)-t(i-1))*states(1)/volReservoir;
    reservoir(i,:)= helpReservoir;
    
end


%% define the two objectyive functions 

gluAdded  = sum(parameters.spikeFlowrate.volumeSpikes.* ... % in gCOD
    parameters.spikeFlowrate.compositionSpikes(:,1)) ...  % added by spike
    + y(end,pos_glu)*y(end,1) + ...     % glucose left over
    + y(1,pos_glu)*y(1,1);    %present in the begining of the batch

propMade = reservoir(end,pos_pro);

%decisionVar
productivity = -reservoir(end,pos_pro)/t(end);
yieldActual = -propMade/gluAdded ;
yield  = max(yieldActual,-0.7);

costFunc(1) = productivity;
costFunc(2) = yield;

parametersOut = parameters; 
end


