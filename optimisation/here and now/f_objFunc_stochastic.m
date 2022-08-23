function [costFunc,yield,productivity] = f_objFunc_stochastic(decisionVar,parameters,xInitial,objective)
% Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% October 2021.Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.
%   This functiom calculates the obective value (10th percentile) of the
%   productivity or yield of proprionate for a stochastic optimisation

%% Select correct data 
if isfield(parameters,'gluSwitch') && parameters.gluSwitch % if glucose also estimated
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
    else
        parameters.Reed.tStart(posReed+1) = startPosNew + spacing(i);
        parameters.Reed.tEnd(posReed+1) = startPosNew + spacing(i) + 2; % stop 2 hours later
        posReed = posReed +1; 
    end
    startPosNew = startPosNew +spacing(i);
end

%% End of sim when last reed has stoped 
tEndSim =parameters.Reed.tEnd(end);

%% %% %%%% %%%% %%%% %%%% %%%% %%
% MONTE CARLO LOOP 

% create distribution matrix 
% load distributions of parameters  
load('bootstrapResultsVeillionella','parameterMatrix','CImatrix') % km_su, Ysu
distVeill  = parameterMatrix;
parVeil = CImatrix(1,:);
load('bootstrapResultsBacillus','parameterMatrix','CImatrix') % km_lac Ylac KI_pro KI_lac kdec_Xlac
distBac = parameterMatrix;
parBac = CImatrix(1,:);

% combined distribution and parameters of interest
distributionMatrix = [distBac, distVeill];
par2names = [parBac,parVeil];
% parameters where uncertainty is evaluated (same oder as the distribution!! important!! )

% find initial values for the montecarlo simulations
pos = zeros(length(par2names),1);  % pre allocation 
for i =1:length(par2names)
    pos(i) = find(strcmp(par2names(i), parameters.parAbb));
end
par2estimateInitial = parameters.parValues(pos);

% run monte carlo functtion 
[~, yield,productivity] = f_montecarlo_here_and_now(par2estimateInitial,par2names,distributionMatrix,xInitial,parameters,tEndSim);
                              
switch objective
    case 'yield'
        costFunc = mean(yield);   % Q: set upper limit??
        %costFunc  = max(costFunc,-0.7);  % so simulation doesn't stretch out
    case 'productivity'
        yieldCnstr = - mean(yield);
        %yieldCnstr = - prctile(yield,10) ;
        if yieldCnstr < 0.63 % extra constraint to yield 
            costFunc = mean(productivity)+1000;
            % costFunc = prctile(productivity,10) +1000;
        else
            costFunc = mean(productivity);
             % costFunc = prctile(productivity,10);
        end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end 



 