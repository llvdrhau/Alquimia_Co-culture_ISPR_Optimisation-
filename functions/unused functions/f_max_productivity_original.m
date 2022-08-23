function [productivity,parameters,t,y] = f_max_productivity(decisionVar,parameters,xInitial,idName,idLength)

%% This script calculates the max productivity of proprionate 
% Lucas Van der Hauwaert
% Department of Chemical Engineering, USC
% May 2021

%% Select correct data from id name...
for lkj=1
    nrVar = length(idName);
    posStart = 1;
    posEnd = 0;
    part = cell(1,nrVar); % preallocation 
    for i = 1:nrVar
        if idLength(i) > 1
            posEnd = posEnd + 3  ;
        else
            posEnd = posStart;
        end
        part{i} =  decisionVar(posStart:posEnd);
        posStart = posEnd+1;
        
    end
  
    if any(strcmp(idName,'timing spike'))
        timeSpike= part{strcmp(idName,'timing spike')};
        % update spikes
        parameters.spikeFlowrate.tSpikes(2:end) = timeSpike;
    end
    if any(strcmp(idName,'timing reed'))
        timeReed= part{strcmp(idName,'timing reed')};
        % update time Reed
        parameters.Reed.tStart = timeReed';
        parameters.Reed.tEnd = timeReed'+2; % stop 2 hours later
    end
    if any(strcmp(idName,'end sim'))
        tEnd= part{strcmp(idName,'end sim')};
    else
        tEnd = 65; % if not given
    end
    if any(strcmp(idName,'spike veil'))
        timeVeillionellaSpike= part{strcmp(idName,'spike veil')};
        % update spikes
        parameters.spikeFlowrate.tSpikes(1) = timeVeillionellaSpike;
    end
    
    if any(strcmp(idName,'concFeed'))
        concFeed= part{strcmp(idName,'concFeed')};
        % update Conc Feed
        parameters.spikeFlowrate.compositionSpikes(2:end,1) = concFeed;
    end
    
end

%% The first spike of glu can not happen before the spike of veillionella 
if parameters.spikeFlowrate.tSpikes(2) < parameters.spikeFlowrate.tSpikes(1)
    parameters.spikeFlowrate.tSpikes(2) = parameters.spikeFlowrate.tSpikes(1) + 0.1;
end 

%% Additional restrictions
% first the spike needs to occure then reed extraction
if any(strcmp(idName,'timing reed')) && any(strcmp(idName,'spike veil'))    
    for i= 1:length(timeReed)
        if timeSpike(i) > timeReed(i)
            timeReed(i) = timeSpike(i)+0.1;    
        end
    end
    parameters.Reed.tStart = timeReed';
    parameters.Reed.tEnd = timeReed'+2; % stop 2 hours later
end

% The last reed extraction needs to happen with in the time frame of the
% experiment
if any(strcmp(idName,'timing reed'))
    for i= 1:length(timeReed)
        if timeReed(i)+2 > tEnd
            timeReed(i) = tEnd-2;
        end
    end
    parameters.Reed.tStart = timeReed';
    parameters.Reed.tEnd = timeReed'+2; % stop 2 hours later
end


%      

%%
nneg = 1:length(xInitial);
options=odeset('RelTol',1e-4,'NonNegative',nneg,'MaxStep',0.1);

tSim = linspace(0,tEnd,300); 
%par2estimate;
[t,y] = ode45(@f_mass_balances,tSim,xInitial,options,parameters);
pos_pro = strcmp(parameters.compoundAbb,'Spro');
pos_lac = strcmp(parameters.compoundAbb,'Slac');
reservoir = zeros(length(t),14); 
volReservoir = 1; %L 
for i = 2:length(t)
    states = y(i,:)';
    [~,~,transferReed]=f_mass_balances(t(i),states,parameters); 
    reservoir(i,:) = reservoir(i-1,:) + transferReed'*(t(i)-t(i-1))*states(1)/volReservoir;
end

idName
decisionVar
productivity = -(y(end,pos_pro)*y(end,1)+reservoir(end,pos_pro))/tEnd %+ reservoir(end,pos_lac)/40  % in gCOD/h y(end,1) end volume
 



