function [resVector,t,y] = f_residual_BV5(par2estimate,tExp,yExp,xInitial,parameters,par2names,caso)

%% This script calculates the objective for parameter estimation
% Lucas Van der Hauwaert
% Department of Chemical Engineering, USC
% May 2021

%% Changing parameters in the model
% km_su

if any(strcmp(par2names,'km_su'))
    km_su=par2estimate(strcmp(par2names,'km_su'));
    parameters.parValues(strcmp(parameters.parAbb,'km_su'))=km_su;
    parameters.parValues(strcmp(parameters.parAbb,'km_ye'))=km_su;  %
     %WHY??
end

% km_lac

if any(strcmp(par2names,'km_lac'))
    km_lac=par2estimate(strcmp(par2names,'km_lac'));
    parameters.parValues(strcmp(parameters.parAbb,'km_lac'))=km_lac;
    parameters.parValues(strcmp(parameters.parAbb,'km_ye'))=km_lac;  %
end

% KS_su
if any(strcmp(par2names,'KS_su'))
    KS_su=par2estimate(strcmp(par2names,'KS_su'));
    parameters.parValues(strcmp(parameters.parAbb,'KS_su'))=KS_su;
end

% KS_yev
if any(strcmp(par2names,'KS_yev'))
    KS_yev=par2estimate(strcmp(par2names,'KS_yev'));
    parameters.parValues(strcmp(parameters.parAbb,'KS_yev'))=KS_yev;
end

% KS_yeb
if any(strcmp(par2names,'KS_yeb'))
    KS_yeb=par2estimate(strcmp(par2names,'KS_yeb'));
    parameters.parValues(strcmp(parameters.parAbb,'KS_yev'))=KS_yeb;
end

% KS_lac
if any(strcmp(par2names,'KS_lac'))
    KS_lac=par2estimate(strcmp(par2names,'KS_lac'));
    parameters.parValues(strcmp(parameters.parAbb,'KS_lac'))=KS_lac;
end

% KI_VFA
if any(strcmp(par2names,'KI_VFA'))
    KI_VFA=par2estimate(strcmp(par2names,'KI_VFA'));
    parameters.parValues(strcmp(parameters.parAbb,'KI_VFA'))=KI_VFA;
end

% KI_lac
if any(strcmp(par2names,'KI_lac'))
    KI_lac=par2estimate(strcmp(par2names,'KI_lac'));
    parameters.parValues(strcmp(parameters.parAbb,'KI_lac'))=KI_lac;
end

% KI_pro
if any(strcmp(par2names,'KI_pro'))
    KI_pro=par2estimate(strcmp(par2names,'KI_pro'));
    parameters.parValues(strcmp(parameters.parAbb,'KI_pro'))=KI_pro;
end


% lag_time
if any(strcmp(par2names,'lag_time'))
    lag_time=par2estimate(strcmp(par2names,'lag_time'));
    parameters.parValues(strcmp(parameters.parAbb,'lag_time'))=lag_time;
end

% decay 
if any(strcmp(par2names,'kdec_Xlac'))
    kdec_Xlac =par2estimate(strcmp(par2names,'kdec_Xlac'));
    parameters.parValues(strcmp(parameters.parAbb,'kdec_Xlac'))=kdec_Xlac;
end

% Yields
if any(ismember(par2names, {'Ysu', 'Ylac','Yye'}))
    if any(strcmp(par2names,'Ysu'))
        Ysu=par2estimate(strcmp(par2names,'Ysu'));
        parameters.parValues(strcmp(parameters.parAbb,'Ysu'))=Ysu;
        parameters.parValues(strcmp(parameters.parAbb,'Yye'))=Ysu;
    end
    
    if any(strcmp(par2names,'Ylac'))
        Ylac=par2estimate(strcmp(par2names,'Ylac'));
        parameters.parValues(strcmp(parameters.parAbb,'Ylac'))=Ylac;
    end
    
    if any(strcmp(par2names,'Yye'))
        Yye=par2estimate(strcmp(par2names,'Yye'));
        parameters.parValues(strcmp(parameters.parAbb,'Yye'))=Yye;
    end
        
    stoM = f_stoichiometric_matrix(parameters.stoTable,parameters);
    parameters.stoM=stoM;
end

% K_YE
if any(strcmp(par2names,'K_YE'))
    K_YE=par2estimate(strcmp(par2names,'K_YE'));
    pos = find(strcmp('Yeast extract',parameters.compounds));
    xInitial(:,pos)=K_YE*xInitial(:,pos);
end

% spike of veill, concentration added 
if any(strcmp(par2names,'spike'))
    spike=par2estimate(strcmp(par2names,'spike'));
    parameters.spikeFlowrate.compositionSpikes(2,11) = spike ;     %Spike  Veillionella
   
end



%% Calculate residual

%% find model results corresponing to time frame and length of sim  
tEnd = tExp(end);% h

nIteration=min(15,tEnd*(24*6));    % one time node per 10 min
tSim=linspace(0,tEnd,100);
[nDataExp,nExp]=size(tExp);
tSim = union(tSim,tExp); 
%indexCell = cell(nDataExp,nExp);   % preallocation, use if more than  one
%dataset
for j=1:nExp
    for i=1:nDataExp
        matrixTime = repmat(tSim,1,length(tExp));
        matrixExpTime = repmat(tExp',length(tSim),1);
        [~,indexExpTime] = min(abs(matrixExpTime-matrixTime));
    end
end

%%
nneg = 1:length(xInitial);
options=odeset('RelTol',1e-4,'NonNegative',nneg,'MaxStep',0.1);
%par2estimate;
[t,y] = ode45(@f_mass_balances,tSim,xInitial,options,parameters);

pos_glu = find(strcmp(parameters.compoundAbb,'Sglu'));
pos_pro = find(strcmp(parameters.compoundAbb,'Spro'));
pos_ace = find(strcmp(parameters.compoundAbb,'Sac'));
pos_lac = find(strcmp(parameters.compoundAbb,'Slac'));
%pos_BM_lac = find(strcmp(parameters.compoundAbb,'Xlac'));


model_results = y(:,[pos_glu,pos_lac,pos_pro,pos_ace]);
residual = (yExp - model_results(indexExpTime,:));
 

switch caso
    case 'normal'
        normal = range(yExp,1);
        residual = bsxfun(@rdivide, residual, normal);  %very low values to be overinfluencial in the residual minimisation.
        resVector = residual(:);
        resVector= norm(resVector)
        par2estimate
    case 'bootstrap'
        resVector = residual(:);
end

