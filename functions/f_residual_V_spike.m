function [resVector] = f_residual_V_spike(par2estimate,tExp,yExp,xInitial,parameters,par2names,caso)
% Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% October 2021.Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.
%    Caculates the residuals of the Veillionella monoculture model

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



%% Calculate residual
% options 
nneg = 1:length(xInitial);
% important put MaxStep to 0.1 otherwise spikes wil be missed
options=odeset('RelTol',1e-4,'NonNegative',nneg,'MaxStep',0.1);
% positions outputs of interest 
pos_pro = find(strcmp(parameters.compoundAbb,'Spro'));
pos_ace = find(strcmp(parameters.compoundAbb,'Sac'));
pos_lac = find(strcmp(parameters.compoundAbb,'Slac'));
pos_BM_lac = find(strcmp(parameters.compoundAbb,'Xlac'));

%get spike data 
spikeData = parameters.spikeFlowrate;
volumes = spikeData.volumeSpikes;
concHelp = spikeData.compositionSpikes;
conc = sum(concHelp(:,1)); %to get rid of the array

%solutions 
[tt,y1]=ode15s(@f_mass_balances,tExp,xInitial(1,:),options,parameters);
[tt,y2]=ode15s(@f_mass_balances,tExp,xInitial(2,:),options,parameters);

model_results_1 = y1(:,[pos_lac,pos_pro,pos_ace ,pos_BM_lac]);
model_results_2 = y2(:,[pos_lac,pos_pro,pos_ace ,pos_BM_lac]);

model_results = [model_results_1 model_results_2];

residual = (model_results - yExp);
residual(isnan(residual))=0;

switch caso
    case 'normal' % nomalise residuals 
        normal = range(yExp,1);
        residual = bsxfun(@rdivide, residual, normal); %very low values to be overinfluencial in the residual minimisation.
        resVector = residual(:);
        
    case 'bootstrap'
        resVector = residual(:);
end

