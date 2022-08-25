function [resVector] = f_residual_bacillus(par2estimate,tExp,yExp,xInitial,parameters,par2names,caso)
% Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% October 2021.Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.
%    This script calculates the residuals for parameter estimation of the
%    bacillus monoculture model

%% Changing parameters in the model 
% km_su

if any(strcmp(par2names,'km_su'))
    km_su=par2estimate(strcmp(par2names,'km_su'));
    parameters.parValues(strcmp(parameters.parAbb,'km_su'))=km_su;
     parameters.parValues(strcmp(parameters.parAbb,'km_ye'))=km_su;
end

% km_lac

if any(strcmp(par2names,'km_lac'))
    km_lac=par2estimate(strcmp(par2names,'km_lac'));
    parameters.parValues(strcmp(parameters.parAbb,'km_lac'))=km_lac;
end

% km_ye

if any(strcmp(par2names,'km_ye'))
    km_ye=par2estimate(strcmp(par2names,'km_ye'));
    parameters.parValues(strcmp(parameters.parAbb,'km_ye'))=km_ye;
end

% KS_su
if any(strcmp(par2names,'KS_su'))
    KS_su=par2estimate(strcmp(par2names,'KS_su'));
    parameters.parValues(strcmp(parameters.parAbb,'KS_su'))=KS_su;
end

% KS_ye
if any(strcmp(par2names,'KS_ye'))
    KS_ye=par2estimate(strcmp(par2names,'KS_ye'));
    parameters.parValues(strcmp(parameters.parAbb,'KS_ye'))=KS_ye;
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

% lag_time
if any(strcmp(par2names,'lag_time'))
    lag_time=par2estimate(strcmp(par2names,'lag_time'));
    parameters.parValues(strcmp(parameters.parAbb,'lag_time'))=lag_time;
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
nneg = 1:length(xInitial); 
options=odeset('RelTol',1e-4,'NonNegative',nneg);

pos_glu = find(strcmp(parameters.compoundAbb,'Sglu'));
pos_lac = find(strcmp(parameters.compoundAbb,'Slac'));
pos_BM_glu = find(strcmp(parameters.compoundAbb,'Xsu'));

[tt,y1]=ode15s(@f_mass_balances,tExp,xInitial(1,:),options,parameters);
[~,y2]=ode15s(@f_mass_balances,tExp,xInitial(2,:),options,parameters);

model_results_1 = y1(:,[pos_glu,pos_lac,pos_BM_glu]);
model_results_2 = y2(:,[pos_glu,pos_lac,pos_BM_glu]);

model_results = [model_results_1 model_results_2];

residual = (model_results - yExp);
residual(isnan(residual))=0;

switch caso
    case 'normal' % normilise residuals 
        normal = range(yExp,1);
        residual = bsxfun(@rdivide, residual, normal);  
        resVector = residual(:);
        
    case 'bootstrap'
        resVector = residual(:);
end

