function [resVector] = f_residual(par2estimate,tExp,yExp,xInitial,parameters,par2names,caso,solver)

%% This script calculates the residuals for parameter estimation
% Alberte Regueira & Lucas Van der Hauwaert
% Department of Chemical Engineering, USC
% October 2020 

if nargin < 8    
solver = 'fconmin';
end
%% Changing parameters in the model 
% km_su
names = fieldnames(parameters);
resVector = [];
for i = 1:length(names)
    par = parameters.(names{i});
    
    if any(strcmp(par2names,'km_su'))
        km_su=par2estimate(strcmp(par2names,'km_su'));
        par.parValues(strcmp(par.parAbb,'km_su'))=km_su;
        %par.parValues(strcmp(par.parAbb,'km_ye'))=km_su; %hun? 
    end
    
    % km_lac
    
    if any(strcmp(par2names,'km_lac'))
        km_lac=par2estimate(strcmp(par2names,'km_lac'));
        par.parValues(strcmp(par.parAbb,'km_lac'))=km_lac;
    end
    
    % km_ye
    
    if any(strcmp(par2names,'km_ye'))
        km_ye=par2estimate(strcmp(par2names,'km_ye'));
        par.parValues(strcmp(par.parAbb,'km_ye'))=km_ye;
    end
    
    % KS_su
    if any(strcmp(par2names,'KS_su'))
        KS_su=par2estimate(strcmp(par2names,'KS_su'));
        par.parValues(strcmp(par.parAbb,'KS_su'))=KS_su;
    end
    
    % KS_ye
    if any(strcmp(par2names,'KS_ye'))
        KS_ye=par2estimate(strcmp(par2names,'KS_ye'));
        par.parValues(strcmp(par.parAbb,'KS_ye'))=KS_ye;
    end
    
    % KS_lac
    if any(strcmp(par2names,'KS_lac'))
        KS_lac=par2estimate(strcmp(par2names,'KS_lac'));
        par.parValues(strcmp(par.parAbb,'KS_lac'))=KS_lac;
    end
    
    % KI_VFA
    if any(strcmp(par2names,'KI_VFA'))
        KI_VFA=par2estimate(strcmp(par2names,'KI_VFA'));
        par.parValues(strcmp(par.parAbb,'KI_VFA'))=KI_VFA;
    end
    
    % lag_time
    if any(strcmp(par2names,'lag_time'))
        lag_time=par2estimate(strcmp(par2names,'lag_time'));
        par.parValues(strcmp(par.parAbb,'lag_time'))=lag_time;
    end
    
    % Yields
    
    if any(ismember(par2names, {'Ysu', 'Ylac','Yye'}))
        if any(strcmp(par2names,'Ysu'))
            Ysu=par2estimate(strcmp(par2names,'Ysu'));
            par.parValues(strcmp(par.parAbb,'Ysu'))=Ysu;
            par.parValues(strcmp(par.parAbb,'Yye'))=Ysu;
        end
        
        if any(strcmp(par2names,'Ylac'))
            Ylac=par2estimate(strcmp(par2names,'Ylac'));
            par.parValues(strcmp(par.parAbb,'Ylac'))=Ylac;
        end
        
        if any(strcmp(par2names,'Yye'))
            Yye=par2estimate(strcmp(par2names,'Yye'));
            par.parValues(strcmp(par.parAbb,'Yye'))=Yye;
        end
        
        stoM = f_stoichiometric_matrix(par.stoTable,par);
        par.stoM=stoM;
    end
    
    % K_YE
    if any(strcmp(par2names,'K_YE'))
        K_YE=par2estimate(strcmp(par2names,'K_YE'));
        pos = find(strcmp('Yeast extract',par.compounds));
        xInitial(:,pos)=K_YE*xInitial(:,pos);
    end
    
    % run simulation 
    nneg = 1:length(xInitial);
    options=odeset('RelTol',1e-4,'NonNegative',nneg);
    [~,y]=ode15s(@f_mass_balances,tExp{i},xInitial(i,:),options,par);
    
    % Calculate residual
    pos_glu = find(strcmp(par.compoundAbb,'Sglu'));
    pos_lac = find(strcmp(par.compoundAbb,'Slac'));
    pos_BM_glu = find(strcmp(par.compoundAbb,'Xsu'));
    
    model_results = y(:,[pos_glu,pos_lac,pos_BM_glu]);
    residual = (model_results - yExp{i});
    
    switch caso
        case 'normal'  % is this to normalise the residual??
            normal = range(yExp{i},1);
            residual = bsxfun(@rdivide, residual, normal);                             %very low values to be overinfluencial in the residual minimisation.
            res = residual(:);
            
        case 'bootstrap'
            res = residual(:);
    end
    
    resVector = [resVector; res]; 
end

if strcmp(solver,'fconmin')
resVector = norm(resVector);
end
