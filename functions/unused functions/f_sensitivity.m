function [dydp,dyNDdp] = f_sensitivity(y,idx,parameterReached,xInitial,parameters,tExp,par2names)
nneg = 1:length(xInitial);
options=odeset('RelTol',1e-4,'NonNegative',nneg);
% options=odeset('RelTol',1e-4);

% sensitvity of outputs 
pos_glu = find(strcmp(parameters.compoundAbb,'Sglu'));
pos_lac = find(strcmp(parameters.compoundAbb,'Slac'));
pos_bm = find(strcmp(parameters.compoundAbb,'Xsu'));
pos_ye = find(strcmp(parameters.compoundAbb,'Sye'));
pos_tan = find(strcmp(parameters.compoundAbb,'SIn'));

idY = [pos_glu pos_lac pos_bm pos_ye pos_tan];

yND = repmat(mean(y),size(tExp,1)*2,1);   % Matrix containing the vector of output mean values repeated t times.

dev = [1 -1];
for i=1:length(idx)                      % Building the diferential sensibility matrix (one column for each parameter)
    for j=1:2
        parIteration = parameterReached;
        parIteration(i) = parIteration(i)*(1+dev(j)*1e-3);
        xInitialUnc = xInitial;
        
        % km_su
        
        if any(strcmp(par2names,'km_su'))
            km_su=parIteration(strcmp(par2names,'km_su'));
            parameters.parValues(strcmp(parameters.parAbb,'km_su'))=km_su;
            parameters.parValues(strcmp(parameters.parAbb,'km_ye'))=km_su;
        end
        
        % km_lac
        
        if any(strcmp(par2names,'km_lac'))
            km_lac=parIteration(strcmp(par2names,'km_lac'));
            parameters.parValues(strcmp(parameters.parAbb,'km_lac'))=km_lac;
        end
        
        % km_ye
        
        if any(strcmp(par2names,'km_ye'))
            km_ye=parIteration(strcmp(par2names,'km_ye'));
            parameters.parValues(strcmp(parameters.parAbb,'km_ye'))=km_ye;
        end
        
        % KS_su
        if any(strcmp(par2names,'KS_su'))
            KS_su=parIteration(strcmp(par2names,'KS_su'));
            parameters.parValues(strcmp(parameters.parAbb,'KS_su'))=KS_su;
        end
        
        % KS_ye
        if any(strcmp(par2names,'KS_ye'))
            KS_ye=parIteration(strcmp(par2names,'KS_ye'));
            parameters.parValues(strcmp(parameters.parAbb,'KS_ye'))=KS_ye;
        end
        
        % KS_lac
        if any(strcmp(par2names,'KS_lac'))
            KS_lac=parIteration(strcmp(par2names,'KS_lac'));
            parameters.parValues(strcmp(parameters.parAbb,'KS_lac'))=KS_lac;
        end
        
        % KI_VFA
        if any(strcmp(par2names,'KI_VFA'))
            KI_VFA=parIteration(strcmp(par2names,'KI_VFA'));
            parameters.parValues(strcmp(parameters.parAbb,'KI_VFA'))=KI_VFA;
        end
        
        % lag_time
        if any(strcmp(par2names,'lag_time'))
            lag_time=parIteration(strcmp(par2names,'lag_time'));
            parameters.parValues(strcmp(parameters.parAbb,'lag_time'))=lag_time;
        end
        
        % Yields
        
        if any(ismember(par2names, {'Ysu', 'Ylac','Yye'}))
            if any(strcmp(par2names,'Ysu'))
                Ysu=parIteration(strcmp(par2names,'Ysu'));
                parameters.parValues(strcmp(parameters.parAbb,'Ysu'))=Ysu;
                parameters.parValues(strcmp(parameters.parAbb,'Yye'))=Ysu;
            end
            
            if any(strcmp(par2names,'Ylac'))
                Ylac=parIteration(strcmp(par2names,'Ylac'));
                parameters.parValues(strcmp(parameters.parAbb,'Ylac'))=Ylac;
            end
            
            if any(strcmp(par2names,'Yye'))
                Yye=parIteration(strcmp(par2names,'Yye'));
                parameters.parValues(strcmp(parameters.parAbb,'Yye'))=Yye;
            end
            
            stoM = f_stoichiometric_matrix(parameters.stoTable,parameters);
            parameters.stoM=stoM;
        end
        
        % K_YE
        if any(strcmp(par2names,'K_YE'))
            K_YE=parIteration(strcmp(par2names,'K_YE'));
            pos = find(strcmp('Yeast extract',parameters.compounds));
            xInitialUnc(:,pos)=K_YE*xInitial(:,pos);
        end
        
        
        
        if j==1
            [~,y1] = ode15s(@f_mass_balances,tExp,xInitialUnc(1,:),options,parameters);
            [~,y2] = ode15s(@f_mass_balances,tExp,xInitialUnc(2,:),options,parameters);
            yPlus=[y2; y1];
        end
        if j==2
            [~,y1] = ode15s(@f_mass_balances,tExp,xInitialUnc(1,:),options,parameters);
            [~,y2] = ode15s(@f_mass_balances,tExp,xInitialUnc(2,:),options,parameters);
            yMinus=[y1; y2];
        end
        
        
    end
    % Building normalised dy/dp matrix
    dp = abs(parameterReached(i) - parIteration(i));
    dy = (yPlus(:,idY)-yMinus(:,idY))/(2*dp);                              % dy/dp
    dyND = dy./yND(:,idY)*parameterReached(i);                         % Normalising by multiplication by parameter/mean(y)
    dydp(:,i) = dy(:);                                                     % Building dy/dp matrix. Each column is one parameters. All the data is expressed as a vector. (data x parameters)
    dyNDdp(:,i) = dyND(:);
end

dydp(isnan(dydp)) = 0;
dyNDdp(isnan(dyNDdp)) = 0;
end

