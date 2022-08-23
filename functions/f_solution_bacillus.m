function [t,y,y1,y2] = f_solution_bacillus(xInitial,parameters,tExp,parameterReached,par2names,caso)
% swaps out parameters according to par2names and parameterReached and runs
% the simulation according to the new parameter values 

 nneg = 1:length(xInitial);
 options=odeset('RelTol',1e-4,'NonNegative',nneg);
%  options=odeset('RelTol',1e-4);

% km_su

if any(strcmp(par2names,'km_su'))
    km_su=parameterReached(strcmp(par2names,'km_su'));
    parameters.parValues(strcmp(parameters.parAbb,'km_su'))=km_su;
      parameters.parValues(strcmp(parameters.parAbb,'km_ye'))=km_su;
end

% km_lac

if any(strcmp(par2names,'km_lac'))
    km_lac=parameterReached(strcmp(par2names,'km_lac'));
    parameters.parValues(strcmp(parameters.parAbb,'km_lac'))=km_lac;
end

% km_ye

if any(strcmp(par2names,'km_ye'))
    km_ye=parameterReached(strcmp(par2names,'km_ye'));
    parameters.parValues(strcmp(parameters.parAbb,'km_ye'))=km_ye;
end

% KS_su
if any(strcmp(par2names,'KS_su'))
    KS_su=parameterReached(strcmp(par2names,'KS_su'));
    parameters.parValues(strcmp(parameters.parAbb,'KS_su'))=KS_su;
end

% KS_ye
if any(strcmp(par2names,'KS_ye'))
    KS_ye=parameterReached(strcmp(par2names,'KS_ye'));
    parameters.parValues(strcmp(parameters.parAbb,'KS_ye'))=KS_ye;
end

% KS_lac
if any(strcmp(par2names,'KS_lac'))
    KS_lac=parameterReached(strcmp(par2names,'KS_lac'));
    parameters.parValues(strcmp(parameters.parAbb,'KS_lac'))=KS_lac;
end

% KI_VFA
if any(strcmp(par2names,'KI_VFA'))
    KI_VFA=parameterReached(strcmp(par2names,'KI_VFA'));
    parameters.parValues(strcmp(parameters.parAbb,'KI_VFA'))=KI_VFA;
end

% lag_time
if any(strcmp(par2names,'lag_time'))
    lag_time=parameterReached(strcmp(par2names,'lag_time'));
    parameters.parValues(strcmp(parameters.parAbb,'lag_time'))=lag_time;
end

if any(ismember(par2names, {'Ysu', 'Ylac','Yye'}))
    if any(strcmp(par2names,'Ysu'))
        Ysu=parameterReached(strcmp(par2names,'Ysu'));
        parameters.parValues(strcmp(parameters.parAbb,'Ysu'))=Ysu;
        parameters.parValues(strcmp(parameters.parAbb,'Yye'))=Ysu;
    end
    
    if any(strcmp(par2names,'Ylac'))
        Ylac=parameterReached(strcmp(par2names,'Ylac'));
        parameters.parValues(strcmp(parameters.parAbb,'Ylac'))=Ylac;
    end
    
    if any(strcmp(par2names,'Yye'))
        Yye=parameterReached(strcmp(par2names,'Yye'));
        parameters.parValues(strcmp(parameters.parAbb,'Yye'))=Yye;
    end
        
    stoM = f_stoichiometric_matrix(parameters.stoTable,parameters);
    parameters.stoM=stoM;
end

% K_YE
if any(strcmp(par2names,'K_YE'))
    K_YE=parameterReached(strcmp(par2names,'K_YE'));
    pos = find(strcmp('Yeast extract',parameters.compounds));
    xInitial(:,pos)=K_YE*xInitial(:,pos);
end

switch caso
    case 'montecarlo'
        t=linspace(0,tExp(end),50);
        t = union(t,tExp);
    case 'normal'
        t=linspace(0,tExp(end),50);
end

[t,y1] = ode15s(@f_mass_balances,t,xInitial(1,:),options,parameters);
[~,y2] = ode15s(@f_mass_balances,t,xInitial(2,:),options,parameters);

y = [y1; y2];

end

