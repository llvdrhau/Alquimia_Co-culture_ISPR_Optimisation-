function [t,y,y1,y2] = f_solution(xInitial,parameters,tExp,parameterReached,par2names,caso)

 nneg = 1:length(xInitial);
 options=odeset('RelTol',1e-4,'NonNegative',nneg,'MaxStep',0.1);
%  options=odeset('RelTol',1e-4);


% change each structure parameters 
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
    parameters.parValues(strcmp(parameters.parAbb,'km_ye'))=km_lac;
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

% KI_lac
if any(strcmp(par2names,'KI_lac'))
    KI_lac=parameterReached(strcmp(par2names,'KI_lac'));
    parameters.parValues(strcmp(parameters.parAbb,'KI_lac'))=KI_lac;
end

% KI_pro
if any(strcmp(par2names,'KI_pro'))
    KI_pro=parameterReached(strcmp(par2names,'KI_pro'));
    parameters.parValues(strcmp(parameters.parAbb,'KI_pro'))=KI_pro;
end


% lag_time
if any(strcmp(par2names,'lag_time'))
    lag_time=parameterReached(strcmp(par2names,'lag_time'));
    parameters.parValues(strcmp(parameters.parAbb,'lag_time'))=lag_time;
end

% decay 
if any(strcmp(par2names,'kdec_Xlac'))
    kdec_Xlac =parameterReached(strcmp(par2names,'kdec_Xlac'));
    parameters.parValues(strcmp(parameters.parAbb,'kdec_Xlac'))=kdec_Xlac;
end

% Yields
if any(ismember(par2names, {'Ysu', 'Ylac','Yye'}))
    if any(strcmp(par2names,'Ysu'))
        Ysu=parameterReached(strcmp(par2names,'Ysu'));
        parameters.parValues(strcmp(parameters.parAbb,'Ysu'))=Ysu;
        %parameters.parValues(strcmp(parameters.parAbb,'Yye'))=Ysu;   %QUESTOIN why make the yield of yeast and glucose the same here?
    end
    
    if any(strcmp(par2names,'Ylac'))
        Ylac=parameterReached(strcmp(par2names,'Ylac'));
        parameters.parValues(strcmp(parameters.parAbb,'Ylac'))=Ylac;
    end
    
    if any(strcmp(par2names,'Yye'))
        Yye=parameterReached(strcmp(par2names,'Yye'));
        parameters.parValues(strcmp(parameters.parAbb,'Yye'))=Yye;
    end
   
    %update Stoichimetric matrix
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
        t=linspace(0,tExp(end),100);
        t = union(t,tExp);
        if parameters.control(2)   % add time stamps of tSpikes
            spikeFlowrate = parameters.spikeFlowrate;
            tSpikes = spikeFlowrate.tSpikes;
            extraTimeFramesUpper = tSpikes +0.2;
            extraTimeFramesLower = tSpikes -0.2;
            extraTimeFrames = [tSpikes, extraTimeFramesLower, extraTimeFramesUpper];
            t = union(t,extraTimeFrames);
        end
             
end

switch parameters.flagCase
    case 'SpikeVeillonella'
   [t,y] = ode15s(@f_mass_balances,t,xInitial,options,parameters);
        
    otherwise
    [t,y] = ode15s(@f_mass_balances,t,xInitial,options,parameters);
    
end

% 
% 
% switch parameters.flagCase
%     case 'SpikeVeillonella'
%        switch caso 
%            case 'normal'
%         %get spike data
%         spikeData = parameters.spikeFlowrate;
%         volumes = spikeData.volumeSpikes;
%         concHelp = spikeData.compositionSpikes;
%         conc = sum(concHelp(:,1)); %to get rid of the array
%         
%         [t1,ys1]=ode15s(@f_mass_balances,tExp(1:4),xInitial(1,:),options,parameters);
%         Initial_Spike = ys1(end,:);
%         Initial_Spike(3) =  Initial_Spike(3) + conc*volumes(1)/1; %lactate spike
%         
%         
%         [t2,ys2]=ode15s(@f_mass_balances,tExp(5:7),Initial_Spike,options,parameters);
%         Initial_Spike = ys2(end,:);
%         Initial_Spike(3) =  Initial_Spike(3) + conc*volumes(2)/1; %lactate spike
%         
%         
%         [t3,yend]=ode15s(@f_mass_balances,tExp(8:end),Initial_Spike,options,parameters);
%         y =  [ys1; ys2; yend];
%         t = [t1; t2; t3];
%         case 'montecarlo'
%             [t,y] = ode15s(@f_mass_balances,t,xInitial,options,parameters);
%        end 
%     otherwise
%         [t,y] = ode15s(@f_mass_balances,t,xInitial,options,parameters);
% end
% 



end

