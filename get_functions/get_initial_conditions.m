function xInitial=get_initial_conditions(parameters,reactor,yExp)
% Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% October 2021.Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.
%    Defines the initial conditions in a reactor
%yExp can be the data of the experiments where the initial conditions can
%be extracted from or a string value specifing a generic initial condition



if ischar(yExp)
    nCompounds=parameters.simulation(1);        %Get the number of compounds
    indexCompounds=parameters.indexCompounds;   %Get the indexes that identifIni the compounds
    
    
    mInitial=reactor.V;     %kg
    
    SsuInitial = 15;         %kg COD m-3
    SlacInitial =0;        %kg COD m-3
    SyeInitial = 0;
    SproInitial = 0;        %kg COD m-3
    SacInitial = 0;     %kg COD m-3
    Sh2Initial = 0;         %kg COD m-3
    SticInitial = 0.1;        %kmole or ??? kg COD m-3
    SInInitial = 0.1;
    XcInitial = 0;                  %kg COD m-3
    XsuInitial = 0.035;         %kg COD m-3
    XlacInitial = 0;        %kg COD m-3
    ScatInitial = 0.005;        %kmole m-3
    SanInitial = 0.040;         %kmole m-3      %kmole m-3
    
else
    flagCase = parameters.flagCase;
    switch flagCase
        case 'batchBacilus'
            
            nCompounds=parameters.simulation(1);        %Get the number of compounds
            indexCompounds=parameters.indexCompounds;   %Get the indexes that identifIni the compounds
            %xInitial=zeros(1,nCompounds);
            %These are the initial conditions. Note that the volume cannot be zero since it would
            %lead to NaN of the concentrations (C = M /V)
            %mInitial= 3.4e6;        %kg
            mInitial=reactor.V;     %kg
            SsuInitial = yExp(1,1);         %kg COD m-3
            SlacInitial = yExp(1,2);        %kg COD m-3
            SyeInitial = 10*1.37;          % kgCOD m-3
            SproInitial = 0;        %kg COD m-3
            SacInitial = 0;         %kg COD m-3
            Sh2Initial = 0;         %kg COD m-3
            SticInitial = 0.1;        %kmole or ??? kg COD m-3
            SInInitial = 0.1;
            XcInitial = 0;          %kg COD m-3
            XsuInitial = yExp(1,3);         %kg COD m-3
            XlacInitial = 0;         %kg COD m-3
            ScatInitial = 0.005;        %kmole m-3
            SanInitial = 0.040;         %kmole m-3      %kmole m-3
            
        case 'batchVellonela'
            
            nCompounds=parameters.simulation(1);        %Get the number of compounds
            indexCompounds=parameters.indexCompounds;   %Get the indexes that identifIni the compounds
            %xInitial=zeros(1,nCompounds);
            %These are the initial conditions. Note that the volume cannot be zero since it would
            %lead to NaN of the concentrations (C = M /V)
            %mInitial= 3.4e6;        %kg
            mInitial=reactor.V;     %kg
            SsuInitial = 0;         %kg COD m-3
            SlacInitial = yExp(1,1); %kg COD m-3
            SyeInitial = 2*1.37;
            SproInitial = yExp(1,2);        %kg COD m-3
            SacInitial = yExp(1,3);     %kg COD m-3
            Sh2Initial = 0;         %kg COD m-3
            SticInitial = 0.1;        %kmole or ??? kg COD m-3
            SInInitial = 0.1;
            XcInitial = 0;          %kg COD m-3
            XsuInitial = 0;         %kg COD m-3
            XlacInitial = yExp(1,4);         %kg COD m-3
            ScatInitial = 0.005;        %kmole m-3
            SanInitial = 0.040;         %kmole m-3      %kmole m-3
            
        case 'SpikeVeillonella'
            
            nCompounds=parameters.simulation(1);        %Get the number of compounds
            indexCompounds=parameters.indexCompounds;   %Get the indexes that identifIni the compounds
            %xInitial=zeros(1,nCompounds);
            %These are the initial conditions. Note that the volume cannot be zero since it would
            %lead to NaN of the concentrations (C = M /V)
            %mInitial= 3.4e6;        %kg
            mInitial=reactor.V;     %kg
            SsuInitial = 0;         %kg COD m-3
            SlacInitial = yExp(1,1); %kg COD m-3
            SyeInitial = 2*1.37;      %kg COD m-3
            SproInitial = yExp(1,2);        %kg COD m-3
            SacInitial = yExp(1,3);     %kg COD m-3
            Sh2Initial = 0;         %kg COD m-3
            SticInitial = 0.1;        %kmole or ??? kg COD m-3
            SInInitial = 0.1;
            XcInitial = 0;          %kg COD m-3
            XsuInitial = 0;         %kg COD m-3
            XlacInitial = yExp(1,4);         %kg COD m-3
            ScatInitial = 0.005;        %kmole m-3
            SanInitial = 0.040;         %kmole m-3      %kmole m-3
            
        case 'CoCulture'
            
            nCompounds=parameters.simulation(1);        %Get the number of compounds
            indexCompounds=parameters.indexCompounds;   %Get the indexes that identifIni the compounds
            %xInitial=zeros(1,nCompounds);
            %These are the initial conditions. Note that the volume cannot be zero since it would
            %lead to NaN of the concentrations (C = M /V)
            %mInitial= 3.4e6;        %kg
            mInitial=reactor.V;     %kg
            
            SsuInitial = yExp(1,1);         %kg COD m-3
            SlacInitial = yExp(1,2);        %kg COD m-3
            %SyeInitial = 10*1.37*0.23;          % kgCOD m-3
            SyeInitial = 0;
            SproInitial = yExp(1,3);        %kg COD m-3
            SacInitial = yExp(1,4);     %kg COD m-3
            Sh2Initial = 0;         %kg COD m-3
            SticInitial = 0.1;        %kmole or ??? kg COD m-3
            SInInitial = 0.1;
            XcInitial = 0;
            
            if yExp(1,5) <= 0   % give initial condition if not available 
                XsuInitial = 0.03;  %kg COD m-3 based on the previous experiments
            else
                XsuInitial = yExp(1,5);  %kg COD m-3
            end
            
                   
            XlacInitial = 0;        %kg COD m-3
            ScatInitial = 0.005;        %kmole m-3
            SanInitial = 0.040;         %kmole m-3      %kmole m-3
            
            
            
    end
end
xInitialAll=[mInitial,SsuInitial,SlacInitial,SyeInitial,SproInitial,SacInitial,Sh2Initial,SticInitial,SInInitial,XcInitial,XsuInitial,XlacInitial,ScatInitial,SanInitial];
xInitial=xInitialAll(indexCompounds);
