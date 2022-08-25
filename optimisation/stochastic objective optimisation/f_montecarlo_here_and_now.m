function [Matrix,yield,productivity] = f_montecarlo_here_and_now(par2estimate,par2names,parameterMatrix,xInitial,parameters,tEndSim)
% Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% October 2021.Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.
%   Monte Carlo here and now: monte carlo simulations for optimisation 
%   under undercertaity 
%   inputs:
%   par2estimate: original parameters values 
%   par2names: cell of parameter names 
%   parameterMatrix: distribution of parameters for sampling 
%   xInitial conditions of simulation
%   parameters: structure for information to run simulation
%   tExp: end of the experment 



 

%for steps to include the spike and REED times steps
tSim1 = linspace(0,tEndSim,1000); %300 steps original, 1000 to reduce errror on reservoir calculations
tSim2 = union(tSim1, parameters.Reed.tStart);
tSim3 = union(tSim2, parameters.Reed.tEnd);
tSim = union(tSim3,  parameters.spikeFlowrate.tSpikes);


nIt = 200;
nPar = length(par2estimate);

%% sampling with inman conover 
Xpar=lhs_empirco(parameterMatrix,nIt);

%For initial values
nVal = 7;
dimensionlessX = lhsdesign(nIt,nVal);    % assume that the initial conditions are uniformly distributed
pos_glu = (strcmp(parameters.compoundAbb,'Sglu'));
pos_lac = (strcmp(parameters.compoundAbb,'Slac'));
pos_pro = (strcmp(parameters.compoundAbb,'Spro'));
pos_ace = (strcmp(parameters.compoundAbb,'Sac'));
pos_bm_su = (strcmp(parameters.compoundAbb,'Xsu'));
pos_bm_lac = (strcmp(parameters.compoundAbb,'Xlac'));
pos_ye = (strcmp(parameters.compoundAbb,'Sye'));
pos_tan = (strcmp(parameters.compoundAbb,'SIn'));


lbVal = [xInitial(:,pos_glu)*0.95  xInitial(:,pos_lac)*1 xInitial(:,pos_pro)*1 ...
    xInitial(:,pos_ace)*1 xInitial(:,pos_bm_su)*0.85 xInitial(:,pos_bm_lac)*0.85 xInitial(:,pos_ye)*0.80];


ubVal = [xInitial(:,pos_glu)*1.05  xInitial(:,pos_lac)*1 xInitial(:,pos_pro)*1 ...
    xInitial(:,pos_ace)*1 xInitial(:,pos_bm_su)*1.15 xInitial(:,pos_bm_lac)*1.15 xInitial(:,pos_ye)*1.20];
diffVal = ubVal - lbVal;

Xval1 = repmat(lbVal(1,:),nIt,1);
Xdiff1 = repmat(diffVal(1,:),nIt,1);
Xval1 = Xval1+Xdiff1.*dimensionlessX;


%preallocation 
matrix1_glu = zeros(length(tSim),nIt);
matrix1_lac = zeros(length(tSim),nIt);
matrix1_pro = zeros(length(tSim),nIt);
matrix1_ace = zeros(length(tSim),nIt);
matrix1_ye = zeros(length(tSim),nIt);
matrix1_bm_su = zeros(length(tSim),nIt);
matrix1_bm_lac = zeros(length(tSim),nIt);
matrix1_tan = zeros(length(tSim),nIt);
productivity = zeros(nIt,1);
yield  = zeros(nIt,1);

for i = 1:nIt
    parSimulation = Xpar(i,:);
    % update structure parameters
    
    % find initial values for the montecarlo simulations
    pos = zeros(nPar,1);  % pre allocation
    for j =1:nPar
        pos(j) = find(strcmp(par2names(j), parameters.parAbb));
    end
    parameters.parValues(pos)= parSimulation;
    
    % update initial conditions 
    xInitial(1,pos_glu) = Xval1(i,1);
    xInitial(1,pos_lac) = Xval1(i,2);
    xInitial(1,pos_pro) = Xval1(i,3);
    xInitial(1,pos_ace) = Xval1(i,4);
    xInitial(1,pos_bm_su) = Xval1(i,5);
    xInitial(1,pos_bm_lac) = Xval1(i,6);
    xInitial(1,pos_ye) = Xval1(i,7);
    
    % run the ode solver
    [t, y] = f_run_co_culture_parameter_dependant(tSim,xInitial,parameters);
    
    % reservoir extraction
    reservoir = zeros(length(t),14);
    helpReservoir = zeros(1,14);
    volReservoir = 1; %L
    % loop to calculate reservoir
    for k = 2:length(t)
        states = y(k,:)';
        [~,~,transferReed]=f_mass_balances(t(k),states,parameters);
        helpReservoir = helpReservoir + transferReed'*(t(k)-t(k-1))*states(1)/volReservoir;
        reservoir(k,:) =helpReservoir;
    end
    
    % calculate cost function for each simulation 
    gluAdded  = sum(parameters.spikeFlowrate.volumeSpikes.* ... % in gCOD
        parameters.spikeFlowrate.compositionSpikes(:,1)) ...  % added by spike
        + y(end,pos_glu)*y(end,1) + ...     % glucose left over
        + y(1,pos_glu)*y(1,1);    %present in the begining of the batch
    
    propMade = reservoir(end,pos_pro);
    
    %calculate yield
    productivity(i) = -reservoir(end,pos_pro)/t(end);
    %yield(i) = propMade/gluAdded ; % YIELD gCOD/gCOD
    yieldActual = -propMade/gluAdded ; % YIELD gCOD/gCOD
    yield(i)  = max(yieldActual,-0.70);  % so simulation doesn't stretch out
    
    % save to matrix structure 
    matrix1_glu(:,i) = y(:,pos_glu);
    matrix1_lac(:,i) = y(:,pos_lac);
    matrix1_pro(:,i) = y(:,pos_pro);
    matrix1_ace(:,i) = y(:,pos_ace);
    matrix1_ye(:,i)  = y(:,pos_ye);
    matrix1_bm_su(:,i) = y(:,pos_bm_su);
    matrix1_bm_lac(:,i) = y(:,pos_bm_lac);
    matrix1_tan(:,i) = y(:,pos_tan);    
end

Matrix.glu = matrix1_glu;
Matrix.lac = matrix1_lac;
Matrix.pro = matrix1_pro;
Matrix.ace = matrix1_ace;
Matrix.ye = matrix1_ye;
Matrix.bm_su = matrix1_bm_su;
Matrix.bm_lac = matrix1_bm_lac;
Matrix.tan = matrix1_tan;



end

