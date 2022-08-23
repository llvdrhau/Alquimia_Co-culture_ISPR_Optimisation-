function [Matrix1, Matrix2 ] = f_montecarlo(par2estimate,parameterMatrix,xInitial,parameters,tExp,par2names,meanPar,sdPar,t)
%% Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% October 2021.Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.
%   runs a montecarlo simulation of the choosen parameters with LHS 
%   (iman an conover (2007))
%   INPUT: par2estimate = uncertain parameters. parameterMatrix =
%   distribution parameters!

if nargin < 9
    switch parameters.flagCase
        case 'batchBacilus'
             t=linspace(0,tExp(end),50);
             t = union(t,tExp);
        case 'SpikeVeillonella'
            t=linspace(0,tExp(end),50);
            t = union(t,tExp);
    end
    
end

caso= 'montecarlo';
nIt = 500;
nPar = length(par2estimate);

%% sampling parameters with inman conover, see file Sampling for connecting scripts
Xpar=lhs_empirco(parameterMatrix,nIt);

%% sampling initial conditions   
% assume that the initial conditions are uniformly distributed

nVal = 7;
dimensionlessX = lhsdesign(nIt,nVal);   
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

Xval2 = repmat(lbVal(2,:),nIt,1);
Xdiff2 = repmat(diffVal(2,:),nIt,1);
Xval2 = Xval2+Xdiff2.*dimensionlessX;

%% run 200 simulations for both data sets (output Matrix 1 and 2)
%preallocation 
matrix1_glu = zeros(length(t),nIt);
matrix1_lac = zeros(length(t),nIt);
matrix1_pro = zeros(length(t),nIt);
matrix1_ace = zeros(length(t),nIt);
matrix1_ye = zeros(length(t),nIt);
matrix1_bm_su = zeros(length(t),nIt);
matrix1_bm_lac = zeros(length(t),nIt);
matrix1_tan = zeros(length(t),nIt);

matrix2_glu = zeros(length(t),nIt);
matrix2_lac = zeros(length(t),nIt);
matrix2_pro = zeros(length(t),nIt);
matrix2_ace = zeros(length(t),nIt);
matrix2_ye = zeros(length(t),nIt);
matrix2_bm_su = zeros(length(t),nIt);
matrix2_bm_lac = zeros(length(t),nIt);
matrix2_tan = zeros(length(t),nIt);

for i = 1:nIt
   % select parameters  and initial conditions form sampling
    parSimulation = Xpar(i,:);
   
    xInitial(1,pos_glu) = Xval1(i,1);
    xInitial(1,pos_lac) = Xval1(i,2);
    xInitial(1,pos_pro) = Xval1(i,3);
    xInitial(1,pos_ace) = Xval1(i,4);
    xInitial(1,pos_bm_su) = Xval1(i,5);
    xInitial(1,pos_bm_lac) = Xval1(i,6);
    xInitial(1,pos_ye) = Xval1(i,7);
    
    xInitial(2,pos_glu) = Xval2(i,1);
    xInitial(2,pos_lac) = Xval2(i,2);
    xInitial(2,pos_ace) = Xval2(i,3);
    xInitial(2,pos_pro) = Xval2(i,4);
    xInitial(2,pos_bm_su) = Xval2(i,5);
    xInitial(2,pos_bm_lac) = Xval2(i,6);
    xInitial(2,pos_ye) = Xval2(i,7);
    % run simulation 
    switch parameters.flagCase
        case 'batchBacilus'
            [~,~,y1,y2] = f_solution_bacillus(xInitial,parameters,tExp,parSimulation,par2names,caso);
        case 'SpikeVeillonella'
            [~,y1] = f_solution(xInitial(1,:),parameters,tExp,parSimulation,par2names,caso);
            [~,y2] = f_solution(xInitial(2,:),parameters,tExp,parSimulation,par2names,caso);
    end
    % save output 
    matrix1_glu(:,i) = y1(:,pos_glu);
    matrix1_lac(:,i) = y1(:,pos_lac);
    matrix1_pro(:,i) = y1(:,pos_pro);
    matrix1_ace(:,i) = y1(:,pos_ace);
    matrix1_ye(:,i)  = y1(:,pos_ye);
    matrix1_bm_su(:,i) = y1(:,pos_bm_su);
    matrix1_bm_lac(:,i) = y1(:,pos_bm_lac);
    matrix1_tan(:,i) = y1(:,pos_tan);

    matrix2_glu(:,i) = y2(:,pos_glu);
    matrix2_lac(:,i) = y2(:,pos_lac);
    matrix2_pro(:,i) = y2(:,pos_pro);
    matrix2_ace(:,i) = y2(:,pos_ace);
    matrix2_ye(:,i) = y2(:,pos_ye);
    matrix2_bm_su(:,i) = y2(:,pos_bm_su);
    matrix2_bm_lac(:,i) = y2(:,pos_bm_lac);
    matrix2_tan(:,i) = y2(:,pos_tan);

    
end
%save all interested elemnts to structure 
Matrix1.glu = matrix1_glu;
Matrix1.lac = matrix1_lac;
Matrix1.pro = matrix1_pro;
Matrix1.ace = matrix1_ace;
Matrix1.ye = matrix1_ye;
Matrix1.bm_su = matrix1_bm_su;
Matrix1.bm_lac = matrix1_bm_lac;
Matrix1.tan = matrix1_tan;

Matrix2.glu = matrix2_glu;
Matrix2.lac = matrix2_lac;
Matrix2.pro = matrix2_pro;
Matrix2.ace = matrix2_ace;
Matrix2.ye = matrix2_ye;
Matrix2.bm_su = matrix2_bm_su;
Matrix2.bm_lac = matrix2_bm_lac;
Matrix2.tan = matrix2_tan;

end

