function [ parameterReported, parameterMatrix, parameterReported_dev,CImatrix ] = ...
    f_bootstrap(yExp,par2estimate,lb,ub,options,tExp,xInitial,parameters,par2names,parameterReached,newBootstrap)
% Script to run the bootstrap method to create robust parameters 
% Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% October 2021.Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.

% input control
if nargin < 11
    newBootstrap = 1;         %Flag (if bootstrap already done, just load the data)
end
%% Bootstrap
caso = 'bootstrap';
switch parameters.flagCase
    case 'batchBacilus'
        [residual] = f_residual_bacillus(parameterReached,tExp,yExp,xInitial,parameters,par2names,caso);
    case 'SpikeVeillonella'
        [residual] = f_residual_V_spike(parameterReached,tExp,yExp,xInitial,parameters,par2names,caso);
end
[n,m] = size(yExp);
residual_m = reshape(residual,n,m);

[rows, columns] = size(yExp);
% now do a bootrap sampling from residuals and re-perform parameter estimation
nboot = 500;               % number of bootstrap iterations

parameterMatrix = zeros(nboot,length(parameterReached));
%optionsLsq = optimset('MaxFunEvals',100,'MaxIter',50);
if newBootstrap
    rng(1) % for reproduceability 
    residualSampled = zeros(rows,columns);
    for k=1:nboot
        disp(['the iteration number is : ',num2str(k)])
        
        for i=1:columns
            if any(isnan(yExp(:,i)))
                posi = find(isfinite(yExp(:,i)));
                iSampling = ceil(length(posi)*rand(length(posi),1));
                for j=1:length(posi)
                    residualSampled(posi(j),i) = residual_m(posi(iSampling(j)),i);
                end
            else
                iSampling = ceil(rows*rand(rows,1));
                    residualSampled(:,i) = residual_m(iSampling,i);
            end
        end
        % random sampling with replacement from measurement errors
        yExpIteration = yExp + residualSampled ;       % artificial data by randomly adding some of the reference residuals
        
        % Do a parameter estimation
        caso='normal';
        par2estimate_boots = par2estimate.*(0.9+0.2*rand(size(par2estimate)));
        
        switch parameters.flagCase
            case 'batchBacilus'
        [parameterIteration,~,~,~,~,~,~] = lsqnonlin(@f_residual_bacillus,par2estimate_boots,lb,ub,options,tExp,yExpIteration,xInitial,parameters,par2names,caso);
            case 'SpikeVeillonella'
        [parameterIteration,~,~,~,~,~,~] = lsqnonlin(@f_residual_V_spike,par2estimate_boots,lb,ub,options,tExp,yExpIteration,xInitial,parameters,par2names,caso);
        end

        parameterMatrix(k,:)=parameterIteration;                % store estimated parameters
    end
else
    switch parameters.flagCase
        case 'batchBacilus'
            load bootstrapResultsBacillus2.mat   %if previously run bootstrap, load results directly
        case 'SpikeVeillonella'
            load bootstrapResultsVeillionella.mat   %if previously run bootstrap, load results directly
    end
end

%Find parameter values from the distribution of parameters (pMatrix)
parameterReported=mean(parameterMatrix);
disp('The mean of distribution of parameters are')
parameterReported_m(1,:) = par2names;
parameterReported_m(2,:) = num2cell(parameterReported);
disp(parameterReported_m)
disp('The std.dev. of distribution of parameters are')
parameterReported_dev(1,:) = par2names;
parameterReported_dev(2,:) = num2cell(std(parameterMatrix));
disp(parameterReported_dev)
disp('The confidence interval (alpha = 0.05) of parameters are')
CImatrix(1,:) = par2names;
CImatrix(2:3,:) = num2cell(prctile(parameterMatrix,[2.5 97.5]));
disp(CImatrix)


end

