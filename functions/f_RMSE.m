function [RMSENor, r2] = f_RMSE(tExp,yExp,ysim,tsim,compounds)
% Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% October 2021.Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.
%   calculates RSME and R2 for evaluation of calibration/validation


%N = length(tExp);
%position of data 
if iscell(ysim)
    yPar = [];
    yExperimental = [];
    Nall = 0;
    for i = 1:length(ysim)
        N = length(tExp{i});
        Nall = Nall + N;
        matrixTime = repmat(tsim{i},1,length(tExp{i}));
        matrixExpTime = repmat(tExp{i}',length(tsim{i}),1);
        [~,pos] = min(abs(matrixExpTime-matrixTime));
        ySimulation = ysim{i};
        yExperimental = [yExperimental; yExp{i}]; 
        yPar = [yPar; ySimulation(pos,:)]; % select data at time t of experiments
    end 

    if size(yPar,2) > 5 % In the Co-Cultures the bacteria are a mix of the two
        yPar(:,5) = yPar(:,end-1) + yPar(:,end);
        yPar = yPar(:,1:5);
    end

    %% RMSE
    RMSE = sqrt(sum((yExperimental-yPar).^2)/N);
    nor = max(yExperimental)-min(yExperimental);
    RMSENor = RMSE./nor;

    %% Regresion coef
    yMean = mean(yExperimental);
    RSS = sqrt(sum((yExperimental-yPar).^2));
    TSS = sqrt(sum((yExperimental-yMean).^2));
    r2 = 1- RSS./TSS;


else
    N = length(tExp);
    matrixTime = repmat(tsim,1,length(tExp));
    matrixExpTime = repmat(tExp',length(tsim),1);
    [~,pos] = min(abs(matrixExpTime-matrixTime));
    yPar =  ysim(pos,:); % select data at time t of experiments

    if size(yPar,2) > 5 % In the Co-Cultures the bacteria are a mix of the two
        yPar(:,5) = yPar(:,end-1) + yPar(:,end);
        yPar = yPar(:,1:5);
    end
    %% RMSE
    RMSE = sqrt(sum((yExp-yPar).^2)/N);
    nor = max(yExp)-min(yExp);
    RMSENor = RMSE./nor;

    %% Regresion coef
    yMean = mean(yExp);
    RSS = sqrt(sum((yExp-yPar).^2));
    TSS = sqrt(sum((yExp-yMean).^2));
    r2 = 1- RSS./TSS;
end



%% print to terminal 
PR(1,:) = compounds;
PR(2,:) = num2cell(RMSENor);
disp('Normalised RMSE for each output is:')
disp(PR)
disp('')

PR2(1,:) = compounds;
PR2(2,:) = num2cell(r2);
disp('r^2 for each output is:')
disp(PR2)
disp('')

end

