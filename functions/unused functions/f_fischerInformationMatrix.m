function FIM = f_fischerInformationMatrix(parameters,parVal,parName,xInitial,control,sigma,sheetData)
%explanation 
%

%% Get number of parameters
par_names = parName;
par_value = parVal;
p = length(par_names);

%% prepare the Q-matrix
Qdiag = 1./sigma;%.^2;
Q = diag(Qdiag); % diag command

%% Get the sensitivity fo each parameter and put it in 3D-matrix
allsens = [];
for i=1:p
    sens = f_forwaardSensitivity(parameters,par_names{i}, par_value(i), xInitial,control,sheetData);
    allsens = cat(3, allsens, sens);                  
end


%% Calculate the FIM 
N = size(sens, 1);
FIM=zeros(p, p);
for i = 1:N
    A = squeeze(allsens(i, 2:end, :));
    FIM = FIM + A'*Q*A;
end

end

