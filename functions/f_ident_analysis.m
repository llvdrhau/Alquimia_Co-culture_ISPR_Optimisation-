function [col,dmsqr_m,sNorm,Combi] = f_ident_analysis(dydp,dyNDdp,idx,par2names)

sAbs = dydp;              %Absolute sensitiviy
sND = dyNDdp;             %Non-dimensional sensitivity

for i=1:length(idx)
    sNorm(:,i) = sND(:,i)./norm(sND(:,i));          %Normalizamos por todas a sensibilidades dun par�metro
end

nsm     = sNorm'*sNorm;                                                    % Informa sobre a sensibilidade dos par�metros entre si (de a� que sexa sim�trica e diagonal 1).
asm     = sAbs'*sAbs;                                                      % Absolute sensitivity matrix, fim Fisher Information Matrix
dtm     = det(asm)^(1/(2*length(idx)));                                    % Determinant index (ind�canos se hai ou non combinaci�ns lineais na matriz de sensibilidades absolutas. Canto maior, mellor)
col     = 1/sqrt(min(eig(nsm)));                                           % Collinearity index
for i =1:length(idx) % for each parameter
    dmsqr(i) = sqrt((sND(:,i)'*sND(:,i))/length(sND(:,i)));                % Suma de sensibilidades dun par�metro. Se � baixo � un par�metro con baixa sensibilidade.
end

dmsqr_m(1,:) = par2names;
dmsqr_m(2,:) = num2cell(dmsqr);

[~,iRank]=sort(dmsqr,'descend');                                           % Ordeno os par�metros por orde de sensibilidade, de maior a menor.
Combi = struct();
for k=2:length(idx);    
    C = combnk(iRank,k);
    eval(horzcat('Combi.num.G',num2str(k),'=C;'));
    eval(horzcat('Combi.letras.L',num2str(k),'=par2names(C);'));
    [nCombinations,~] = size(C);
    for i=1:nCombinations
        iSelected = C(i,:);
        nsm     = sNorm(:,iSelected)'*sNorm(:,iSelected);                  % normalized sensitivity matrix for the selected parameters
        col(k,i)     = 1/sqrt(min(eig(nsm)));                              % collinearity index for the selected parameters
    end
end


end

