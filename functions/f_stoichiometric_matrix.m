function stoM = f_stoichiometric_matrix(stoTable,parameters)
% Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% October 2021.Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.
%    Stoichiometric matrix evaluation  

switch parameters.kineticFlag
    case 'CoCulture'
        parAbb = parameters.parAbb;
        parValues = parameters.parValues;
        %Identify parameters
        Cxi = parValues(strcmp(parAbb,'Cxi'));
        Csu = parValues(strcmp(parAbb,'Csu'));
        Cye = parValues(strcmp(parAbb,'Cye'));
        Caa = parValues(strcmp(parAbb,'Caa'));
        Nye = parValues(strcmp(parAbb,'Nye'));
        Nbac = parValues(strcmp(parAbb,'Nbac'));
        Nxc = parValues(strcmp(parAbb,'Nxc'));
        Cpro = parValues(strcmp(parAbb,'Cpro'));
        Clac = parValues(strcmp(parAbb,'Clac'));
        Cxc = parValues(strcmp(parAbb,'Cxc'));
        Cac = parValues(strcmp(parAbb,'Cac'));
        Cbac = parValues(strcmp(parAbb,'Cbac'));
        Ysu = parValues(strcmp(parAbb,'Ysu'));
        Yyev = parValues(strcmp(parAbb,'Yyev'));
        Yyeb = parValues(strcmp(parAbb,'Yyeb'));
        Ylac = parValues(strcmp(parAbb,'Ylac'));
        
    otherwise
        parAbb = parameters.parAbb;
        parValues = parameters.parValues;
        %Identify parameters
        Cxi = parValues(strcmp(parAbb,'Cxi'));
        Csu = parValues(strcmp(parAbb,'Csu'));
        Cye = parValues(strcmp(parAbb,'Cye'));
        Caa = parValues(strcmp(parAbb,'Caa'));
        Nye = parValues(strcmp(parAbb,'Nye'));
        Nbac = parValues(strcmp(parAbb,'Nbac'));
        Nxc = parValues(strcmp(parAbb,'Nxc'));
        Cpro = parValues(strcmp(parAbb,'Cpro'));
        Clac = parValues(strcmp(parAbb,'Clac'));
        Cxc = parValues(strcmp(parAbb,'Cxc'));
        Cac = parValues(strcmp(parAbb,'Cac'));
        Cbac = parValues(strcmp(parAbb,'Cbac'));
        Ysu = parValues(strcmp(parAbb,'Ysu'));
        Yye = parValues(strcmp(parAbb,'Yye'));
        Ylac = parValues(strcmp(parAbb,'Ylac'));
end
stoMatrixCell=stoTable(2:end,2:end-1);

%Replace parameters by their value
[np,nc]=size(stoMatrixCell);
stoM=zeros(np,nc);
for i=1:np
    for j=1:nc
        elementCell=stoMatrixCell{i,j};
        if isa(elementCell,'double')
            stoM(i,j) = elementCell;
        else
            stoM(i,j)=eval(elementCell);
        end
    end
end
