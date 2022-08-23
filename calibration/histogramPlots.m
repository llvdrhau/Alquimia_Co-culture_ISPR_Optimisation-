
%% Veilionella plots
load('bootstrapResultsVeillionella')
nParB = size(parameterMatrix,2);
parNamesLx = {'km_{su}' 'km_{ye}' 'km_{lac}' 'KS_{su}' 'KS_{ye}' 'KS_{lac}','Y_{su}','Y_{ye}','Y_{lac}','KI_{VFA}' 'lag_time','K_{YE}', 'kdec_{Xlac}', 'KI_{pro}', 'KI_{lac}'};
parUnits = {'(h^{-1})','(gCOD_{X}/gCOD_{lac})', '(gCOD_{pro}/L^{-1})','(gCOD_{lac}/L^{-1})','(h^{-1})'};
idx = [3 9 14 15 13]; %Parameters to be determined from the list above
parLatex = parNamesLx(idx);

for i=1:nParB
    %figure
    subplot(3,3,i)
    set(gcf,'color','w'); %background color
    histogram(parameterMatrix(:,i))
    txtX = [parLatex{i} parUnits{i}];
    xlabel(txtX)
    ylabel('Frequency')
end

%% Blacilus plots
load('bootstrapResultsBacillus2')
nPar = size(parameterMatrix,2);
parNamesLx = {'km_{su}' 'km_{ye}' 'km_{lac}' 'KS_{su}' 'KS_{ye}' 'KS_{lac}','Y_{su}','Y_{ye}','Y_{lac}','KI_{VFA}' 'lag_time','K_{YE}', 'kdec_{Xlac}', 'KI_{pro}', 'KI_{lac}'};
parUnits = {'(h^{-1})','(gCOD_{X}/gCOD_{glu})', '(-)'};
idx = [1 7];
parLatex = parNamesLx(idx);

for j = 1:2
    subplot(3,3,5+j)
    set(gcf,'color','w'); %background color
    histogram(parameterMatrix(:,j))
    txtX = [parLatex{j} parUnits{j}];
    xlabel(txtX)
    ylabel('Frequency')
end


