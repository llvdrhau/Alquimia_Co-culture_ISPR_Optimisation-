function [] = plotSensitivity(y,idx,parameterReached,xInitial,parameters,tExp,par2names, stateNames)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

allsens = [];
par_names = stateNames;
par_value = par2names;
p = length(idx);

for i=1:p
    [~,sens] = f_sensitivity(y,idx(i),parameterReached,xInitial,parameters,tExp,par2names, stateNames);
    allsens = cat(3, allsens, sens);                  
end

%% plot sensitivities
states = parameters.compounds(2:end);
t = allsens(:,1,1);
poslegend = zeros(1,length(par_value));

for k = 1:4
    y = squeeze(allsens(:,k,:));
    subplot(2,2,k)
    hold on
    for i = 1:length(par_value)
        yState = y(:,i);
        plot(t,yState)
    end
end
legend(par_names)
title(stateNames{k})
%poslegend = zeros(1,length(par_value));
hold off
end



%% 
% StartRow = 1;
% endRow  = 22;
% time = 1:22;
% 
% for i = 1:length(stateNames)
%     subplot(1,length(stateNames),i)
% 
%     for j = 1:length(par2names)
%         dydp(StartRow:endRow,j)
%         plot(time, dydp(StartRow:endRow,j))
%           
%     end
%     StartRow = StartRow + 22;
%     endRow  = endRow + 22;
%     legend(par2names)
%     title(stateNames{i}) 
% end
