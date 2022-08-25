% Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% October 2021.Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.

% this script is to check out the growth rate of v. creseti in function of 
% the inhibition 
saveName = 'paretosearchOpt2_theOne'; % load saved data
load(saveName,'parameters','xInitial')

x3 = [169.2741    4.5953   18.5137    2.0500    2.7480    2.0500 ...
    1.9471    2.0500    1.7875]; % decesion variables that reult in strategy 3 

[~,parametersOut,t,y] = multiObjectiveFunc(x3,parameters,xInitial);
%% recall outputs to caculate inhibition and growth rate 

inhibition = zeros(length(t),1); 
growthrate = zeros(length(t),1); 

for i = 2:length(t)    
    states = y(i,:)';
    [~,~,~,Ivfa,growthRateLac]=f_mass_balances(t(i),states,parametersOut); 
    inhibition(i)= Ivfa;
    growthrate(i) = growthRateLac;
end
yyaxis left
plot(t,growthrate)
yyaxis right
plot(t,inhibition)
xlabel('time')
legend('growthrate','inhibition')
