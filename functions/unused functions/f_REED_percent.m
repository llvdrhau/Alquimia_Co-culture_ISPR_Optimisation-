function [transfer] = f_REED_percent(states,parameters,t)
%Determins the transfer of mass to the accumulation reservoir 
%   Detailed explanation goes here

transfer = 0*states;

Reed = parameters.Reed; 
nReed = length(Reed.tStart); 

tStart = Reed.tStart;
tEnd = Reed.tEnd;
RateMatrix = Reed.RateMatrix;
ExtractionRate = zeros(nReed,length(RateMatrix));

for i=1:nReed
   switchReed = -0.5*(1+tanh(100*(tStart(i)-t)))+0.5*(1+tanh(100*(tEnd(i)-t)));
   %ExtractionRate(i,:) = switchReed*states(2:end)'.*RateMatrix(i,:).*states(1)/(tEnd(i) - tStart(i));
   ExtractionRate(i,:) = switchReed*states(2:end)'.*RateMatrix(i,:)/(tEnd(i) - tStart(i));
end

rTransfer = sum(ExtractionRate); 


compoundAbb=parameters.compoundAbb;
transfer(strcmp(compoundAbb,'Sglu'))=rTransfer(1);           % Glucose
transfer(strcmp(compoundAbb,'Slac'))=rTransfer(2);           % Lactate
transfer(strcmp(compoundAbb,'Sye'))=rTransfer(3);            % Yeast extract
transfer(strcmp(compoundAbb,'Spro'))=rTransfer(4);           % Propionate
transfer(strcmp(compoundAbb,'Sac'))=rTransfer(5);            % Acetate
transfer(strcmp(compoundAbb,'Sh2'))=rTransfer(6);            % H2
transfer(strcmp(compoundAbb,'Stic'))=rTransfer(7);           % TIC
transfer(strcmp(compoundAbb,'SIn'))=rTransfer(8);            % TAN
transfer(strcmp(compoundAbb,'Xc'))=rTransfer(9);             % Composites
transfer(strcmp(compoundAbb,'Xsu'))=rTransfer(10);           % Sugar degraders
transfer(strcmp(compoundAbb,'Xlac'))=rTransfer(11);          % Lactate degraders
transfer(strcmp(compoundAbb,'Scat'))=rTransfer(11);          % Cations 
transfer(strcmp(compoundAbb,'San'))=rTransfer(11);           % Anions 


end
