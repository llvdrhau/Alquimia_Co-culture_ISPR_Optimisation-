function [transfer] = f_REED(states,parameters,t)
% Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% October 2021.Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.
%   Determins the transfer of mass to the accumulation reservoir due to
%   reverse enhanced electrodialisis at time t 


transfer = 0*states;

Reed = parameters.Reed; 
nReed = length(Reed.tStart); % amount of reed processes 

tStart = Reed.tStart;
tEnd = Reed.tEnd;
RateMatrix = Reed.RateMatrix;
ExtractionRate = zeros(nReed,length(RateMatrix));


minConc = ones(1,length(RateMatrix))*0.1;
switchConc = (tanh(30*(states(2:end)'-minConc))+1)/2;
% to make the extraction zero if their is nothing to extract 
% 0,1 is the threshold: this is to ensure that extraction is not
% taken into acount when ploting the reed reservoir 


for i=1:nReed % loop over to see if t falls in the activation period 
   switchReed = -0.5*(1+tanh(100*(tStart(i)-t)))+0.5*(1+tanh(100*(tEnd(i)-t)));
   ExtractionRate(i,:) = switchReed*switchConc.*RateMatrix(i,:); 
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
