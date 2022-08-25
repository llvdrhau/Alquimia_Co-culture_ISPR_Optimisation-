function [Fspikes,inSpikes]=f_spikes_flowrate(t, parameters)
% Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% October 2021.Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.
%    Returns relevant data to simulate the addition of spikes to the reactor

nCompounds=parameters.nStates-1;              %Get the number of compounds

spikeFlowrate = parameters.spikeFlowrate;

tSpikes = spikeFlowrate.tSpikes;
volumeSpikes = spikeFlowrate.volumeSpikes;
compositionSpikes = spikeFlowrate.compositionSpikes;

nSpikes = length(tSpikes);
flowrate = zeros(1,nSpikes);
contributionSpike = zeros(nSpikes,nCompounds);
spikeHalfLength = 0.1; %h

for i=1:nSpikes
   tStart = tSpikes(i)- spikeHalfLength;
   tEnd = tSpikes(i)+ spikeHalfLength;
   switchSpike = -0.5*(1+tanh(100*(tStart-t)))+0.5*(1+tanh(100*(tEnd-t)));
   flowrate(i) = switchSpike*volumeSpikes(i)/(2*spikeHalfLength);
   contributionSpike(i,:) = flowrate(i)*compositionSpikes(i,:);
end

Fspikes = sum(flowrate);
%inSpikes = zeros(nCompounds,1);
totalSpikes = sum(contributionSpike)/(1e-9+Fspikes);
inSpikes = totalSpikes';


