function [spikeFlowrate] = get_spike_data(sheetname)
%Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% May 2020. Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.
%   Retrives data for spikes from excel file spikeData.xlsx 


[data,Names]=xlsread('spikeData.xlsx',sheetname);

tSpikes=data(:,strcmp(Names,'time'));
spikeFlowrate.tSpikes = tSpikes;

spikeFlowrate.volumeSpikes = data(:,(strcmp(Names,'Vol')))/1000; % in L  
% voor moest het moelijk doen vol = data(:,logical(sum(strcmp(Names,'Vol'))))/1000
spikeFlowrate.compositionSpikes = data(1:end,3:end);
end

