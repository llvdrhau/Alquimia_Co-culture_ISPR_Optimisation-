function [ReedData] = get_reed_data(sheetname)
% Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% October 2021.Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.
%   Retrives data for spikes from excel file ReedData.xlsx 

[data,Names]=xlsread('ReedData.xlsx',sheetname);

tStart = data(:,strcmp(Names,'tStart'));
ReedData.tStart = tStart;

tEnd = data(:,strcmp(Names,'tEnd'));
ReedData.tEnd = tEnd;

RateMatrix = data(:,3:end);
ReedData.RateMatrix = RateMatrix;

end

