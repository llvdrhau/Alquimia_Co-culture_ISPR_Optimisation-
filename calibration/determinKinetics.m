% %Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
% May 2020. Please contact lucas.vanderhauwaert@usc.es if you
% intend to use this code.
% runs different kinentic models of the veillionella spike experiment V3/4
% to dertermin the best fit and if the parameters of this simumaltion is
% applicable to the co-culture experiments 

clear 
clc
%% run Veillionella estimation ON PRO-LAC INHIBITION
kinetics = 'seperateProLac';
f_run_veillionella_spike_estimation(kinetics);


%% less effective inhibition models 
%% run Veillionella estimation ON TOTAL VFA inhibition 
kinetics = 'totalVFA';
f_run_veillionella_spike_estimation(kinetics);

%% run Veillionella estimation ON LAC inhibition 
kinetics = 'proOnly';
f_run_veillionella_spike_estimation(kinetics);

%% run Veillionella estimation ON PRO inhibition 
kinetics = 'lacOnly';
f_run_veillionella_spike_estimation(kinetics);





