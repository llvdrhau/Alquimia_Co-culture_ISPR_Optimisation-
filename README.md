# Co-culture fermentation model and optimisation 

Lucas Van der Hauwaert. University of Santiago de Compostela. Spain
October 2021.

Please contact lucas.vanderhauwaert@usc.es if you intend to use this code.

Code to develop a model of a co-culture reactor producing proprionate 
while implementing In-situ product removal (ISPR) by reverse enhanced 
electro dialisis   


folders: 
* excel files: contain the excel file containing experimental data, 
extraction data, spike data and information for the structure "parameters" 

* get_functions: functions to retrive data from excel files 

* functions: functions to run the various scripts throught all the folders

* kinetics: function holding the kinetic expresions used in the mass
 balances 

* Sampling: functions to perform a latin hyper cube sampling incorporating 
correlations between parameters according to the methods of 
Iman-Conover (2007). 

* Calibration: runs parameter estimation, Monte Carlo simuilations and the
bootstrap methode for the two batch models

* Validation: comibing parameters for both models and validating them on
 experiment BV5
 
* optimisation: scripts running single, multiple and stochastic 
optimisation 


!! 
Please run RUN_ME_FIRST to connect all the folders to the working directory
!! 

bib

Iman, Ronald L, and W J Conover. 2007. “Communications in Statistics - Simulation and A Distribution-Free Approach to Inducing Rank Correlation 
among Input Variables,” no. September 2012: 37–41.

Van der hauwaert et al.(2022)
DOI: XXXXXXXXX



