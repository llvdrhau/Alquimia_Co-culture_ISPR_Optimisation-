# Co-culture Fermentation Model and Optimization

**Lucas Van der Hauwaert, University of Santiago de Compostela, Spain**  
*October 2021*

For inquiries or to request the use of this code, please contact: [lucas.vanderhauwaert@usc.es](mailto:lucas.vanderhauwaert@usc.es)

## Description
This code develops a model for a co-culture reactor designed for proprionate production, incorporating In-situ Product Removal (ISPR) through reverse enhanced electro-dialysis.

## Contents

1. **Excel Files**: Contains data files for experiments, extraction, and spikes, along with information for the "parameters" structure.
2. **Get_Functions**: Functions to retrieve data from Excel files.
3. **Functions**: Essential functions for running scripts across all folders.
4. **Kinetics**: Contains kinetic expressions used in mass balances.
5. **Sampling**: Implements Latin Hypercube Sampling with correlations based on Iman-Conover's method (2007).
6. **Calibration**: Executes parameter estimation, Monte Carlo simulations, and the bootstrap method for batch models.
7. **Validation**: Combines parameters for both models, validating them against experiment BV5.
8. **Optimisation**: Includes scripts for single, multiple, and stochastic optimization.

## Important Note
Please run *RUN_ME_FIRST* to connect all folders to the working directory.

## Bibliography
- Iman, Ronald L, and W J Conover. (2007). "A Distribution-Free Approach to Inducing Rank Correlation among Input Variables." Communications in Statistics - Simulation and Computation, no. September 2012: 37–41.

## Publication 
- Van der Hauwaert, L., Regueira, A., Selder, L., Zeng, A. P., & Mauricio-Iglesias, M. (2022). "Optimising bioreactor processes with in-situ product removal using mathematical programming: A case study for propionate production." Computers & Chemical Engineering, 168, 108059.  [DOI](https://doi.org/10.1016/j.compchemeng.2022.108059)

## Acknowledgments
This work was supported by the ALQUIMIA project (PID2019-110993RJ-I00), funded by the Agencia Estatal de Investigación under the Programa Retos de la sociedad, modalidad Jovenes investigadores, convocatoria 2019.


