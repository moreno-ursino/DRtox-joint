# Implementation of the method to determine the maximal tolerated dose-regimen using PK/PD

This is the code for the paper "Bayesian modeling of a bivariate toxicity outcome for early phase oncology trials evaluating dose regimens" authored by Emma Gerard, Sarah Zohar, Christelle Lorenzato, Moreno Ursino and Marie-Karelle Riviere. 

## Prerequisite

PK/PD estimation is performed via `R` using the software Monolix (http://lixoft.com/products/monolix/) through the package `lixoftConnectors`.

Bayesian analysis is performed with Stan via the package `rstan` (https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).


## Implementation

### Tools folder

The tools folder contains all generic functions used for the simulations:
  - `pk_functions.R` contains all PK/PD functions and functions to generate PK/PD data  
  - `simu_trials_PKPD.R` contains the functions to simulate the data under the CRM design and estimate the PKPD models with Monolix
  - `simu_trials_stat.R` contains the functions to extract PK/PD estimates and fit the statistical models on each dataset to estimate the probability of toxicity at each dose regimen
  - `CRM_2_param` contains the functions perform the 2 parameter-CRM
  
### Scenario folder

The scenario folder contains programs to define the toxicity scenarios:
  1. `run_Rmax.R`: Inside the Rmax folder, simulation of the PD response of a large number of patients (using `PK_PD_param.R`). It generates the individual Rmax `Rmax_local.Rdata`
  2. `run_sc.R`: Definition of the toxicity scenarios from the data generated with `run_Rmax.R`. It generates the `scenario.Rdata` file in the scenario folder containing all infomation on the scenarios and in particular the true probabilities of DLT, CRS and DLTo.
  
### Run

We will develop the steps to run the simulations:

  1. Change the `path_simus` variable in `run_simu_PKPD.R` and `run_simu_stat.R` to the path of the folder `code`
  2. Define the PK/PD parameters in `PK_PD_param.R`
  3. Run the PK/PD estimation with `run_simu_PKPD.R`. For each trial, data and estimated PK/PD parameters are written in different folders.
  4. Run the statistical estimation with `run_simu_stat.R`. Results are written in a `.Rdata` file in the results folder of the design x scenario chosen.

The code on one dataset is available in the folder ex_one_trial.
  
## Disclaimer

Programs are provided "as is" without any warranty.


For any issue, please contact emma.g025@gmail.com
