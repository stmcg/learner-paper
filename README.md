# LEARNER: A Transfer Learning Method for Low-Rank Matrix Estimation

This repository contains the code for the simulation study and data application in the manuscript ["LEARNER: A Transfer Learning Method for Low-Rank Matrix Estimation"](
https://doi.org/10.48550/arXiv.2412.20605) by Sean McGrath, Cenhao Zhu, Min Guo, and Rui Duan.

---

## Simulation Study

The folder `simulations` contains the code for the simulation study. Below is a description of the contents of this folder.

### Prerequisites

The code requires the following R packages to be installed: 

- **doParallel** (version 1.0.17)
- **doRNG** (version 1.8.6)
- **foreach** (version 1.5.2)
- **learner** (version 0.2.0)
- **MASS** (version 7.3-60.2)
- **ScreeNOT** (version 0.1.0)
- **RColorBrewer** (version 1.1-3)

The package version numbers listed above were used in the analyses in the manuscript.

### File Structure

#### 1. Helper Files

* ``helper.R``: Helper functions for running all simulations
*  ``helper_ext.R``: Additional helper function for running the simulations with an external dataset for LEARNER

#### 2. Running Simulations

##### Independent Noise Scenarios

* ``runsim-same.R``: Runs simulation scenarios with high similarity in the latent spaces
* ``runsim-dif-moderate.R``: Runs simulation scenarios with moderate similarity in the latent spaces
* ``runsim-dif_verylow.R``: Runs simulation scenarios with low similarity in the latent spaces

##### Correlated Noise Scenarios

* Files ending with ``-cor`` (e.g., ``runsim-same-cor.R``): Runs simulation scenarios with correlation of 0.1
* Files ending with ``-cor2``: Runs simulation scenarios with correlation of 0.25
* Files ending with ``-cor3``: Runs simulation scenarios with correlation of 0.5

##### Correlated Noise + External Dataset Scenarios

* Files ending with ``-cor-ext`` (e.g., ``runsim-same-cor-ext.R``): Runs simulation scenarios with an external dataset

The `results` folder contains the output from running these simulations. Note that these files typically require between 8 to 30 hours to run when parallelized across 10 CPU cores.

#### 3. Analyzing Results

* ``analyze_results.R``: Creates the figures summarizing the simulation results


---

## Data Application

The folder `application` contains the code for the simulation study. Below is a description of the contents of this folder.

### Prerequisites

The code for the data processing requires the following R packages to be installed: 

- **data.table** (version 1.16.0)  
- **genio** (version 1.1.2)  

The code for the data analysis requires the following R packages to be installed:

- **doParallel** (version 1.0.17)  
- **foreach** (version 1.5.2)  
- **lattice** (version 0.22-6)  
- **learner** (version 0.2.0)  
- **MASS** (version 7.3-60.2)  
- **RColorBrewer** (version 1.1-3)  
- **ScreeNOT** (version 0.1.0)  
- **softImpute** (version 1.4-1)  

The package version numbers listed above were used in the analyses in the manuscript.

### File Structure

#### Main Folder

- `key.RData`: Contains a mapping of the column names of `Y_0` and `Y_1` to phenotype names and ICD-10 codes.  
- `ScreeNOT.R`: Applies ScreeNOT to `Y_0` and `Y_1`.  
- `analyze_latent_spaces.R`: Performs analyses of the latent spaces based on the source-only, target-only, LEARNER, and D-LEARNER estimates. The analyses based on the LEARNER and D-LEARNER estimates require running the analyses in the subfolder `apply_learner_dlearner`.

#### Subfolder: `data_processing`

The folder contains the code for the data processing. The source data can be downloaded from the [NBDC Human Database](https://humandbs.dbcls.jp/en/hum0197-v3-220). After downloading the source data, the data processing steps should be run in the following order: 

1. `matrix_generation.R`: Performs variant screening and other data processing steps to help construct the matrices `Y_0` and `Y_1`.
2. `eur_filter.py` and `filter_new_eur.py`: Computes z-scores from the beta coefficients and standard errors in the European population based on the output of `matrix_generation.R`.
3.  `bbj_filter.py` and `filter_new_bbj.py`: Computes z-scores from the beta coefficients and standard errors in the BioBank Japan population based on the output of `matrix_generation.R`.
4.  `create_analytic_datasets.R`: Performs final data processing tasks (e.g., row and column naming and re-ordering) to create the analytic data sets `Y_0` and `Y_1`.

Note that step 1 requires approximately one week to run.

#### Subfolder: `apply_learner_dlearner`

This folder corresponds to the application of LEARNER and D-LEARNER in Section 4.2 of the manuscript.

- `dlearner.R`: Applies D-LEARNER. Results are saved in `DLEARNER-estimate-BBJ.RData`.  
- `find-lambda.R`: Selects hyperparameters for LEARNER by cross-validation. Results are saved in `lambda_small.RData`.  
- `learner.R`: Applies LEARNER with the selected hyperparameters. Results are saved in `learner.RData`.  
- `analyze-results.R`: Plots the results of the hyperparameter selection analyses and plots the convergence of LEARNER.  

#### Subfolder: `cross-validation`

This folder corresponds to the cross-validation analyses in Section 4.3 of the manuscript.

- `dlearner-svd.R`: Applies D-LEARNER and the target-only SVD method. D-LEARNER results are saved in `dlearner_setX.RData` and target-only SVD results are saved in `hard_setX.RData` for training set `X` (`X=1,2,3,4,5`).
- `find-lambda-helper.R`: Helper file for selecting hyperparameters for LEARNER.  
- `find-lambda-setX.R`: Selects hyperparameters for LEARNER in training set `X` (`X=1,2,3,4,5`). Results are saved in `lambda_small_setX.RData`.  
- `learner.R`: Applies LEARNER with the selected hyperparameters in each test dataset. Results are saved in `learner_setX.RData` for training set `X` (`X=1,2,3,4,5`).  
- `heldout-mse.R`: Computes the mean squared error values in the test datasets.  

The selection of hyperparameters for LEARNER took between 3 and 4 hours to run when parallelized across 13 CPU cores. The application of LEARNER took up to 3 hours to run on a single CPU core.

