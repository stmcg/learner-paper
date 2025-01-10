# LEARNER: A Transfer Learning Method for Low-Rank Matrix Estimation

This repository contains the code for the simulation study and data application in the manuscript ["LEARNER: A Transfer Learning Method for Low-Rank Matrix Estimation"](
https://doi.org/10.48550/arXiv.2412.20605) by Sean McGrath, Cenhao Zhu, Min Guo, and Rui Duan.

---

## Simulations

The folder `Simulations` contains the code for the simulation study. Below is a description of the contents of this folder.

### Pre-requisites

The code requires the following R packages to be installed: 

- **doParallel** (version 1.0.17)
- **doRNG** (version 1.8.6)
- **foreach** (version 1.5.2)
- **learner** (version 0.1.0)
- **MASS** (version 7.3-60.2)
- **ScreeNOT** (version 0.1.0)
- **RColorBrewer** (version 1.1-3)

The package version numbers listed above were used in the analyses in the manuscript.

### File structure

#### 1. Helper Files

* ``helper.R``: Helper functions for running all simulations
*  ``helper_ext.R``: Additional helper function for running the simulations with an external dataset for LEARNER

#### 2. Running Simulations

##### Independent Noise Scenarios

* ``main-same.R``: Runs simulation scenarios with high similarity in the latent spaces
* ``main-dif-moderate.R``: Runs simulation scenarios with moderate similarity in the latent spaces
* ``main-dif_verylow.R``: Runs simulation scenarios with low similarity in the latent spaces

##### Correlated Noise Scenarios

* Files ending with ``-cor`` (e.g., ``main-same-cor.R``): Runs simulation scenarios with correlation of 0.1
* Files ending with ``-cor2``: Runs simulation scenarios with correlation of 0.25
* Files ending with ``-cor3``: Runs simulation scenarios with correlation of 0.5

##### Correlated Noise + External Dataset Scenarios

* Files ending with ``-cor-ext`` (e.g., ``main-same-cor-ext.R``): Runs simulation scenarios with an external dataset

The `results` folder contains the output from running these simulations. Note that these files typically require between 8 to 30 hours to run when parallelized across 10 CPU cores.

#### 3. Analyzing results

* ``analyze_results.R``: Creates the figures summarizing the simulation results


---

## Data Application

The folder `Application` contains the code for the simulation study. Below is a description of the contents of this folder.

### Pre-requisites

The code requires the following R packages to be installed:

- **doParallel** (version 1.0.17)  
- **foreach** (version 1.5.2)  
- **lattice** (version 0.22-6)  
- **learner** (version 0.1.0)  
- **MASS** (version 7.3-60.2)  
- **RColorBrewer** (version 1.1-3)  
- **ScreeNOT** (version 0.1.0)  
- **softImpute** (version 1.4-1)  

The package version numbers listed above were used in the analyses in the manuscript.

### File Structure

