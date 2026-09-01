# Within-host viral dynamics of respiratory infections using RELAB data

This repository contains the code associated with the manuscript:

**"Using routine community testing data to quantify viral dynamics in symptomatic SARS-CoV-2, influenza, and RSV infections, France, 2024–2026"**

## Overview

This study uses community-based RT-PCR data from the French RELAB network to characterize the within-host viral dynamics of four major respiratory viruses:

- SARS-CoV-2
- Influenza A virus (IAV)
- Influenza B virus (IBV)
- Respiratory syncytial virus (RSV)

Within-host viral kinetics were characterized using a mechanistic target-cell-limited model fitted to routine cycle threshold (Ct) data.

The repository contains the code used for data preparation, Monolix model fitting, simulation of viral trajectories, analysis of viral kinetic characteristics, and generation of the manuscript figures and tables.


## Repository structure

```text
RELAB---viral-dynamics/
│
├── 01_data_processing/
│
├── 02_monolix/
│   ├── COVID/
│   ├── IAV/
│   ├── IBV/
│   └── RSV/
│
├── 03_results_analysis/
│   │
│   ├── 01_simulate_COVID_trajectories.R
│   ├── 02_simulate_IAV_trajectories.R
│   ├── 03_simulate_IBV_trajectories.R
│   ├── 04_simulate_RSV_trajectories.R
│   │
│   ├── model_outputs/
│   │   ├── COVID/
│   │   ├── IAV/
│   │   ├── IBV/
│   │   └── RSV/
│   │
│   └── results/
│       └── simulations/
│           ├── COVID/
│           ├── IAV/
│           ├── IBV/
│           └── RSV/
│
└── README.md
```


## Directory description

### `01_data_processing/`

Scripts used to prepare and format the RELAB data for the viral kinetic analyses.

Individual-level RELAB data are not included in the repository.


### `02_monolix/`

Monolix model files used to estimate the population parameters of the within-host viral dynamics models.

Models are organized separately for:

- SARS-CoV-2 (`COVID`)
- Influenza A virus (`IAV`)
- Influenza B virus (`IBV`)
- Respiratory syncytial virus (`RSV`)


### `03_results_analysis/`

R scripts used to simulate viral trajectories from the estimated Monolix model parameters and derive the results presented in the manuscript.

The main trajectory simulation scripts are organized by virus:

- `01_simulate_COVID_trajectories.R`
- `02_simulate_IAV_trajectories.R`
- `03_simulate_IBV_trajectories.R`
- `04_simulate_RSV_trajectories.R`

The `model_outputs/` directory contains selected Monolix outputs required for downstream analyses and simulations.

The `results/simulations/` directory is used to store simulation outputs and is organized separately for each virus.

Large simulation files are not tracked in the public repository.


## Data availability

Individual-level RELAB data are not publicly available in this repository.

Please refer to the manuscript for further information regarding data access and availability.


## Software

Analyses were performed using:

- **R**
- **Monolix**

Viral kinetic models were fitted using Monolix, while simulations, statistical analyses, data processing, and figure generation were performed in R.


## Reproducibility

The analysis workflow broadly consists of:

1. Preparation of RELAB RT-PCR data.
2. Estimation of within-host viral kinetic parameters using Monolix.
3. Simulation of viral trajectories using the estimated population parameters.
4. Calculation of viral kinetic characteristics and covariate effects.
5. Generation of manuscript figures and tables.

Individual-level RELAB data and large simulation files are not included in the repository.


## Citation

If you use this code, please cite the associated manuscript:

**Mulas L, et al. "Using routine community testing data to quantify viral dynamics in symptomatic SARS-CoV-2, influenza, and RSV infections, France, 2024–2026."**

Citation details will be updated upon publication.


## Contact

**Laura Mulas**

Université Paris Cité / Inserm  
IAME – Infection, Antimicrobials, Modelling, Evolution  

For questions regarding the code or analyses, please contact the corresponding repository author.
