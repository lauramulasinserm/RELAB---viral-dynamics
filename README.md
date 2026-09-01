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
├── 01_data_management/
│   ├── 01_data_cleaning.R
│   ├── 02_EEQ_normalization.R
│   ├── 03_data_preparation_for_monolix.R
│   ├── 04_characteristics_table.R
│   └── README.md
│
├── 02_monolix/
│   ├── model_monolix_COVID.mlxtran
│   ├── model_monolix_IAV.mlxtran
│   ├── model_monolix_IBV.mlxtran
│   ├── model_monolix_RSV.mlxtran
│   ├── modele_TV_COVID.txt
│   ├── modele_TV_IAV.txt
│   ├── modele_TV_IBV.txt
│   ├── modele_TV_RSV.txt
│   └── README.md
│
├── 03_results/
│   ├── analyses/
│   │   ├── covariate_effects/
│   │   ├── vaccine_effectiveness/
│   │   └── virus_comparisons/
│   │
│   ├── model_assessment/
│   │   ├── simulate_COVID_with_residual_error.R
│   │   ├── simulate_IAV_with_residual_error.R
│   │   ├── simulate_IBV_with_residual_error.R
│   │   ├── simulate_RSV_with_residual_error.R
│   │   └── README.md
│   │
│   ├── simulations/
│   │   ├── simulate_COVID_trajectories.R
│   │   ├── simulate_IAV_trajectories.R
│   │   ├── simulate_IBV_trajectories.R
│   │   ├── simulate_RSV_trajectories.R
│   │   └── README.md
│   │
│   └── README.md
│
├── .gitignore
└── README.md
```


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

For questions regarding the code or analyses, please contact the corresponding repository author.
