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
