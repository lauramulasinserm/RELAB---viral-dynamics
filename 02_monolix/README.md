# Monolix modelling

This directory contains the Monolix model files used to characterize within-host viral dynamics of the four respiratory viruses included in the RELAB study:

- SARS-CoV-2
- Influenza A virus (IAV)
- Influenza B virus (IBV)
- Respiratory syncytial virus (RSV)

Individual-level RELAB data are not publicly available and are therefore not included in this repository.

## Model structure

Viral kinetics were described using a target-cell-limited mechanistic model.

For each virus, two files are provided:

- `model_monolix_*.mlxtran`: Monolix project file specifying the statistical model, parameter estimation, inter-individual variability, covariate effects, observation model, and estimation settings.
- `modele_TV_*.txt`: structural model defining the within-host viral dynamics equations.

The corresponding files are:

| Virus | Monolix model | Structural model |
|---|---|---|
| SARS-CoV-2 | `model_monolix_COVID.mlxtran` | `modele_TV_COVID.txt` |
| Influenza A | `model_monolix_IAV.mlxtran` | `modele_TV_IAV.txt` |
| Influenza B | `model_monolix_IBV.mlxtran` | `modele_TV_IBV.txt` |
| RSV | `model_monolix_RSV.mlxtran` | `modele_TV_RSV.txt` |

## Model estimation

Models were fitted using the Stochastic Approximation Expectation-Maximization (SAEM) algorithm implemented in Monolix.

The `.mlxtran` files contain the final models used for the analyses presented in the manuscript, including:

- fixed-effect parameter estimates;
- inter-individual variability;
- residual error model;
- retained covariate effects;
- parameter constraints and fixed parameters;
- Fisher Information Matrix and likelihood settings;
- model diagnostic settings.

The Fisher Information Matrix and log-likelihood were computed using linearization after model convergence.

## Data

The input datasets required by the `.mlxtran` files are generated from the RELAB data-management workflow available in the `01_data_management` directory.

Because the RELAB dataset contains individual-level clinical and virological information, these input datasets are not distributed in this public repository.

## Software

Models were fitted using **Monolix Suite 2023R1**.
