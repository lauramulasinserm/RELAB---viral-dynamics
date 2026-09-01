# Monolix modelling

This directory contains the code used to prepare the RELAB datasets and fit the within-host viral dynamics models presented in the manuscript.

- SARS-CoV-2
- Influenza A virus (IAV)
- Influenza B virus (IBV)
- Respiratory syncytial virus (RSV)

Individual-level RELAB data are not publicly available and are therefore not included in this repository.

## Monolix models

Models were fitted using the SAEM algorithm in Monolix. The structural model, observation model, inter-individual variability, covariate effects, and parameter settings are specified in the corresponding `.mlxtran` files.

## Software

Analyses were performed using Monolix Suite 2023R1.
