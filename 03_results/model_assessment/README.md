# VPC-like model evaluation

R scripts used to evaluate model predictions by comparing observed and simulated Ct distributions.

Simulations include inter-individual variability and residual error and are compared with observed Ct values across time since symptom onset.

## Scripts

- `simulate_COVID_with_residual_error.R`: simulation of Ct trajectories including residual error
- `RELAB_simu_Ct_distrib.R`: comparison of observed and simulated Ct distributions for SARS-CoV-2, IAV, IBV, and RSV
