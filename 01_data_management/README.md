# Data management

This directory contains the R scripts used to clean, harmonize, and prepare the RELAB data for the analyses presented in the manuscript.

Because RELAB contains individual-level clinical and virological data, the raw and processed datasets are not included in this public repository.

## Workflow

The scripts should be run in the following order:

### 1. `01_data_cleaning.R`

Cleans and harmonizes the raw RELAB datasets from the 2024/25 and 2025/26 epidemic seasons.

Main steps include:

- removal of duplicate observations;
- harmonization of dates, age, sex, and PCR results;
- quality control of Ct values;
- harmonization of symptom information;
- vaccination-status processing;
- definition of infection episodes;
- reconstruction of time since symptom onset;
- harmonization of influenza typing and PCR technique names;
- combination of the two epidemic seasons.

Output:

`data/processed/RELAB_24_26_clean_start.csv`

### 2. `02_EEQ_normalization.R`

Uses External Quality Evaluation (EEQ) data to estimate PCR technique- and virus-specific Ct normalization rules.

For each virus and PCR technique, the mean Ct deviation from the overall EEQ mean is calculated across high, medium, and low viral-load samples.

These correction factors are then applied to RELAB Ct values.

Outputs:

- `data/processed/EEQ_Ct_normalization_rules.csv`
- normalized RELAB dataset used for downstream analyses.

### 3. `03_data_preparation_for_monolix.R`

Creates the virus-specific datasets used for within-host viral dynamics modelling in Monolix.

Separate datasets are generated for:

- SARS-CoV-2;
- influenza A virus (IAV);
- influenza B virus (IBV);
- respiratory syncytial virus (RSV).

The script defines the observation variable, censoring status, epidemic period, season, vaccination status, and other covariates required by the Monolix models.

### 4. `04_characteristics_table.R`

Creates the descriptive table of the study population by virus.

The table includes:

- epidemic season;
- sex;
- age category;
- vaccination status;
- symptoms;
- co-infections;
- laboratory;
- number of tests per infection;
- presence of a negative follow-up test.

It also generates vaccination coverage summaries by virus and epidemic season.

## Data availability

Individual-level RELAB data are not publicly available because they contain potentially sensitive clinical and virological information.

The code provided here documents the complete data-processing workflow used in the study.

## Software

Data management and descriptive analyses were performed in R.
