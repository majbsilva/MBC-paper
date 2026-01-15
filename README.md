
# MBC paper

<!-- badges: start -->
<!-- badges: end -->

# Competing Risk Analysis in Metastatic Breast Cancer

This repository contains the R scripts used to perform competing risk analyses
for the study submitted to *Clinical Breast Cancer*. The analytical pipeline
implements Fine & Gray subdistribution hazard models, cause-specific Cox models,
model discrimination using time-dependent concordance index (C-index),
calibration analysis, and nomogram construction.

## Data Sources and Availability
The analyses were conducted using secondary data obtained from:

- Surveillance, Epidemiology, and End Results (SEER) Program (USA)
- São Paulo Oncocenter Foundation (FOSP) Hospital-Based Cancer Registry (Brazil)

Due to data use agreements and registry policies, the original datasets cannot
be redistributed. Researchers interested in reproducing the analyses must obtain
access directly from the data providers:

- SEER data access: https://seer.cancer.gov/data/access.html  
- FOSP public interface: https://fosp.saude.sp.gov.br/fosp/direto-ao-dado/

After access approval, datasets must be harmonized according to the variable
definitions described in the manuscript and saved as `mbc.csv`.


## Reproducibility
The script is fully executable once the harmonized dataset is used. All analyses are deterministic and reproducible using the fixed random seed specified in the code.

## Software Requirements
- R version ≥ 4.3.0
- Required packages: cmprsk, riskRegression, pec, rms, survival, dplyr, readr

## License
This repository is distributed under the MIT License.

