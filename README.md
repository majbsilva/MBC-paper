
# Survival Analysis of Male Breast Cancer

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
# Competing Risk Analysis and Nomogram Construction

This repository contains the source code for the study **"Survival Analysis of Male Breast Cancer Accounting for Competing Risk: a Multicohort Study Integrating SEER and Hospital-based Data"**, submitted to *Clinical Breast Cancer*.

The repository is organized into a reproducible pipeline that performs:
1.  **Data Preprocessing:** Cleaning and harmonization of raw data from Brazilian (FOSP) and American (SEER) registries.
2.  **Statistical Analysis:** Implementation of Fine & Gray subdistribution hazard models, cause-specific Cox models, and model validation.
3.  **Visualization:** Generation of calibration plots and a prognostic nomogram.

## 📊 Data Sources and Access

This study utilizes secondary data. Due to data use agreements, **we cannot redistribute the raw datasets**. Researchers must obtain access directly from the providers:

1.  **SEER (Surveillance, Epidemiology, and End Results):**
    * Request access at: [SEER Data Access](https://seer.cancer.gov/data/access.html)
    * Download the case files and save as `big_seer.csv`.

2.  **FOSP (Fundação Oncocentro de São Paulo):**
    * Public access at: [FOSP Dados Públicos](https://fosp.saude.sp.gov.br/fosp/direto-ao-dado/)
    * Download the Hospital Cancer Registry `.dbf` files.
    * *Note:* The FOSP script also requires the IBGE municipality codes (`RELATORIO_DTB_BRASIL_2024_MUNICIPIOS.ods`).

## 🚀 Usage Instructions

The analysis pipeline is divided into three steps. Run the scripts in the following order:

### Step 1: FOSP Data Preparation
Run `01_preprocessing_fosp.R`.
* **Input:** Raw `.dbf` files and Municipality ODS.
* **Process:** Merges files, filters for breast cancer (C50), cleans clinical variables, and calculates geographical distances using `geobr`.
* **Output:** `fosp_consolidated.csv`

### Step 2: SEER Data Preparation
Run `02_preprocessing_seer.R`.
* **Input:** `big_seer.csv`.
* **Process:** Filters for male patients, harmonizes staging systems (AJCC 6th, 7th, and EOD 2018), and standardizes histology codes.
* **Output:** `seer_consolidated_male.csv`

### Step 3: Competing Risk Analysis & Nomogram
Run `03_analysis_mbc.R`.
* **Input:** The cleaned CSV files from Steps 1 and 2.
* **Process:**
    * Univariate and Multivariate Fine & Gray Models.
    * Comparison with Cause-Specific Cox Regression.
    * Model discrimination (C-index) and Calibration Plots (3, 5, and 8 years).
    * Construction of the Prognostic Nomogram.

## 📦 Software Requirements

The code was developed using **R version ≥ 4.3.0**. Please ensure the following packages are installed:

**Data Manipulation & Geolocation:**
```r
install.packages(c("tidyverse", "janitor", "readODS", "foreign", "geobr", "sf", "geosphere", "stringr", "stringi"))


## Reproducibility
The script is fully executable once the harmonized dataset is used. All analyses are deterministic and reproducible using the fixed random seed specified in the code.

## Software Requirements
- R version ≥ 4.3.0
- Required packages: cmprsk, riskRegression, pec, rms, survival, dplyr, readr

## License
This repository is distributed under the MIT License.

