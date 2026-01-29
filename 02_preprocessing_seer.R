# ==============================================================================
# TITLE: SEER Data Preprocessing and Cleaning Pipeline (Male Breast Cancer)
# DESCRIPTION: 
#   This script processes data from the Surveillance, Epidemiology, and End 
#   Results (SEER) program. It filters for male breast cancer, harmonizes 
#   staging systems across different eras (AJCC 6th, 7th, EOD 2018), standardizes 
#   histology and treatment variables, and performs data cleaning for survival analysis.
#
# INPUT: 
#   - Raw CSV file 'big_seer.csv' in 'data/raw/'
# OUTPUT: 
#   - Cleaned CSV file in 'output/'
# ==============================================================================

# 1. SETUP & LIBRARIES ---------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, janitor, stringr, stringi, readr)

# Set paths
PATH_RAW <- "data/raw"
PATH_OUT <- "output"

# Ensure output directory exists
if(!dir.exists(PATH_OUT)) dir.create(PATH_OUT)

# 2. DATA INGESTION ------------------------------------------------------------
message(">>> Loading SEER dataset...")

# Adjust filename if necessary
seer_file <- file.path(PATH_RAW, "big_seer.csv")

if(!file.exists(seer_file)) stop("File 'big_seer.csv' not found in ", PATH_RAW)

big_seer <- read_csv(seer_file, show_col_types = FALSE)

# 3. INITIAL FILTERING ---------------------------------------------------------
message(">>> Applying inclusion criteria...")

# Criteria: Breast Site, Positive Histology, Primary Tumor, Male
elegiveis <- big_seer %>%
  filter(
    `Site recode ICD-O-3/WHO 2008` == "Breast",
    `Diagnostic Confirmation` == "Positive histology",
    `Primary by international rules` == "Yes",
    `Sex` == "Male"
  ) %>%
  janitor::clean_names()

message(">>> Eligible cases after initial filter: ", nrow(elegiveis))

# 4. VARIABLE SELECTION --------------------------------------------------------
vars_interesse <- c(
  # Demographics
  "patient_id", "age_recode_with_single_ages_and_90", "race_recode_white_black_other",
  "marital_status_at_diagnosis", "median_household_income_inflation_adj_to_2022",
  "rural_urban_continuum_code", "race_and_origin_recode_nhw_nhb_nhaian_nhapi_hispanic_no_total",
  "type_of_reporting_source",
  
  # Staging (Multi-era)
  "year_of_diagnosis",
  "breast_adjusted_ajcc_6th_stage_1988_2015", "derived_seer_cmb_stg_grp_2016_2017",
  "derived_eod_2018_stage_group_2018", "x7th_edition_stage_group_recode_2016_2017",
  
  # TNM Components
  "breast_adjusted_ajcc_6th_t_1988_2015", "breast_adjusted_ajcc_6th_n_1988_2015", "breast_adjusted_ajcc_6th_m_1988_2015",
  "derived_eod_2018_t_2018", "derived_eod_2018_n_2018", "derived_eod_2018_m_2018",
  
  # Histology & Grade
  "diagnostic_confirmation", "histologic_type_icd_o_3", "icd_o_3_hist_behav_malignant",
  "breast_subtype_2010", "grade_clinical_2018", "grade_pathological_2018", "grade_recode_thru_2017",
  
  # Biomarkers
  "er_status_recode_breast_cancer_1990", "estrogen_receptor_summary_2018",
  "pr_status_recode_breast_cancer_1990", "progesterone_receptor_summary_2018",
  "her2_overall_summary_2018", "derived_her2_recode_2010",
  
  # Tumor Characteristics
  "laterality", "tumor_size_over_time_recode_1988", "regional_nodes_positive_1988",
  "record_number_recode", "site_recode_icd_o_3_who_2008", "primary_by_international_rules",
  
  # Treatment
  "surgery_of_oth_reg_dis_sites_1998_2002", "rx_summ_surg_prim_site_1998",
  "rx_summ_reg_ln_examined_1998_2002", "rx_summ_scope_reg_ln_sur_2003",
  "rx_summ_surg_oth_reg_dis_2003", "rx_summ_surg_rad_seq",
  "rx_summ_systemic_sur_seq_2007", "reason_no_cancer_directed_surgery",
  "radiation_recode", "chemotherapy_recode_yes_no_unk",
  "response_to_neoadjuvant_therapy_recode_2010",
  
  # Temporal & Outcomes
  "year_of_death_recode", "year_of_follow_up_recode", "survival_months",
  "survival_months_flag", "time_from_diagnosis_to_treatment_in_days_recode",
  "seer_cause_specific_death_classification", "vital_status_recode_study_cutoff_used",
  "seer_other_cause_of_death_classification"
)

elegiveis <- elegiveis %>% select(any_of(vars_interesse))

# 5. DATA CLEANING & TRANSFORMATION --------------------------------------------
message(">>> Processing clinical variables...")

# Valid Histology Codes (AJCC)
codigos_validos <- c(
  "8022", "8032", "8035", "8041", "8070", "8200", "8201", "8211",
  "8246", "8290", "8314", "8315", "8410", "8430", "8480", "8500",
  "8502", "8503", "8504", "8507", "8509", "8510", "8513", "8520",
  "8525", "8530", "8540", "8570", "8571", "8572", "8574", "8575",
  "8982", "8983", "8000", "8010", "8140", "8255", "8401", "8501",
  "8521", "8522", "8523", "8524", "8541", "8543"
)

elegiveis_clean <- elegiveis %>%
  # --- Age Cleaning ---
  mutate(
    age = str_remove(age_recode_with_single_ages_and_90, " years"),
    age = if_else(age == "90+", "90", age),
    age = as.numeric(age)
  ) %>%
  
  # --- Histology Filter ---
  mutate(
    codigo_histologico = str_extract(icd_o_3_hist_behav_malignant, "^\\d{4}"),
    histologia_valida = if_else(codigo_histologico %in% codigos_validos, "Sim", "Não"),
    histologic_type = str_remove(icd_o_3_hist_behav_malignant, "^\\d{4}/\\d:\\s*")
  ) %>%
  filter(histologia_valida == "Sim") %>%
  
  # --- Survival Months ---
  mutate(
    survival_months = str_remove_all(survival_months, "^0+"),
    survival_months = if_else(survival_months == "", "0", survival_months),
    survival_months = as.numeric(survival_months)
  ) %>%
  filter(survival_months > 0) %>% # Minimum 1 month follow-up
  
  # --- Race/Ethnicity (Custom Logic) ---
  mutate(
    race_ethnicity = case_when(
      race_and_origin_recode_nhw_nhb_nhaian_nhapi_hispanic_no_total == "Non-Hispanic White" ~ "White",
      race_and_origin_recode_nhw_nhb_nhaian_nhapi_hispanic_no_total == "Non-Hispanic Black" ~ "Black",
      race_and_origin_recode_nhw_nhb_nhaian_nhapi_hispanic_no_total == "Non-Hispanic Asian or Pacific Islander" ~ "Asian/Pacific Islander",
      race_and_origin_recode_nhw_nhb_nhaian_nhapi_hispanic_no_total == "Hispanic (All Races)" ~ "Hispanic",
      
      # Fallback to general race code
      is.na(race_and_origin_recode_nhw_nhb_nhaian_nhapi_hispanic_no_total) & race_recode_white_black_other == "White" ~ "White",
      is.na(race_and_origin_recode_nhw_nhb_nhaian_nhapi_hispanic_no_total) & race_recode_white_black_other == "Black" ~ "Black",
      TRUE ~ "Other/Unknown"
    )
  ) %>%
  rename(residence_region = rural_urban_continuum_code)

# 6. TREATMENT VARIABLES -------------------------------------------------------
message(">>> Standardizing treatment data...")

elegiveis_treat <- elegiveis_clean %>%
  # --- Surgery Status ---
  mutate(
    surgery_recode_desc = case_when(
      rx_summ_surg_prim_site_1998 == "00" ~ "No surgery",
      rx_summ_surg_prim_site_1998 %in% c(as.character(10:19)) ~ "Tumor destruction",
      rx_summ_surg_prim_site_1998 %in% c(as.character(20:80)) ~ "Resection",
      rx_summ_surg_prim_site_1998 == "90" ~ "Surgery, NOS",
      rx_summ_surg_prim_site_1998 == "99" ~ "Unknown",
      TRUE ~ NA_character_
    ),
    surgery_status = case_when(
      # Confirmed by code
      surgery_recode_desc %in% c("Resection", "Surgery, NOS") ~ "YES",
      # Inferred by sequence
      str_detect(rx_summ_systemic_sur_seq_2007, "surgery") | str_detect(rx_summ_surg_rad_seq, "surgery") ~ "YES",
      # Confirmed absence
      surgery_recode_desc == "No surgery" ~ "NO",
      TRUE ~ NA_character_
    )
  ) %>%
  # Filter out death certificate only / autopsy based on surgery field
  filter(!rx_summ_surg_prim_site_1998 %in% c("99", "00") | type_of_reporting_source != "Autopsy Only") %>%
  
  # --- Radiotherapy Status ---
  mutate(
    radiotherapy_status = case_when(
      str_detect(rx_summ_surg_rad_seq, "Radiation") ~ "YES",
      radiation_recode %in% c("Beam radiation", "Radioisotopes", "Radiation, NOS") ~ "YES",
      radiation_recode == "None/Unknown" ~ "NO",
      TRUE ~ NA_character_
    )
  ) %>%
  
  # --- Chemotherapy Status ---
  mutate(
    chemotherapy_status = case_when(
      chemotherapy_recode_yes_no_unk == "Yes" ~ "YES",
      chemotherapy_recode_yes_no_unk == "No/Unknown" ~ "NO",
      TRUE ~ NA_character_
    )
  )

# 7. STAGING AND TUMOR CHARACTERISTICS -----------------------------------------
message(">>> Harmonizing Tumor Stage (Cross-Era)...")

elegiveis_stage <- elegiveis_treat %>%
  mutate(
    # Harmonized Stage (AJCC 6th, 7th, EOD 2018)
    tumor_stage = case_when(
      # AJCC 6th (<= 2015)
      year_of_diagnosis <= 2015 & breast_adjusted_ajcc_6th_stage_1988_2015 == "0" ~ "0",
      year_of_diagnosis <= 2015 & breast_adjusted_ajcc_6th_stage_1988_2015 == "I" ~ "I",
      year_of_diagnosis <= 2015 & breast_adjusted_ajcc_6th_stage_1988_2015 %in% c("IIA", "IIB") ~ "II",
      year_of_diagnosis <= 2015 & str_detect(breast_adjusted_ajcc_6th_stage_1988_2015, "III") ~ "III",
      year_of_diagnosis <= 2015 & breast_adjusted_ajcc_6th_stage_1988_2015 == "IV" ~ "IV",
      
      # SEER Combined (2016-2017)
      year_of_diagnosis %in% 2016:2017 & derived_seer_cmb_stg_grp_2016_2017 == "0" ~ "0",
      year_of_diagnosis %in% 2016:2017 & str_detect(derived_seer_cmb_stg_grp_2016_2017, "^1") ~ "I",
      year_of_diagnosis %in% 2016:2017 & str_detect(derived_seer_cmb_stg_grp_2016_2017, "^2") ~ "II",
      year_of_diagnosis %in% 2016:2017 & str_detect(derived_seer_cmb_stg_grp_2016_2017, "^3") ~ "III",
      year_of_diagnosis %in% 2016:2017 & str_detect(derived_seer_cmb_stg_grp_2016_2017, "^4") ~ "IV",
      
      # EOD 2018 (>= 2018)
      year_of_diagnosis >= 2018 & derived_eod_2018_stage_group_2018 == "0" ~ "0",
      year_of_diagnosis >= 2018 & str_detect(derived_eod_2018_stage_group_2018, "^1") ~ "I",
      year_of_diagnosis >= 2018 & str_detect(derived_eod_2018_stage_group_2018, "^2") ~ "II",
      year_of_diagnosis >= 2018 & str_detect(derived_eod_2018_stage_group_2018, "^3") ~ "III",
      year_of_diagnosis >= 2018 & derived_eod_2018_stage_group_2018 == "4" ~ "IV",
      
      TRUE ~ NA_character_
    )
  ) %>%
  # Filter out In Situ ("0") and Unknown
  filter(tumor_stage %in% c("I", "II", "III", "IV")) %>%
  
  # --- Histologic Grade (Harmonized) ---
  mutate(
    histologic_grade = case_when(
      year_of_diagnosis <= 2017 & str_detect(grade_recode_thru_2017, "Grade I") ~ "1",
      year_of_diagnosis <= 2017 & str_detect(grade_recode_thru_2017, "Grade II") ~ "2",
      year_of_diagnosis <= 2017 & str_detect(grade_recode_thru_2017, "Grade III") ~ "3",
      year_of_diagnosis >= 2018 & grade_pathological_2018 %in% c("1", "2", "3") ~ grade_pathological_2018,
      TRUE ~ NA_character_
    )
  ) %>%
  
  # --- Biomarkers ---
  mutate(
    er_status = ifelse(er_status_recode_breast_cancer_1990 %in% c("Positive", "Negative"), er_status_recode_breast_cancer_1990, NA),
    pr_status = ifelse(pr_status_recode_breast_cancer_1990 %in% c("Positive", "Negative"), pr_status_recode_breast_cancer_1990, NA),
    her2_status = ifelse(derived_her2_recode_2010 %in% c("Positive", "Negative"), derived_her2_recode_2010, NA)
  ) %>%
  
  # --- TNM Extraction (Simplified) ---
  mutate(
    T_raw = coalesce(breast_adjusted_ajcc_6th_t_1988_2015, derived_eod_2018_t_2018),
    N_raw = coalesce(breast_adjusted_ajcc_6th_n_1988_2015, derived_eod_2018_n_2018),
    M_raw = coalesce(breast_adjusted_ajcc_6th_m_1988_2015, derived_eod_2018_m_2018),
    
    T_final = str_extract(T_raw, "^T[0-4]"),
    N_final = str_extract(N_raw, "^N[0-3]"),
    M_final = str_extract(M_raw, "^M[0-1]")
  )

# 8. OUTCOMES AND FINAL FORMATTING ---------------------------------------------
message(">>> Finalizing outcomes and formatting...")

seer_final <- elegiveis_stage %>%
  mutate(
    # Time to Treatment
    TTS = as.numeric(str_extract(time_from_diagnosis_to_treatment_in_days_recode, "\\d+")),
    TTS_CAT60 = ifelse(TTS >= 60, ">= 60 days", "< 60 days"),
    
    # Event Indicator (1 = Dead, 0 = Alive)
    event_indicator = if_else(vital_status_recode_study_cutoff_used == "Dead", 1, 0),
    
    # Tumor Size (Cleanup 990+)
    tumor_size_code = as.numeric(str_extract(tumor_size_over_time_recode_1988, "^\\d+")),
    tumor_size = case_when(
      tumor_size_code > 989 ~ NA_real_,
      TRUE ~ tumor_size_code
    ),
    
    # Lymph Nodes
    Lymph_nodes_positive = case_when(
      str_detect(regional_nodes_positive_1988, "^9[5-9]") ~ NA_real_,
      TRUE ~ as.numeric(regional_nodes_positive_1988)
    )
  ) %>%
  # Select final columns
  select(
    patient_id, year_of_diagnosis, age, race_ethnicity,
    marital_status_at_diagnosis, residence_region,
    histologic_type, histologic_grade,
    er_status, pr_status, her2_status,
    tumor_stage, T_final, N_final, M_final,
    tumor_size, Lymph_nodes_positive,
    surgery_status, chemotherapy_status, radiotherapy_status,
    laterality, TTS, TTS_CAT60,
    survival_months, event_indicator
  ) %>%
  # String Cleaning (lowercase, trim)
  mutate(across(where(is.character), ~ str_to_lower(str_squish(.x)))) %>%
  mutate(across(where(is.character), ~ stri_trans_general(.x, "Latin-ASCII")))

# 9. EXPORT --------------------------------------------------------------------
file_name <- paste0("seer_consolidated_male_", Sys.Date(), ".csv")
write.csv(seer_final, file.path(PATH_OUT, file_name), row.names = FALSE)

message(">>> SEER Processing complete! File saved: ", file.path(PATH_OUT, file_name))