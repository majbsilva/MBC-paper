# ==============================================================================
# TITLE: FOSP Data Preprocessing and Cleaning Pipeline
# DESCRIPTION: 
#   This script processes hospital-based cancer registry data from FOSP (São Paulo).
#   It merges raw DBF files, filters for breast cancer (C50), standardizes 
#   variables (staging, histology, molecular subtypes), calculates time-to-event 
#   intervals (TTS, TCS, TTR, Survival), and computes geographical distances.
#
# INPUT: 
#   - Raw .dbf files in 'data/raw/'
#   - Municipality codes .ods file in 'data/auxiliary/'
# OUTPUT: 
#   - Cleaned CSV file in 'output/'
# ==============================================================================

# 1. SETUP & LIBRARIES ---------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readr, foreign, janitor, dplyr, tidyverse, stringr, readODS, geobr, sf, geosphere)

# Set paths (Adjust if necessary)
PATH_RAW <- "data/raw"
PATH_AUX <- "data/auxiliary"
PATH_OUT <- "output"

# Ensure output directory exists
if(!dir.exists(PATH_OUT)) dir.create(PATH_OUT)

# 2. DATA INGESTION (DBF MERGE) ------------------------------------------------
message(">>> Loading and merging DBF files...")

dbf_files <- list.files(path = PATH_RAW, pattern = "\\.dbf$", full.names = TRUE, ignore.case = TRUE)

if(length(dbf_files) == 0) stop("No .dbf files found in ", PATH_RAW)

fosp_list <- lapply(dbf_files, read.dbf, as.is = TRUE)
fosp_raw  <- do.call(rbind, fosp_list)

message(">>> Total Raw Records: ", nrow(fosp_raw))

# Optional: Save raw combined checkpoint
# write.csv(fosp_raw, file.path(PATH_OUT, "fosp_combined_raw.csv"), row.names = FALSE)

# 3. FILTERING ELIGIBLE CASES (BREAST CANCER) ----------------------------------
message(">>> Filtering eligible C50 cases...")

# Inclusion Criteria:
# 1. Topography: C50 (Breast)
# 2. Age >= 18
# 3. Diagnosis Basis: 3 (Microscopic confirmation)
# 4. Treatment started (DTTRAT is not NA)
# 5. Invasive tumors only (Stage I-IV) - Filtered later but prepared here

fosp_c50 <- fosp_raw %>%
  filter(TOPOGRUP == "C50") %>%
  mutate(patient_id = row_number()) # Create ID before further filtering

message(">>> Total C50 Cases: ", nrow(fosp_c50))

# 4. DATA CLEANING & TRANSFORMATION --------------------------------------------
message(">>> Processing variables and dates...")

fosp_clean <- fosp_c50 %>%
  # --- Date Formatting ---
  mutate(across(c(DTCONSULT, DTDIAG, DTTRAT, DTULTINFO, DTRECIDIVA), as.Date)) %>%
  
  # --- Basic Filters (Age, Treatment, Microscopic Conf) ---
  filter(
    !is.na(DTTRAT),
    IDADE >= 18,
    BASEDIAG == "3"
  ) %>%
  
  # --- Demographics ---
  mutate(
    SEXO = recode(SEXO, '1' = "Male", '2' = "Female"),
    education_level = case_when(
      ESCOLARI == 1 ~ "Illiterate",
      ESCOLARI == 2 ~ "Primary Incomplete",
      ESCOLARI == 3 ~ "Primary Complete",
      ESCOLARI == 4 ~ "High School",
      ESCOLARI == 5 ~ "Higher Education",
      ESCOLARI == 9 ~ "Unknown",
      TRUE ~ NA_character_
    )
  ) %>%
  
  # --- Time Intervals (Delays & Survival) ---
  mutate(
    # TTS: Time to Start Treatment (Diagnosis to Treatment)
    TTS = as.numeric(difftime(DTTRAT, DTDIAG, units = "days")),
    TTS_status = case_when(
      TTS < 0 ~ "Inconsistent",
      TTS == 0 ~ "Same-day diagnosis and treatment",
      TTS <= 60 ~ "Recommended window",
      TTS > 60 ~ "Delayed treatment",
      TRUE ~ NA_character_
    ),
    TTS_CAT60 = ifelse(TTS >= 60, ">= 60 days", "< 60 days"),
    
    # TCS: Time form Consultation to Treatment
    TCS = as.numeric(difftime(DTTRAT, DTCONSULT, units = "days")),
    TCS_status = case_when(
      TCS < 0 ~ "Inconsistent",
      TCS == 0 ~ "Same-day treatment",
      TCS <= 60 ~ "Recommended window",
      TCS > 60 ~ "Delayed treatment",
      TRUE ~ NA_character_
    ),
    TCS_CAT60 = ifelse(TCS >= 60, ">= 60 days", "< 60 days"),
    
    # Survival
    SURV_DAYS = as.numeric(DTULTINFO - DTDIAG),
    survival_months = round(SURV_DAYS / 30.44, 3)
  ) %>%
  filter(SURV_DAYS >= 0) # Remove inconsistent survival times

# 5. CLINICAL & PATHOLOGICAL VARIABLES -----------------------------------------
message(">>> Extracting clinical features (Stage, Histology, Biomarkers)...")

# Define Valid Malignant Histology Codes (AJCC)
valid_histo_codes <- c(
  "80223", "80323", "80353", "80413", "80703", "82003", "82013", "82113",
  "82463", "8290/", "83143", "83153", "84103", "84303", "84803", "85003",
  "85023", "85033", "85043", "85073", "85093", "85103", "85133", "85203",
  "85253", "85303", "85403", "85703", "85713", "85723", "85743", "85753",
  "89823", "89833", "80003", "80103", "81403", "82553", "84013", "85013",
  "85213", "85223", "85233", "85243", "85413", "85433"
)

fosp_clinical <- fosp_clean %>%
  # Filter Valid Histology & Invasive Stages
  filter(CIDO %in% valid_histo_codes) %>%
  filter(ECGRUP %in% c("I", "II", "III", "IV")) %>%
  
  # --- Histological Description ---
  mutate(
    HISTDESC = case_when(
      CIDO == "85003" ~ "Invasive ductal carcinoma, NOS",
      CIDO == "85203" ~ "Invasive lobular carcinoma",
      CIDO == "85223" ~ "Infiltrating duct and lobular carcinoma",
      CIDO == "80003" ~ "Neoplasm, malignant",
      CIDO == "80103" ~ "Carcinoma, NOS",
      TRUE ~ "Other/Unclassified" # Simplified for brevity, add full list if needed
    )
  ) %>%
  
  # --- Regex Extraction for Grade & Biomarkers (from 'OUTRACLA') ---
  mutate(
    # Histological Grade
    GRAUHIST_raw = str_match(OUTRACLA, regex("GRAU\\s*([123]|I{1,3})\\b", ignore_case = TRUE))[,2],
    histologic_grade = case_when(
      GRAUHIST_raw %in% c("1", "I") ~ 1L,
      GRAUHIST_raw %in% c("2", "II") ~ 2L,
      GRAUHIST_raw %in% c("3", "III") ~ 3L,
      TRUE ~ NA_integer_
    ),
    
    # ER Status
    .er_clean = str_replace_all(toupper(OUTRACLA), "[- ]", ""),
    er_status = case_when(
      str_detect(.er_clean, "ERPOS|ESTROGPOS") ~ "Positive",
      str_detect(.er_clean, "ERNEG|ESTROGENEG") ~ "Negative",
      TRUE ~ NA_character_
    ),
    
    # PR Status
    .pr_clean = str_replace_all(toupper(OUTRACLA), "[- ]", ""),
    pr_status = case_when(
      str_detect(.pr_clean, "PROG(?:EST)?POS") ~ "Positive",
      str_detect(.pr_clean, "PROG(?:EST)?NEG") ~ "Negative",
      TRUE ~ NA_character_
    ),
    
    # HER2 Status
    .her2_clean = str_replace_all(toupper(OUTRACLA), "[- ]", ""),
    her2_status = case_when(
      str_detect(.her2_clean, "HER2POS") ~ "Positive",
      str_detect(.her2_clean, "HER2NEG") ~ "Negative",
      TRUE ~ NA_character_
    )
  ) %>%
  select(-starts_with("."), -GRAUHIST_raw) # Cleanup helper columns

# 6. TREATMENT & OUTCOMES ------------------------------------------------------
message(">>> Standardizing treatments and outcomes...")

fosp_final_vars <- fosp_clinical %>%
  filter(NENHUMANT == 1) %>% # No prior treatment outside hospital
  mutate(
    # Event Indicator (0=Censor, 1=Cancer Death, 2=Other Death)
    event_indicator2 = case_when(
      ULTINFO %in% c(1, 2) ~ 0L, # Alive
      ULTINFO == 3 ~ 1L,         # Dead (Cancer)
      ULTINFO == 4 ~ 2L,         # Dead (Other)
      TRUE ~ NA_integer_
    ),
    
    # Recurrence
    TTR = as.numeric(difftime(DTRECIDIVA, DTTRAT, units = "days")),
    TTR_months = round(TTR / 30.44, 2),
    recurrence_status = ifelse(!is.na(TTR), "Recurrence", "No recurrence"),
    
    # Treatment Modalities
    surgery_status = ifelse(CIRURGIA == 1, "YES", "NO"),
    chemotherapy_status = ifelse(QUIMIO == 1, "YES", "NO"),
    radiotherapy_status = ifelse(RADIO == 1, "YES", "NO"),
    hormone_therapy_status = ifelse(HORMONIO == 1, "YES", "NO"),
    
    # Laterality
    laterality = recode_factor(as.character(LATERALI), 
                               `1`="right", `2`="left", `3`="bilateral", `8`="not applicable")
  )

# 7. GEOSPATIAL ANALYSIS (DISTANCES) -------------------------------------------
message(">>> Performing geospatial analysis (Municipality distances)...")

# Load Auxiliary Municipality Data (ODS)
# Ensure this file is in 'data/auxiliary/'
ods_file <- list.files(PATH_AUX, pattern = "RELATORIO_DTB_.*\\.ods", full.names = TRUE)[1]

if(!is.na(ods_file)) {
  muni_names <- read_ods(ods_file, sheet = 1, skip = 6) %>% 
    select(`Código Município Completo`, Nome_Município) %>%
    mutate(IBGE_CODE = as.character(`Código Município Completo`))
  
  # Join names
  fosp_geo <- fosp_final_vars %>%
    mutate(IBGE = as.character(IBGE), IBGEATEN = as.character(IBGEATEN)) %>%
    left_join(muni_names, by = c("IBGE" = "IBGE_CODE")) %>%
    rename(municipio_residencia = Nome_Município) %>%
    left_join(muni_names, by = c("IBGEATEN" = "IBGE_CODE")) %>%
    rename(municipio_tratamento = Nome_Município) %>%
    filter(!is.na(municipio_residencia))
  
  # Calculate Coordinates & Distance using geobr (Requires Internet)
  # Wrap in tryCatch to allow offline execution
  tryCatch({
    message("   Downloading municipality coordinates via geobr...")
    muni_coords <- read_municipality(year = 2020, showProgress = FALSE) %>%
      st_centroid() %>%
      mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2]) %>%
      st_drop_geometry() %>%
      mutate(code_muni = as.character(code_muni)) %>%
      select(code_muni, lon, lat)
    
    # Finalização da lógica de geolocalização e cálculo de distância
    fosp_geo <- fosp_geo %>%
      left_join(muni_coords, by = c("IBGE" = "code_muni")) %>%
      rename(lon_resid = lon, lat_resid = lat) %>%
      left_join(muni_coords, by = c("IBGEATEN" = "code_muni")) %>%
      rename(lon_trat = lon, lat_trat = lat) %>%
      mutate(
        distancia_km = distHaversine(
          cbind(lon_resid, lat_resid),
          cbind(lon_trat, lat_trat)
        ) / 1000
      )
  }, error = function(e) {
    message("   Error in geolocation: ", e$message)
    fosp_geo$distancia_km <- NA
  })
} 