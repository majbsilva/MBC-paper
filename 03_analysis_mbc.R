# ==============================================================================
# TITLE: Survival Analysis of Male Breast Cancer - Competing Risks & Nomogram
#
# DESCRIPTION:
# This script performs a survival analysis using the Fine-Gray subdistribution
# hazard model to account for competing risks. It compares the results with a
# cause-specific Cox model, validates the model (C-index, Calibration), and
# generates a prognostic nomogram.
#
# METHODOLOGY:
# 1. Fine-Gray Model (cmprsk) for subdistribution hazard ratios.
# 2. Model Validation using time-dependent C-index and calibration plots.
# 3. Nomogram construction using the weighted Cox regression approach (rms).
#
# REQUIRED DATA STRUCTURE (Input CSV):
# - survival_days: Time to event or censoring (numeric).
# - event_indicator: Standard censoring (0=Censored, 1=Event).
# - event_indicator2: Competing risks (0=Censored, 1=Event of Interest, 2=Competing Event).
# - Covariates: age_group, residence_region, tumor_stage, etc.
#
# DEPENDENCIES: cmprsk, riskRegression, rms, pROC, dplyr, ggplot2, pec
# ==============================================================================

# 1. SETUP AND PACKAGE LOADING -------------------------------------------------
required_packages <- c("cmprsk", "riskRegression", "rms", "pROC", 
                       "ggplot2", "gridExtra", "dplyr", "pec", "readr")

# Check and load packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
lapply(required_packages, library, character.only = TRUE)

# 2. DATA LOADING AND PREPROCESSING --------------------------------------------
# Load dataset (Replace 'data/dataset.csv' with your actual file path)
# Ensure the file follows the structure defined in the header.
mbc_raw <- read_csv('data/mbc_seer.csv', show_col_types = FALSE)

# Select variables of interest
vars_included <- c("survival_days", "event_indicator", "event_indicator2", 
                   "age_group", "residence_region", "histologic_grade",
                   "er_status", "pr_status", "her2_status", "tumor_stage",
                   "chemotherapy_status", "radiotherapy_status", 
                   "surgery_status", "TTS_CAT60")

# Create complete case dataset
mbc_complete <- mbc_raw[, vars_included]
mbc_complete <- na.omit(mbc_complete)

cat("Original observations:", nrow(mbc_raw), "\n")
cat("Complete cases included:", nrow(mbc_complete), "\n")
cat("Missing observations excluded:", nrow(mbc_raw) - nrow(mbc_complete), "\n")

# 3. DATA SPLITTING (TRAIN/TEST) -----------------------------------------------
set.seed(123) # Fixed seed for reproducibility
train_idx <- sample(1:nrow(mbc_complete), size = 0.7 * nrow(mbc_complete))

train_data <- mbc_complete[train_idx, ]
test_data  <- mbc_complete[-train_idx, ]

cat("\nTraining Set:", nrow(train_data), "- Events of Interest:", sum(train_data$event_indicator2 == 1), "\n")
cat("Testing Set:", nrow(test_data), "- Events of Interest:", sum(test_data$event_indicator2 == 1), "\n")

# 4. UNIVARIATE FINE-GRAY ANALYSIS ---------------------------------------------
# Analyze each covariate individually to screen for inclusion (p < 0.20 usually)

covariates_list <- c('age_group', 'residence_region', 'histologic_grade',
                     'er_status', 'pr_status', 'her2_status', 'tumor_stage',
                     'chemotherapy_status', 'radiotherapy_status', 
                     'surgery_status', 'TTS_CAT60')

cat("\n========== UNIVARIATE FINE-GRAY RESULTS ==========\n")

# Loop through covariates for cleaner output
for (var in covariates_list) {
  # Create model matrix for the specific variable
  formula_uni <- as.formula(paste("~", var))
  X_uni <- model.matrix(formula_uni, data = train_data)[, -1, drop = FALSE]
  
  fg_uni <- crr(
    ftime = train_data$survival_days,
    fstatus = train_data$event_indicator2,
    cov1 = X_uni,
    failcode = 1, # Event of interest (Cancer Death)
    cencode = 0   # Censored
  )
  
  cat("\nVariable:", var, "\n")
  print(summary(fg_uni))
}

# 5. MULTIVARIATE FINE-GRAY MODEL ----------------------------------------------
# Define the model matrix for selected variables (based on univariate screening)

# Prepare Train Matrix
X_train <- model.matrix(
  ~ age_group + 
    chemotherapy_status + 
    radiotherapy_status + 
    tumor_stage +
    her2_status +
    histologic_grade +
    pr_status,
  data = train_data
)[, -1] # Remove intercept

# Prepare Test Matrix (for later validation)
X_test <- model.matrix(
  ~ age_group + 
    chemotherapy_status + 
    radiotherapy_status + 
    tumor_stage +
    her2_status +
    histologic_grade +
    pr_status,
  data = test_data
)[, -1]

# Fit Multivariate Model
fg_multi <- crr(
  ftime = train_data$survival_days,
  fstatus = train_data$event_indicator2,
  cov1 = X_train,
  failcode = 1,
  cencode = 0
)

cat("\n========== MULTIVARIATE FINE-GRAY MODEL SUMMARY ==========\n")
summary(fg_multi)

# Extract Coefficients for reporting
coef_fg <- fg_multi$coef
hr_fg <- exp(coef_fg)
se_fg <- sqrt(diag(fg_multi$var))
ci_lower <- exp(coef_fg - 1.96 * se_fg)
ci_upper <- exp(coef_fg + 1.96 * se_fg)
p_values <- 2 * (1 - pnorm(abs(coef_fg / se_fg)))

results_fg <- data.frame(
  Variable = colnames(X_train),
  SHR = round(hr_fg, 3),
  CI_lower = round(ci_lower, 3),
  CI_upper = round(ci_upper, 3),
  p_value = round(p_values, 4)
)
print(results_fg)

# 6. COMPARISON: COX (CAUSE-SPECIFIC) VS FINE-GRAY -----------------------------

# Fit Cause-Specific Cox Model
cox_train <- coxph(
  Surv(survival_days, event_indicator) ~ 
    age_group + 
    chemotherapy_status + 
    radiotherapy_status + 
    tumor_stage +
    her2_status +
    histologic_grade +
    pr_status,
  data = train_data
)

# Extract Cox Results
cox_summary <- summary(cox_train)
cox_hr <- exp(coef(cox_train))
cox_p_values <- cox_summary$coefficients[, 5] 

# Create Comparison Table
fg_shr <- exp(fg_multi$coef)
fg_p_values <- summary(fg_multi)$coef[, 5]

comparison <- data.frame(
  Variable = names(coef(cox_train)),
  Cox_HR = round(cox_hr, 3),
  FG_SHR = round(fg_shr, 3),
  Cox_P = round(cox_p_values, 4),
  FG_P = round(fg_p_values, 4),
  Diff_Est = round(cox_hr - fg_shr, 3) # Difference in risk estimation
)

cat("\n========== COMPARISON: COX VS FINE-GRAY ==========\n")
print(comparison)

# 7. MODEL VALIDATION (C-INDEX & PREDICTION) -----------------------------------

# Helper function for Fine-Gray Predictions (CIF)
predict_crr <- function(model, newdata, times) {
  predictions <- list()
  for(t in times) {
    lp <- as.vector(newdata %*% model$coef)
    time_idx <- which(model$uftime <= t)
    
    if(length(time_idx) == 0) {
      cif <- rep(0, nrow(newdata))
    } else {
      baseline_cif <- model$bfitj[max(time_idx)]
      # Formula: CIF = 1 - (1 - baseline_CIF)^exp(lp)
      cif <- 1 - (1 - baseline_cif)^exp(lp)
    }
    predictions[[as.character(t)]] <- cif
  }
  return(predictions)
}

# Define Time Horizons (3, 5, 8 years in days)
times_pred <- c(3*365.25, 5*365.25, 8*365.25)

# Calculate C-Index for Test Set
pred_test <- predict_crr(fg_multi, X_test, times_pred)

cindex_results <- data.frame(
  Time_Horizon = c("3 Years", "5 Years", "8 Years"),
  C_Index = NA
)

for(i in 1:length(times_pred)) {
  # Create binary event for time t
  event_at_time <- ifelse(
    test_data$event_indicator2 == 1 & test_data$survival_days <= times_pred[i],
    1, 0
  )
  if(sum(event_at_time) > 0) {
    roc_obj <- pROC::roc(event_at_time, pred_test[[i]], quiet = TRUE)
    cindex_results$C_Index[i] <- round(as.numeric(roc_obj$auc), 3)
  }
}

cat("\n========== C-INDEX PERFORMANCE (TEST SET) ==========\n")
print(cindex_results)

# 8. CALIBRATION PLOTS ---------------------------------------------------------
# Note: pec/riskRegression require a Cause-Specific Cox object structure to 
# plot calibration for Competing Risks correctly.

formula_full <- Hist(survival_days, event_indicator2) ~ age_group + 
  chemotherapy_status + 
  radiotherapy_status + 
  tumor_stage +
  her2_status +
  histologic_grade +
  pr_status

# Fit CSC object for calibration function
fg_csc <- CSC(
  formula = formula_full,
  data = train_data,
  cause = 1 
)

# Plot Calibration for 3, 5, and 8 years (Training Data example)
# To plot for Test Data, change data = test_data

par(mfrow = c(1,3)) # 3 plots in one row

horizons <- c(3, 5, 8) * 365.25

for(h in horizons) {
  pec::calPlot(
    list("FineGray_CSC" = fg_csc),
    formula = formula_full,
    data = test_data,        # Validating on Test Data
    cause = 1,
    times = h,
    splitMethod = "none",
    main = paste("Calibration at", round(h/365.25), "Years")
  )
}
par(mfrow = c(1,1)) # Reset plot layout

# 9. NOMOGRAM CONSTRUCTION -----------------------------------------------------
# Based on Fine & Gray weights (Pseudo-observations method)

# A. Create Weighted Dataset
# The finegray() function creates a dataset where individuals with competing events
# remain in the risk set with decreasing weights.
fg_data <- finegray(Surv(survival_days, as.factor(event_indicator2)) ~ ., 
                    data = train_data, 
                    etype = 1) # 1 = Event of Interest

# B. Configure Datadist for the new weighted dataset
dd <- datadist(fg_data)
options(datadist = 'dd')

# C. Fit Weighted Cox Model
# We use the weights (fgwt) generated by the finegray function.
# This approximates the Subdistribution Hazard in a format rms understands.
fit_fg_nomo <- cph(Surv(fgstart, fgstop, fgstatus) ~ age_group + 
                     chemotherapy_status + 
                     radiotherapy_status + 
                     tumor_stage + 
                     her2_status + 
                     histologic_grade + 
                     pr_status,
                   data = fg_data, 
                   weight = fgwt,  # CRITICAL STEP: Apply weights
                   x = TRUE, y = TRUE, surv = TRUE)

# D. Define Survival/Risk Functions
surv_fg <- Survival(fit_fg_nomo)

# Convert Survival to Risk (Risk = 1 - Survival)
cif3_fg <- function(x) 1 - surv_fg(times_pred[1], x)
cif5_fg <- function(x) 1 - surv_fg(times_pred[2], x)
cif8_fg <- function(x) 1 - surv_fg(times_pred[3], x)

# E. Generate Nomogram Object
nom_final_fg <- nomogram(fit_fg_nomo, 
                         fun = list(cif3_fg, cif5_fg, cif8_fg),
                         funlabel = c("3-Year Risk", "5-Year Risk", "8-Year Risk"),
                         fun.at = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6), 
                         lp = FALSE,
                         vnames =