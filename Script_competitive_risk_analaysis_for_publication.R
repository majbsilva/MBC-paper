# ==============================================================================
# Competing Risk Analysis Using Fine & Gray Models
# ==============================================================================
#
# Description:
# This script implements a full competing risk analysis pipeline using
# Fine & Gray subdistribution hazard models, including:
#  - Data preprocessing
#  - Univariate and multivariable modeling
#  - Comparison with cause-specific Cox models
#  - Model discrimination using time-dependent C-index
#  - Calibration under a competing risk framework
#  - Nomogram construction 
#
# Data availability:
# The dataset used in this analysis originates from third-party sources
# (SEER Program and FOSP Hospital-Based Cancer Registry) and cannot be
# redistributed. Instructions for data access are provided in the README file.
#
# Software:
# R >= 4.3.0
#
# Author:
# Marcelo José Barbosa Silva
#
# License:
# MIT License
# ==============================================================================


## ============================================================================
## 0. PACKAGES
## ============================================================================

library(cmprsk)
library(riskRegression)
library(pec)
library(rms)
library(survival)
library(dplyr)
library(readr)


## ============================================================================
## 1. DATA INPUT AND PREPARATION
## ============================================================================

# The dataset must be placed in the 'data/' directory after being obtained
# directly from the original registries and harmonized as described
# in the manuscript.

data_path <- "data/mbc.csv"

if (!file.exists(data_path)) {
  stop("Dataset not found. Place 'mbc.csv' in the 'data/' directory.")
}

mbc_raw <- read_csv(data_path)

vars_required <- c(
  "survival_days",
  "event_indicator",    # cause-specific event
  "event_indicator2",   # competing risks indicator
  "age",
  "age_group",
  "residence_region",
  "tumor_stage",
  "chemotherapy_status",
  "radiotherapy_status",
  "surgery_status",
  "TTS_CAT60"
)

mbc <- mbc_raw[, vars_required] |> na.omit()

cat("Original observations:", nrow(mbc_raw), "\n")
cat("Complete cases:", nrow(mbc), "\n")
cat("Excluded due to missingness:", nrow(mbc_raw) - nrow(mbc), "\n")


## ============================================================================
## 2. TRAIN / TEST SPLIT
## ============================================================================

# The split allows internal validation of model performance.

set.seed(123)

train_idx <- sample(seq_len(nrow(mbc)), size = 0.7 * nrow(mbc))
train_data <- mbc[train_idx, ]
test_data  <- mbc[-train_idx, ]

cat("\nTraining set:", nrow(train_data),
    "- Events:", sum(train_data$event_indicator2 == 1), "\n")
cat("Testing set:", nrow(test_data),
    "- Events:", sum(test_data$event_indicator2 == 1), "\n")


## ============================================================================
## 3. UNIVARIATE FINE & GRAY MODELS
## ============================================================================

# Each covariate is evaluated independently to explore its association
# with the cumulative incidence of the event of interest.

sink("output/fine_gray_univariate_results.txt")

run_fg_uni <- function(var, data) {
  X <- model.matrix(as.formula(paste("~", var)), data)[, -1, drop = FALSE]
  crr(
    ftime    = data$survival_days,
    fstatus  = data$event_indicator2,
    cov1     = X,
    failcode = 1,
    cencode  = 0
  )
}

covariates_uni <- c(
  "age_group",
  "chemotherapy_status",
  "surgery_status",
  "radiotherapy_status",
  "residence_region",
  "tumor_stage",
  "TTS_CAT60"
)

for (v in covariates_uni) {
  cat("\n--- Univariate Fine & Gray:", v, "---\n")
  print(summary(run_fg_uni(v, train_data)))
}

sink()


## ============================================================================
## 4. MULTIVARIABLE FINE & GRAY MODEL
## ============================================================================

# Covariates are selected a priori based on clinical relevance.

X_train <- model.matrix(
  ~ age_group + chemotherapy_status + surgery_status +
    residence_region + tumor_stage,
  data = train_data
)[, -1]

X_test <- model.matrix(
  ~ age_group + chemotherapy_status + surgery_status +
    residence_region + tumor_stage,
  data = test_data
)[, -1]

fg_multi <- crr(
  ftime    = train_data$survival_days,
  fstatus  = train_data$event_indicator2,
  cov1     = X_train,
  failcode = 1,
  cencode  = 0
)

cat("\n========== MULTIVARIABLE FINE & GRAY MODEL ==========\n")
print(summary(fg_multi))


## ============================================================================
## 5. COMPARISON WITH CAUSE-SPECIFIC COX MODEL
## ============================================================================

# The Cox model is fitted for illustrative comparison of risk interpretation.

cox_model <- coxph(
  Surv(survival_days, event_indicator) ~
    age_group + chemotherapy_status + surgery_status +
    residence_region + tumor_stage,
  data = train_data
)

comparison <- data.frame(
  Variable = names(coef(cox_model)),
  Cox_HR   = round(exp(coef(cox_model)), 3),
  FG_SHR   = round(exp(fg_multi$coef), 3),
  Cox_P   = round(summary(cox_model)$coefficients[, 5], 4),
  FG_P    = round(summary(fg_multi)$coef[, 5], 4)
)

cat("\n========== COX VS FINE & GRAY ==========\n")
print(comparison)


## ============================================================================
## 6. MODEL DISCRIMINATION AND CALIBRATION FRAMEWORK
## ============================================================================

# The Fine & Gray model from cmprsk does not expose standardized prediction
# methods compatible with survival-based performance metrics. Therefore,
# a cause-specific Cox (CSC) model is fitted using the riskRegression framework.
#
# The Hist() interface explicitly accommodates multiple competing events and
# enables estimation of cumulative incidence functions required for
# discrimination and calibration analyses.

library(riskRegression)

formula_hist <- Hist(survival_days, event_indicator2) ~
  age_group + chemotherapy_status + surgery_status +
  residence_region + tumor_stage

fg_csc <- CSC(
  formula = formula_hist,
  data    = train_data,
  cause   = 1
)


## ============================================================================
## 7. TIME-DEPENDENT C-INDEX (TEST SET)
## ============================================================================

# Discrimination is assessed using the time-dependent concordance index (C-index),
# accounting for censoring and competing risks.

library(pec)

times_cindex <- c(3, 5, 8) * 365.25

cindex_obj <- pec::cindex(
  object = fg_csc,
  formula = formula_hist,
  data = test_data,
  eval.times = times_cindex,
  cause = 1,
  splitMethod = "none"
)

cindex_results <- data.frame(
  Time = c("3 years", "5 years", "8 years"),
  C_index = round(cindex_obj$AppCindex, 3)
)

cat("\n========== TIME-DEPENDENT C-INDEX (TEST SET) ==========\n")
print(cindex_results)


## ============================================================================
## 8. CALIBRATION
## ============================================================================

# Calibration curves compare predicted and observed cumulative incidence
# probabilities under a competing risk framework.

calPlot(
  list("Training set" = fg_csc),
  formula = formula_hist,
  data = train_data,
  cause = 1,
  times = 3 * 365.25,
  splitMethod = "none"
)

calPlot(
  list("Testing set" = fg_csc),
  formula = formula_hist,
  data = test_data,
  cause = 1,
  times = 3 * 365.25,
  splitMethod = "none"
)


## ============================================================================
## 9. NOMOGRAM (CAUSE-SPECIFIC MODEL)
## ============================================================================

# The nomogram provides individualized cause-specific risk estimates and is
# intended for clinical illustration purposes.

label(train_data$age_group)           <- "Age group"
label(train_data$chemotherapy_status) <- "Chemotherapy"
label(train_data$surgery_status)      <- "Surgery"
label(train_data$residence_region)    <- "Population size"
label(train_data$tumor_stage)         <- "Tumor stage"

dd <- datadist(train_data)
options(datadist = "dd")

fit_nomo <- cph(
  Surv(survival_days, event_indicator) ~
    age_group + chemotherapy_status + surgery_status +
    residence_region + tumor_stage,
  data = train_data,
  x = TRUE, y = TRUE, surv = TRUE
)

surv_fn <- Survival(fit_nomo)

nom <- nomogram(
  fit_nomo,
  fun = list(
    function(x) 1 - surv_fn(3 * 365.25, x),
    function(x) 1 - surv_fn(5 * 365.25, x),
    function(x) 1 - surv_fn(8 * 365.25, x)
  ),
  funlabel = c("3-year risk", "5-year risk", "8-year risk")
)

png("output/nomogram_fosp.png", width = 3500, height = 2500, res = 300)
plot(nom)
dev.off()
