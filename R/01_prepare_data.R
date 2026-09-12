# ============================================================
# 01_prepare_data.R
# Triglyceride-glucose (TyG) index and incident stroke in
# middle-aged and older Chinese adults with chronic diseases
# ------------------------------------------------------------
# Purpose : Build the analysis dataset from CHARLS Wave 1 (2011)
#           baseline data, derive the TyG index, and apply
#           the inclusion / exclusion criteria described in
#           the manuscript.
# Input   : CHARLS_2011.csv  (not redistributed; see README)
# Output  : analysis_cohort.rds
# ============================================================

# ------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

set.seed(20260912)

# Path to the CHARLS baseline extract. Obtained by the user from
# http://charls.pku.edu.cn/en/ after registration. The exact
# variable names below follow the CHARLS harmonised release.
input_file  <- "data/CHARLS_2011.csv"
output_file <- "data/analysis_cohort.rds"

# ------------------------------------------------------------
# 1. Read raw data
# ------------------------------------------------------------
raw <- read_csv(input_file, show_col_types = FALSE)

# ------------------------------------------------------------
# 2. Derive the TyG index
# ------------------------------------------------------------
# TyG = ln[ fasting triglycerides (mg/dL) x fasting glucose (mg/dL) / 2 ]
#
# CHARLS reports triglycerides and glucose in mg/dL in the
# harmonised blood-biomarker file (bl_tg, bl_glu). If your extract
# stores them in mmol/L, convert first:
#   TG  (mg/dL) = TG  (mmol/L) x 88.57
#   GLU (mg/dL) = GLU (mmol/L) x 18.02
#
# A guard is included so the script fails loudly rather than
# silently producing implausible values.

dat <- raw %>%
  mutate(
    tg_mgdl  = bl_tg,
    glu_mgdl = bl_glu,
    TyG      = log(tg_mgdl * glu_mgdl / 2)
  )

stopifnot(
  "TyG values outside the plausible 4-15 range - check units" =
    all(dat$TyG > 4 & dat$TyG < 15, na.rm = TRUE)
)

# ------------------------------------------------------------
# 3. Chronic disease indicator
# ------------------------------------------------------------
# At least one self-reported physician-diagnosed chronic condition
# among the 13 conditions listed in the manuscript.
chronic_vars <- c(
  "hibpe",   # hypertension
  "diabe",   # diabetes
  "dyslipe", # dyslipidaemia
  "cancre",  # cancer (excluding minor skin cancers)
  "lunge",   # chronic lung disease
  "livere",  # liver disease
  "hearte",  # heart disease
  "kidneye", # kidney disease
  "digeste", # digestive disease
  "psyche",  # emotional / psychiatric disorders
  "memrye",  # memory-related disorders
  "arthre",  # arthritis or rheumatism
  "asthmae"  # asthma
)

missing_all_chronic <- dat %>%
  select(all_of(chronic_vars)) %>%
  is.na() %>%
  apply(1, all)

dat <- dat %>%
  mutate(
    n_chronic = rowSums(select(., all_of(chronic_vars)) == 1, na.rm = TRUE),
    any_chronic = if_else(missing_all_chronic, NA_real_, as.numeric(n_chronic >= 1))
  )

# ------------------------------------------------------------
# 4. Follow-up time and event indicator
# ------------------------------------------------------------
# Incident stroke is defined as a first self-reported
# physician-diagnosed stroke during Waves 2-5 (2013-2020).
# Follow-up time is expressed in years.

dat <- dat %>%
  mutate(
    stroke_event = stroke_y,
    followup_yrs = stroke_time / 365.25
  )

# ------------------------------------------------------------
# 5. Apply inclusion / exclusion criteria
# ------------------------------------------------------------
# Order of exclusions as reported in Figure 1 of the manuscript:
#   17,708 assessed
#     -648  age <45 or missing age
#     -5,500 missing chronic disease history
#     -167  prevalent stroke at baseline
#     -158  history of cancer
#     -3,567 unavailable TyG data
#     -452  prior stroke, missing data, or loss to follow-up
#   = 7,216 analysed

n_assessed <- nrow(dat)

step1 <- dat %>% filter(!is.na(age), age >= 45)
step2 <- step1 %>% filter(!is.na(any_chronic))
step3 <- step2 %>% filter(is.na(stroke) | stroke != 1)   # no prevalent stroke
step4 <- step3 %>% filter(is.na(cancre) | cancre != 1)   # no cancer history
step5 <- step4 %>% filter(!is.na(TyG))                   # TyG available
step6 <- step5 %>% filter(!is.na(stroke_event),
                          !is.na(followup_yrs),
                          followup_yrs > 0)

cohort <- step6

cat("Participants assessed :", n_assessed, "\n")
cat("After age exclusion   :", nrow(step1), "\n")
cat("After chronic disease :", nrow(step2), "\n")
cat("After prevalent stroke:", nrow(step3), "\n")
cat("After cancer exclusion:", nrow(step4), "\n")
cat("After TyG available   :", nrow(step5), "\n")
cat("Final eligible cohort :", nrow(cohort), "\n")
cat("Incident strokes      :", sum(cohort$stroke_event), "\n")

# ------------------------------------------------------------
# 6. Covariates used in the fully adjusted model
# ------------------------------------------------------------
# Model 2 covariates, matching the manuscript:
#   age, sex, area of residence, alcohol consumption, smoking
#   status, annual income category, BMI group, C-reactive
#   protein, platelet count, systolic and diastolic blood pressure

model2_covariates <- c(
  "age", "gender", "rural",
  "drinking1_new", "smoking1_new",
  "income_total_cs", "bmi_cs_1",
  "bl_crp", "bl_plt", "systo", "diasto"
)

cohort <- cohort %>%
  mutate(
    complete_case = if_else(
      rowSums(is.na(select(., all_of(model2_covariates)))) == 0,
      TRUE, FALSE
    ),
    # TyG tertiles are defined within the complete-case sample so
    # that Table 2 and Figure 2 describe the same participants as
    # the adjusted models.
    tyg_tertile = if_else(
      complete_case,
      ntile(TyG[complete_case], 3),
      NA_integer_
    )
  )

cat("\nComplete-case sample  :", sum(cohort$complete_case), "\n")
cat("Events in complete case:", sum(cohort$stroke_event[cohort$complete_case]), "\n")

# ------------------------------------------------------------
# 7. Save
# ------------------------------------------------------------
saveRDS(cohort, output_file)
cat("\nSaved:", output_file, "\n")
