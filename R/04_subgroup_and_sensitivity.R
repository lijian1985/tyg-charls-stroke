# ============================================================
# 04_subgroup_and_sensitivity.R
# Subgroup analyses, interaction tests and sensitivity analyses
# ------------------------------------------------------------
# Input  : data/analysis_cohort.rds
# Output : results/subgroup_analysis.csv
#          results/sensitivity_analysis.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(readr)
})

cohort <- readRDS("data/analysis_cohort.rds")
cc     <- cohort %>% filter(complete_case)

covars <- "age + gender + rural + drinking1_new + smoking1_new +
           income_total_cs + bmi_cs_1 + bl_crp + bl_plt + systo + diasto"

# ------------------------------------------------------------
# 1. Subgroup analyses
# ------------------------------------------------------------
# Prespecified subgroups. Age is dichotomised at 60 years as in
# the manuscript.
cc <- cc %>%
  mutate(
    age_group = factor(if_else(age < 60, "<60", ">=60"), levels = c("<60", ">=60")),
    hypertension = factor(hypertension_lbl),
    diabetes     = factor(diabetes_lbl),
    dyslipidemia = factor(dyslipidemia_lbl)
  )

subgroups <- list(
  "Age"             = "age_group",
  "Sex"             = "gender",
  "Area of residence" = "rural",
  "Smoking"         = "smoking1_new",
  "Alcohol"         = "drinking1_new",
  "Hypertension"    = "hypertension",
  "Diabetes"        = "diabetes",
  "Dyslipidaemia"   = "dyslipidemia"
)

subgroup_rows <- list()

for (label in names(subgroups)) {
  v <- subgroups[[label]]
  levs <- levels(factor(cc[[v]]))

  for (lv in levs) {
    sub <- cc %>% filter(as.character(.data[[v]]) == lv)
    if (nrow(sub) < 50 || sum(sub$stroke_event) < 10) next

    fit <- coxph(
      as.formula(paste("Surv(followup_yrs, stroke_event) ~ TyG +", covars)),
      data = sub
    )
    ci <- summary(fit)$conf.int
    co <- summary(fit)$coefficients

    subgroup_rows[[length(subgroup_rows) + 1]] <- data.frame(
      subgroup = label, level = lv,
      n = nrow(sub), events = sum(sub$stroke_event),
      HR = ci["TyG", "exp(coef)"],
      lower = ci["TyG", "lower .95"],
      upper = ci["TyG", "upper .95"],
      p = co["TyG", "Pr(>|z|)"]
    )
  }
}

subgroup_df <- bind_rows(subgroup_rows)

# ------------------------------------------------------------
# 2. Interaction tests
# ------------------------------------------------------------
# Likelihood ratio test comparing a model with the exposure-by-
# subgroup interaction term against the model without it.
interaction_rows <- list()

for (label in names(subgroups)) {
  v <- subgroups[[label]]

  f0 <- as.formula(paste("Surv(followup_yrs, stroke_event) ~ TyG +", v, "+", covars))
  f1 <- as.formula(paste("Surv(followup_yrs, stroke_event) ~ TyG *", v, "+", covars))

  m0 <- coxph(f0, data = cc)
  m1 <- coxph(f1, data = cc)

  df_diff <- length(coef(m1)) - length(coef(m0))
  lr <- 2 * (m1$loglik[2] - m0$loglik[2])
  p  <- pchisq(lr, df = df_diff, lower.tail = FALSE)

  interaction_rows[[length(interaction_rows) + 1]] <- data.frame(
    subgroup = label, chi_sq = lr, df = df_diff, p_interaction = p
  )
}

interaction_df <- bind_rows(interaction_rows)

cat("== Subgroup HRs (per 1-unit TyG, fully adjusted) ==\n")
for (i in seq_len(nrow(subgroup_df))) {
  with(subgroup_df[i, ], cat(sprintf(
    "  %-18s %-6s n = %4d, ev = %3d, HR = %.2f (%.2f-%.2f)\n",
    subgroup, level, n, events, HR, lower, upper)))
}

cat("\n== Interaction tests ==\n")
for (i in seq_len(nrow(interaction_df))) {
  with(interaction_df[i, ], cat(sprintf(
    "  %-18s LRT chi2 = %6.3f, df = %d, P = %.4f\n",
    subgroup, chi_sq, df, p_interaction)))
}

# ------------------------------------------------------------
# 3. Sensitivity analyses
# ------------------------------------------------------------
sens_rows <- list()

add_sens <- function(label, data, note) {
  fit <- coxph(
    as.formula(paste("Surv(followup_yrs, stroke_event) ~ TyG +", covars)),
    data = data
  )
  ci <- summary(fit)$conf.int
  co <- summary(fit)$coefficients
  sens_rows[[length(sens_rows) + 1]] <<- data.frame(
    analysis = label,
    n = nrow(data), events = sum(data$stroke_event),
    HR = ci["TyG", "exp(coef)"],
    lower = ci["TyG", "lower .95"],
    upper = ci["TyG", "upper .95"],
    p = co["TyG", "Pr(>|z|)"],
    note = note
  )
}

# (a) primary
add_sens("Primary (complete case)", cc, "reference")

# (b) excluding follow-up < 2 years
add_sens("Follow-up >= 2 years",
         cc %>% filter(followup_yrs >= 2),
         "excludes early events / short follow-up")

# (c) excluding events within the first 2 years
add_sens("Excluding events in first 2 years",
         cc %>% filter(!(stroke_event == 1 & followup_yrs < 2)),
         "addresses reverse causation")

# (d) excluding sparse low-TyG tail
add_sens("Excluding TyG < 7.0",
         cc %>% filter(TyG >= 7.0),
         "sparse tail removed")
add_sens("Excluding TyG < 7.5",
         cc %>% filter(TyG >= 7.5),
         "sparse tail removed")

sens_df <- bind_rows(sens_rows)

cat("\n== Sensitivity analyses ==\n")
for (i in seq_len(nrow(sens_df))) {
  with(sens_df[i, ], cat(sprintf(
    "  %-34s n = %4d, ev = %3d, HR = %.2f (%.2f-%.2f), P = %.4f\n",
    analysis, n, events, HR, lower, upper, p)))
}

# ------------------------------------------------------------
# 4. Save
# ------------------------------------------------------------
write_csv(subgroup_df,   "results/subgroup_analysis.csv")
write_csv(interaction_df, "results/interaction_tests.csv")
write_csv(sens_df,       "results/sensitivity_analysis.csv")

cat("\nSaved: results/subgroup_analysis.csv\n")
cat("Saved: results/interaction_tests.csv\n")
cat("Saved: results/sensitivity_analysis.csv\n")
