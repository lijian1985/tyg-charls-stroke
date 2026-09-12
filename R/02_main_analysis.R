# ============================================================
# 02_main_analysis.R
# Primary Cox models, tertile analysis and trend test
# ------------------------------------------------------------
# Input  : data/analysis_cohort.rds
# Output : results/table2_association.csv
#          results/table2_association.txt
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(readr)
})

cohort <- readRDS("data/analysis_cohort.rds")

# Complete-case analytic sample used for all adjusted models
cc <- cohort %>% filter(complete_case)

cat("Analytic sample N =", nrow(cc),
    "| events =", sum(cc$stroke_event), "\n\n")

# ------------------------------------------------------------
# Model definitions
# ------------------------------------------------------------
# Unadjusted : no covariates
# Model 1    : age, sex, area of residence
# Model 2    : Model 1 + alcohol, smoking, income category,
#              BMI group, CRP, platelet count, SBP, DBP

f_unadj <- Surv(followup_yrs, stroke_event) ~ TyG
f_m1    <- Surv(followup_yrs, stroke_event) ~ TyG + age + gender + rural
f_m2    <- Surv(followup_yrs, stroke_event) ~ TyG + age + gender + rural +
             drinking1_new + smoking1_new + income_total_cs + bmi_cs_1 +
             bl_crp + bl_plt + systo + diasto

models <- list(
  Unadjusted = f_unadj,
  `Model 1`  = f_m1,
  `Model 2`  = f_m2
)

# ------------------------------------------------------------
# 1. Continuous TyG: HR per 1-unit increase
# ------------------------------------------------------------
extract_hr <- function(fit, term) {
  s <- summary(fit)$coefficients
  ci <- summary(fit)$conf.int
  data.frame(
    term      = term,
    HR        = ci[term, "exp(coef)"],
    lower     = ci[term, "lower .95"],
    upper     = ci[term, "upper .95"],
    p         = s[term, "Pr(>|z|)"],
    row.names = NULL
  )
}

cont_results <- lapply(names(models), function(nm) {
  fit <- coxph(models[[nm]], data = cc)
  cbind(model = nm, extract_hr(fit, "TyG"))
}) %>% bind_rows()

cat("== HR per 1-unit increase in TyG ==\n")
for (i in seq_len(nrow(cont_results))) {
  with(cont_results[i, ], cat(sprintf(
    "  %-11s HR = %.2f (%.2f-%.2f), P = %.4f\n",
    model, HR, lower, upper, p)))
}

# ------------------------------------------------------------
# 2. TyG tertiles (low tertile as reference)
# ------------------------------------------------------------
cc <- cc %>% mutate(tyg_tertile = factor(tyg_tertile, levels = c(1, 2, 3)))

f_unadj_t <- Surv(followup_yrs, stroke_event) ~ tyg_tertile
f_m1_t    <- Surv(followup_yrs, stroke_event) ~ tyg_tertile + age + gender + rural
f_m2_t    <- Surv(followup_yrs, stroke_event) ~ tyg_tertile + age + gender + rural +
               drinking1_new + smoking1_new + income_total_cs + bmi_cs_1 +
               bl_crp + bl_plt + systo + diasto

tert_models <- list(
  Unadjusted = f_unadj_t,
  `Model 1`  = f_m1_t,
  `Model 2`  = f_m2_t
)

tert_results <- lapply(names(tert_models), function(nm) {
  fit <- coxph(tert_models[[nm]], data = cc)
  r <- extract_hr(fit, "tyg_tertile2")
  r <- rbind(r, extract_hr(fit, "tyg_tertile3"))
  r$level <- c("Middle vs Low", "High vs Low")
  cbind(model = nm, r)
}) %>% bind_rows()

cat("\n== TyG tertiles ==\n")
for (i in seq_len(nrow(tert_results))) {
  with(tert_results[i, ], cat(sprintf(
    "  %-11s %-14s HR = %.2f (%.2f-%.2f), P = %.4f\n",
    model, level, HR, lower, upper, p)))
}

# ------------------------------------------------------------
# 3. P for trend (tertile as a continuous score, per tertile step)
# ------------------------------------------------------------
trend_results <- lapply(names(tert_models), function(nm) {
  f <- update(tert_models[[nm]], . ~ . - tyg_tertile + as.numeric(tyg_tertile))
  fit <- coxph(f, data = cc)
  cbind(model = nm,
        extract_hr(fit, "as.numeric(tyg_tertile)"),
        level = "Trend per tertile")
}) %>% bind_rows()

cat("\n== P for trend ==\n")
for (i in seq_len(nrow(trend_results))) {
  with(trend_results[i, ], cat(sprintf(
    "  %-11s HR = %.2f (%.2f-%.2f), P = %.4f\n",
    model, HR, lower, upper, p)))
}

# ------------------------------------------------------------
# 4. Incidence rates by tertile
# ------------------------------------------------------------
rates <- cc %>%
  group_by(tyg_tertile) %>%
  summarise(
    n        = n(),
    events   = sum(stroke_event),
    person_y = sum(followup_yrs),
    rate_1000py = 1000 * sum(stroke_event) / sum(followup_yrs),
    .groups = "drop"
  )

cat("\n== Incidence rate per 1,000 person-years ==\n")
print(as.data.frame(rates))

# ------------------------------------------------------------
# 5. Save
# ------------------------------------------------------------
all_results <- bind_rows(
  cbind(cont_results, level = "Per 1-unit increase") %>% select(model, level, term, HR, lower, upper, p),
  tert_results %>% select(model, level, term, HR, lower, upper, p),
  trend_results %>% select(model, level, term, HR, lower, upper, p)
)

dir.create("results", showWarnings = FALSE)
write_csv(all_results, "results/table2_association.csv")
write_csv(rates, "results/incidence_rates_by_tertile.csv")

cat("\nSaved: results/table2_association.csv\n")
cat("Saved: results/incidence_rates_by_tertile.csv\n")
