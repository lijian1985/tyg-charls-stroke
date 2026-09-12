# ============================================================
# 03_threshold_analysis.R
# Restricted cubic splines and two-piecewise Cox regression
# ------------------------------------------------------------
# Input  : data/analysis_cohort.rds
# Output : results/table3_threshold.csv
#          results/rcs_curve.csv
# ============================================================
#
# NOTE ON INTERPRETATION
# ----------------------
# Two separate tests are reported and must not be conflated:
#
#   (a) RCS log-likelihood ratio test for non-linearity:
#       compares a 5-knot restricted cubic spline against the
#       linear model. This tests whether the dose-response
#       curve deviates from linearity overall.
#
#   (b) Two-piecewise model log-likelihood ratio test:
#       compares a broken-line (2-segment) model against the
#       linear model. This tests whether allowing a change in
#       slope at the inflection point improves fit.
#
# They are reported separately in the manuscript because they
# answer different questions and can disagree.

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(rms)
  library(readr)
})

options(datadist = NULL)

cohort <- readRDS("data/analysis_cohort.rds")
cc     <- cohort %>% filter(complete_case)

# The two-piecewise search is performed over this grid.
# The grid is restricted to the interior of the distribution so
# that neither segment is estimated from a handful of
# observations alone.
knot_grid <- seq(7.5, 10.0, by = 0.1)

# ------------------------------------------------------------
# 1. Restricted cubic spline (5 knots) vs linear model
# ------------------------------------------------------------
dd <- datadist(cc)
options(datadist = "dd")

f_lin <- cph(Surv(followup_yrs, stroke_event) ~ TyG + age + gender + rural +
               drinking1_new + smoking1_new + income_total_cs + bmi_cs_1 +
               bl_crp + bl_plt + systo + diasto,
             data = cc, x = TRUE, y = TRUE, surv = TRUE)

f_rcs <- update(f_lin, . ~ . - TyG + rcs(TyG, 5))

ll_lin <- f_lin$loglik[2]
ll_rcs <- f_rcs$loglik[2]

# Non-linearity is tested with the extra degrees of freedom used
# by the spline terms.
df_diff <- length(coef(f_rcs)) - length(coef(f_lin))
p_nonlin <- pchisq(2 * (ll_rcs - ll_lin), df = df_diff, lower.tail = FALSE)

cat("== Restricted cubic spline ==\n")
cat(sprintf("  Linear model log-likelihood : %.4f\n", ll_lin))
cat(sprintf("  RCS (5 knots) log-likelihood: %.4f\n", ll_rcs))
cat(sprintf("  Extra df                    : %d\n", df_diff))
cat(sprintf("  P for non-linearity         : %.4f\n\n", p_nonlin))

# Predicted HR curve relative to a reference value.
ref_tyg <- 8.3
TyG_seq <- seq(quantile(cc$TyG, 0.01), quantile(cc$TyG, 0.99), length.out = 200)

pred <- Predict(f_rcs, TyG = TyG_seq, ref.zero = TRUE, fun = exp)
pred <- as.data.frame(pred)
pred <- pred %>% mutate(TyG = TyG + ref_tyg)

dir.create("results", showWarnings = FALSE)
write_csv(pred, "results/rcs_curve.csv")
cat("Saved: results/rcs_curve.csv\n\n")

# ------------------------------------------------------------
# 2. Two-piecewise Cox regression
# ------------------------------------------------------------
# Segment 1: (TyG - k) for TyG <  k, 0 otherwise
# Segment 2: (TyG - k) for TyG >= k, 0 otherwise
#
# The slopes of the two segments are estimated simultaneously in a
# single Cox model, so both segments are adjusted for the same
# covariates and the log-likelihoods are directly comparable.

fit_piecewise <- function(k) {
  d <- cc
  d$seg1 <- ifelse(d$TyG <  k, d$TyG - k, 0)
  d$seg2 <- ifelse(d$TyG >= k, d$TyG - k, 0)

  fit <- coxph(
    Surv(followup_yrs, stroke_event) ~ seg1 + seg2 +
      age + gender + rural + drinking1_new + smoking1_new +
      income_total_cs + bmi_cs_1 + bl_crp + bl_plt + systo + diasto,
    data = d
  )
  list(k = k, fit = fit, loglik = fit$loglik[2])
}

grid_fits <- lapply(knot_grid, fit_piecewise)
logliks   <- sapply(grid_fits, function(x) x$loglik)
best      <- grid_fits[[which.max(logliks)]]
k_opt     <- best$k

cat("== Two-piecewise Cox regression ==\n")
cat(sprintf("  Optimal inflection point: TyG = %.1f\n", k_opt))

# Compare the two-piecewise model against the linear model on the
# same sample. The linear model is refitted here without the rms
# wrapper so the log-likelihoods come from identical fitting code.
f_lin_plain <- coxph(
  Surv(followup_yrs, stroke_event) ~ TyG + age + gender + rural +
    drinking1_new + smoking1_new + income_total_cs + bmi_cs_1 +
    bl_crp + bl_plt + systo + diasto,
  data = cc
)

ll_pw   <- best$loglik
ll_flat <- f_lin_plain$loglik[2]
lr_stat <- 2 * (ll_pw - ll_flat)
p_lrt   <- pchisq(lr_stat, df = 1, lower.tail = FALSE)

cat(sprintf("  Two-piecewise log-likelihood: %.4f\n", ll_pw))
cat(sprintf("  Linear log-likelihood       : %.4f\n", ll_flat))
cat(sprintf("  LRT chi-square (df = 1)     : %.3f, P = %.4f\n\n",
            lr_stat, p_lrt))

s  <- summary(best$fit)
ci <- summary(best$fit)$conf.int

seg_rows <- data.frame(
  segment = c("Below threshold", "Above threshold"),
  HR      = c(ci["seg1", "exp(coef)"], ci["seg2", "exp(coef)"]),
  lower   = c(ci["seg1", "lower .95"],  ci["seg2", "lower .95"]),
  upper   = c(ci["seg1", "upper .95"],  ci["seg2", "upper .95"]),
  p       = c(s$coefficients["seg1", "Pr(>|z|)"],
              s$coefficients["seg2", "Pr(>|z|)"])
)

print(seg_rows)

# ------------------------------------------------------------
# 3. Test for a difference in slopes between segments
# ------------------------------------------------------------
# This is a separate Wald test and its P value is NOT the same as
# the log-likelihood ratio test in section 2. Both are reported in
# the manuscript.
V <- vcov(best$fit)
diff_beta <- coef(best$fit)["seg2"] - coef(best$fit)["seg1"]
se_diff   <- sqrt(V["seg1", "seg1"] + V["seg2", "seg2"] - 2 * V["seg1", "seg2"])
z_diff    <- diff_beta / se_diff
p_diff    <- 2 * pnorm(-abs(z_diff))

cat(sprintf("\nDifference in slopes: z = %.3f, P = %.4f\n", z_diff, p_diff))

# ------------------------------------------------------------
# 4. Distribution of the low-TyG tail
# ------------------------------------------------------------
cat("\n== Low-TyG tail (complete-case sample) ==\n")
for (cut in c(7.0, 7.5, 8.0, 8.3)) {
  sub <- cc %>% filter(TyG < cut)
  cat(sprintf("  TyG < %.1f : n = %5d, events = %4d\n",
              cut, nrow(sub), sum(sub$stroke_event)))
}

# ------------------------------------------------------------
# 5. Save
# ------------------------------------------------------------
out <- rbind(
  data.frame(test = "RCS non-linearity (5 knots)",
             statistic = 2 * (ll_rcs - ll_lin), df = df_diff,
             p = p_nonlin, note = "vs linear model"),
  data.frame(test = "Two-piecewise LRT",
             statistic = lr_stat, df = 1,
             p = p_lrt, note = sprintf("inflection at TyG = %.1f", k_opt)),
  data.frame(test = "Difference in segment slopes",
             statistic = z_diff, df = NA,
             p = p_diff, note = "Wald z-test")
)

write_csv(out, "results/table3_threshold.csv")
write_csv(seg_rows, "results/two_piecewise_segments.csv")

cat("\nSaved: results/table3_threshold.csv\n")
cat("Saved: results/two_piecewise_segments.csv\n")
