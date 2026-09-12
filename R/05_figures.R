# ============================================================
# 05_figures.R
# Reproduces Figures 2-4 and the supplementary TyG distribution
# ------------------------------------------------------------
# Input  : data/analysis_cohort.rds
#          results/rcs_curve.csv
# Output : figures/*.pdf  and  figures/*.png
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(readr)
})

cohort <- readRDS("data/analysis_cohort.rds")
cc     <- cohort %>% filter(complete_case)
dir.create("figures", showWarnings = FALSE)

# ------------------------------------------------------------
# Figure 2. Kaplan-Meier curves by TyG tertile
# ------------------------------------------------------------
cc <- cc %>%
  mutate(tyg_tertile = factor(tyg_tertile, levels = c(1, 2, 3),
                              labels = c("Low", "Middle", "High")))

km_fit <- survfit(Surv(followup_yrs, stroke_event) ~ tyg_tertile, data = cc)

p_km <- ggsurvplot(
  km_fit,
  data        = cc,
  risk.table  = TRUE,
  pval        = TRUE,
  conf.int    = FALSE,
  xlab        = "Follow-up time (years)",
  ylab        = "Stroke-free survival probability",
  legend.title = "TyG tertile",
  legend.labs  = c("Low", "Middle", "High"),
  palette      = c("#4C72B0", "#55A868", "#C44E52"),
  break.time.by = 2,
  ggtheme     = theme_classic(base_size = 12),
  tables.height = 0.22
)

ggsave("figures/Figure2_KM_TyG_tertiles.pdf", plot = p_km$plot, width = 7, height = 5)
ggsave("figures/Figure2_KM_TyG_tertiles.png", plot = p_km$plot, width = 7, height = 5, dpi = 600)

# ------------------------------------------------------------
# Figure 3. Restricted cubic spline
# ------------------------------------------------------------
rcs <- read_csv("results/rcs_curve.csv", show_col_types = FALSE)
rcs <- rcs %>% filter(!is.na(yhat))

p_rcs <- ggplot(rcs, aes(x = TyG, y = yhat)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.18, fill = "#4C72B0") +
  geom_line(linewidth = 0.9, colour = "#1F4E79") +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey45") +
  geom_vline(xintercept = 8.3, linetype = 3, colour = "#C44E52") +
  scale_y_log10() +
  labs(x = "TyG index", y = "Hazard ratio (95% CI) for incident stroke") +
  theme_classic(base_size = 12)

ggsave("figures/Figure3_TyG_RCS.pdf", plot = p_rcs, width = 7, height = 5)
ggsave("figures/Figure3_TyG_RCS.png", plot = p_rcs, width = 7, height = 5, dpi = 600)

# ------------------------------------------------------------
# Figure 4. Subgroup forest plot
# ------------------------------------------------------------
sg <- read_csv("results/subgroup_analysis.csv", show_col_types = FALSE)

sg <- sg %>%
  mutate(label = paste0(subgroup, ": ", level)) %>%
  arrange(HR)

p_forest <- ggplot(sg, aes(x = HR, y = reorder(label, HR))) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey45") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.22, colour = "#1F4E79") +
  geom_point(size = 2.1, colour = "#1F4E79") +
  scale_x_log10() +
  labs(x = "Hazard ratio per 1-unit increase in TyG (95% CI)", y = NULL) +
  theme_classic(base_size = 11)

ggsave("figures/Figure4_Subgroup_Forest.pdf", plot = p_forest, width = 8, height = 6)
ggsave("figures/Figure4_Subgroup_Forest.png", plot = p_forest, width = 8, height = 6, dpi = 600)

# ------------------------------------------------------------
# Supplementary figure. TyG distribution with events per bin
# ------------------------------------------------------------
bins <- cc %>%
  mutate(bin = cut(TyG, breaks = seq(7, 13, by = 0.5), right = FALSE)) %>%
  group_by(bin) %>%
  summarise(n = n(), events = sum(stroke_event), .groups = "drop") %>%
  mutate(mid = as.numeric(sub("\\[(.*),.*", "\\1", bin)) + 0.25)

p_dist <- ggplot(bins, aes(x = mid)) +
  geom_col(aes(y = n), fill = "grey80", width = 0.45) +
  geom_point(aes(y = events * 20), colour = "#C44E52", size = 2) +
  scale_y_continuous(
    name = "Participants per bin",
    sec.axis = sec_axis(~ . / 20, name = "Stroke events per bin")
  ) +
  labs(x = "TyG index") +
  theme_classic(base_size = 12)

ggsave("figures/Supplementary_Figure_TyG_Distribution.pdf", plot = p_dist, width = 7, height = 5)
ggsave("figures/Supplementary_Figure_TyG_Distribution.png", plot = p_dist, width = 7, height = 5, dpi = 600)

cat("Figures written to figures/\n")
