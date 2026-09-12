# TyG index and incident stroke in Chinese adults with chronic diseases

Analysis code for:

> **Triglyceride-Glucose Index and Risk of Incident Stroke in Middle-Aged and Older
> Chinese Adults with Chronic Diseases: A Prospective Cohort Study from the China
> Health and Retirement Longitudinal Study**

Corresponding author: Jian Li, Department of Neurology, Beijing Chao Yang Hospital,
Capital Medical University, Beijing, China.

---

## Data availability

This study is a secondary analysis of the **China Health and Retirement Longitudinal
Study (CHARLS)**, a publicly available, de-identified dataset.

- CHARLS website: <http://charls.pku.edu.cn/en/>
- Data access: CHARLS data are released to registered users free of charge upon
  application. Users must agree to the CHARLS terms of use.
- **Individual-level CHARLS data are not redistributed in this repository**, in
  accordance with the CHARLS data use agreement.

To reproduce the analysis, download the 2011 baseline and 2013-2020 follow-up
waves from the CHARLS website, build the extract described in
`R/01_prepare_data.R`, and place it at `data/CHARLS_2011.csv`.

---

## Repository structure

```
.
├── README.md
├── R/
│   ├── 01_prepare_data.R            # cohort construction, TyG derivation
│   ├── 02_main_analysis.R           # Cox models, tertiles, trend test
│   ├── 03_threshold_analysis.R      # RCS and two-piecewise Cox regression
│   ├── 04_subgroup_and_sensitivity.R
│   └── 05_figures.R                 # Figures 2-4 and supplementary figure
├── data/                            # not tracked (see .gitignore)
├── results/                         # machine-readable outputs
└── figures/                         # generated figures
```

---

## Software

- R version 4.3.1 or later
- Required packages:

```r
install.packages(c(
  "dplyr", "readr", "survival", "rms",
  "survminer", "ggplot2"
))
```

---

## How to run

```r
setwd("path/to/repository")

# 1. Build the analysis cohort (requires CHARLS data; see above)
source("R/01_prepare_data.R")

# 2. Primary models -> Table 2
source("R/02_main_analysis.R")

# 3. Threshold analysis -> Table 3
source("R/03_threshold_analysis.R")

# 4. Subgroups and sensitivity analyses
source("R/04_subgroup_and_sensitivity.R")

# 5. Figures
source("R/05_figures.R")
```

All outputs are written to `results/` and `figures/`.

---

## Key definitions

### TyG index

```
TyG = ln[ fasting triglycerides (mg/dL) × fasting glucose (mg/dL) / 2 ]
```

Measured once, at baseline (CHARLS Wave 1, 2011-2012).

### Cohort

Adults aged >= 45 years with at least one self-reported physician-diagnosed chronic
condition and no history of stroke at baseline.

| Step | N |
|---|---|
| Assessed at baseline | 17,708 |
| Excluded: age < 45 or missing | 648 |
| Excluded: missing chronic disease history | 5,500 |
| Excluded: prevalent stroke at baseline | 167 |
| Excluded: history of cancer | 158 |
| Excluded: TyG unavailable | 3,567 |
| Excluded: other missing data / loss to follow-up | 452 |
| **Eligible cohort** | **7,216** |
| Excluded: missing adjustment covariates | 1,181 |
| **Complete-case analytic sample** | **6,035** (677 incident strokes) |

### Outcomes

Incident stroke: first self-reported physician-diagnosed stroke during Waves 2-5
(2013-2020). Follow-up time in years.

### Covariates in the fully adjusted model (Model 2)

Age, sex, area of residence, alcohol consumption, smoking status, annual income
category, BMI group, C-reactive protein, platelet count, systolic blood pressure,
diastolic blood pressure.

---

## Important note on the threshold analysis

Two separate tests are reported in the manuscript and **must not be conflated**:

1. **RCS log-likelihood ratio test for non-linearity** — compares a 5-knot
   restricted cubic spline against the linear model. In this dataset this test
   was **not statistically significant (P = 0.068)**.

2. **Two-piecewise model log-likelihood ratio test** — compares a broken-line
   (two-segment) model against the linear model. This test **was significant**
   (P = 0.018 in R).

These tests answer different questions and can disagree. The two-piecewise model
identifies a better-fitting change in slope, but this does not by itself establish
a non-linear dose-response relationship.

**The low-TyG tail is sparse.** In the complete-case sample only 5 participants
had TyG < 7.0 (0 events) and 66 had TyG < 7.5 (1 event). The steep slope below
the inflection point is therefore estimated from very few observations and wide
confidence intervals. The threshold finding should be treated as
hypothesis-generating and requires external validation.

---

## Reproducibility notes

- R was used for the analyses reported in the manuscript (via EmpowerStats).
  This repository provides a clean, self-contained re-implementation in base R
  and the `survival` / `rms` packages.
- Small numerical differences between this code and the manuscript values may
  arise from differences in the tertile assignment routine, the knot placement
  algorithm, and optimiser defaults. Direction, magnitude and significance of all
  reported associations are unchanged.
- The inflection point search is performed over the grid TyG = 7.5 to 10.0 in
  steps of 0.1, restricted to the interior of the distribution so that neither
  segment is estimated from a handful of observations alone.

---

## Citation

If you use this code, please cite the paper above and this repository:

```
Li J, Li Q, Zhou L. Analysis code for "Triglyceride-Glucose Index and Risk of
Incident Stroke in Middle-Aged and Older Chinese Adults with Chronic Diseases".
Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX
```

## Licence

MIT
