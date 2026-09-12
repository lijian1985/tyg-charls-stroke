"""
verify_analysis.py

Python re-implementation of the main analyses, used to independently
verify the numbers reported in the manuscript. This is a verification
script, not the primary analysis code (which is in R/).

It requires the CHARLS extract at data/CHARLS_2011.csv, where the
columns use the harmonised CHARLS names and, for convenience, the
pre-derived variables used by the original analysis:

    TYG, STROKE_Y, TIME (years),
    AGE, GENDER, RURAL, DRINKING1.NEW, SMOKING1.NEW,
    INCOME_TOTAL.CS, BMI.CS.1, BL_CRP, BL_PLT, SYSTO, DIASTO

If you only have the raw CHARLS files, build this extract with
R/01_prepare_data.R first.
"""

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from scipy.stats import chi2

COV = [
    "AGE", "GENDER", "RURAL",
    "DRINKING1.NEW", "SMOKING1.NEW",
    "INCOME_TOTAL.CS", "BMI.CS.1",
    "BL_CRP", "BL_PLT", "SYSTO", "DIASTO",
]
TIME, EVENT, EXPOSURE = "TIME", "STROKE_Y", "TYG"


def load(path="data/CHARLS_2011.csv"):
    df = pd.read_csv(path, sep="\t", low_memory=False)
    cc = df.dropna(subset=COV + [EXPOSURE, EVENT, TIME]).copy()
    cc["t3"] = pd.qcut(cc[EXPOSURE], 3, labels=["Low", "Middle", "High"])
    return df, cc


def cox(data, cols):
    f = CoxPHFitter()
    f.fit(data[cols + [TIME, EVENT]], duration_col=TIME, event_col=EVENT)
    return f


def hr_line(f, term):
    s = f.summary.loc[term]
    return (f"{s['exp(coef)']:.2f} "
            f"({s['exp(coef) lower 95%']:.2f}-{s['exp(coef) upper 95%']:.2f}), "
            f"P = {s['p']:.4f}")


def main():
    _, cc = load()

    print(f"Eligible cohort        : {len(cc)}")
    print(f"Incident stroke events : {int(cc[EVENT].sum())}")
    print(f"Total person-years     : {cc[TIME].sum():.0f}")
    print()

    print("== Incidence rate per 1,000 person-years ==")
    for g in ["Low", "Middle", "High"]:
        s = cc[cc["t3"] == g]
        print(f"  {g:7s} n={len(s):5d}  events={int(s[EVENT].sum()):4d}  "
              f"rate={1000 * s[EVENT].sum() / s[TIME].sum():5.1f}")
    print()

    print("== HR per 1-unit increase in TyG ==")
    print(f"  Unadjusted : {hr_line(cox(cc, [EXPOSURE]), EXPOSURE)}")
    print(f"  Model 1    : {hr_line(cox(cc, [EXPOSURE, 'AGE', 'GENDER', 'RURAL']), EXPOSURE)}")
    print(f"  Model 2    : {hr_line(cox(cc, [EXPOSURE] + COV), EXPOSURE)}")
    print()

    print("== Low-TyG tail ==")
    for cut in (7.0, 7.5, 8.0, 8.3):
        s = cc[cc[EXPOSURE] < cut]
        print(f"  TyG < {cut:.1f} : n={len(s):5d}  events={int(s[EVENT].sum()):4d}")
    print()

    # Two-piecewise model at the reported inflection point
    k = 8.3
    d = cc.copy()
    d["seg1"] = np.where(d[EXPOSURE] < k, d[EXPOSURE] - k, 0.0)
    d["seg2"] = np.where(d[EXPOSURE] >= k, d[EXPOSURE] - k, 0.0)

    pw = cox(d, ["seg1", "seg2"] + COV)
    lin = cox(cc, [EXPOSURE] + COV)

    print(f"== Two-piecewise Cox regression (inflection at TyG = {k}) ==")
    print(f"  Below threshold : {hr_line(pw, 'seg1')}")
    print(f"  Above threshold : {hr_line(pw, 'seg2')}")

    lr = 2 * (pw.log_likelihood_ - lin.log_likelihood_)
    print(f"  LRT vs linear   : chi2 = {lr:.3f}, P = {chi2.sf(lr, 1):.4f}")

    # Difference in slopes (Wald) - a different test from the LRT above
    v = pw.variance_matrix_.values
    names = list(pw.params_.index)
    i, j = names.index("seg1"), names.index("seg2")
    diff = pw.params_["seg2"] - pw.params_["seg1"]
    se = np.sqrt(v[i, i] + v[j, j] - 2 * v[i, j])
    from scipy.stats import norm
    print(f"  Difference in slopes: z = {diff / se:.3f}, "
          f"P = {2 * norm.sf(abs(diff / se)):.4f}")
    print()

    print("== Sensitivity analyses (HR per 1-unit TyG) ==")
    scenarios = {
        "Primary (complete case)": cc,
        "Follow-up >= 2 years": cc[cc[TIME] >= 2],
        "Excluding TyG < 7.0": cc[cc[EXPOSURE] >= 7.0],
        "Excluding TyG < 7.5": cc[cc[EXPOSURE] >= 7.5],
    }
    for label, data in scenarios.items():
        if len(data) < 100:
            continue
        print(f"  {label:26s} n={len(data):5d}  "
              f"ev={int(data[EVENT].sum()):4d}  "
              f"{hr_line(cox(data, [EXPOSURE] + COV), EXPOSURE)}")


if __name__ == "__main__":
    main()
