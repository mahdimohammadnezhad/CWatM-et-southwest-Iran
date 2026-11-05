"""Agreement between basin-mean monthly CWatM ET and the reference datasets (2000-2019).

Outputs (results/):
  table1_agreement.csv        N, r, RMSE, PBIAS, NSE, KGE and upper-decile PBIAS (Table 1)
  pbias_by_class.csv          PBIAS and mean difference by reference quartile and calendar quarter
  monthly_climatology.csv     mean and standard deviation of each calendar month over the paired months
  mean_adjustment.csv         additive offsets of GLDAS-Noah and SMAP L4 (inter-annual comparison only)
  weighting_sensitivity.csv   Table 1 statistics with polygon-fraction weights
  product_series_monthly.csv  product series on calendar months
  cwatm_basin_monthly.csv     retrospective basin-mean series (mask and polygon weights)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from etcompare import cwatm, metrics, products  # noqa: E402
from etcompare.config import EVALUATION, RESULTS  # noqa: E402

RESULTS.mkdir(exist_ok=True)
start, end = (pd.Period(x, "M") for x in EVALUATION)

sim = cwatm.historical_basin_mean("w")
sim_poly = cwatm.historical_basin_mean("w_poly")
prod = products.load_series()


def paired(ref, s):
    common = ref.dropna().index.intersection(s.index)
    common = common[(common >= start) & (common <= end)]
    return s[common].to_numpy(), ref[common].to_numpy(), common


table, classes, clim, sens, adjust = [], [], [], [], []
for name in products.PRODUCTS:
    s, o, months = paired(prod[name], sim)
    row = {"dataset": name, "start": str(months.min()), "end": str(months.max()), **metrics.scores(s, o)}
    row["PBIAS_upper_decile"], row["N_upper_decile"] = metrics.upper_decile_pbias(s, o)
    row["mean_ref"], row["mean_cwatm"] = o.mean(), s.mean()
    table.append(row)
    classes.append(metrics.by_class(s, o, months.month).assign(dataset=name))
    c = pd.DataFrame({"month": months.month, "ref": o, "cwatm": s}).groupby("month").agg(["mean", "std"])
    c.columns = ["ref_mean", "ref_sd", "cwatm_mean", "cwatm_sd"]
    clim.append(c.reset_index().assign(dataset=name))
    sp, op, _ = paired(prod[name], sim_poly)
    sens.append({"dataset": name, **{k: v for k, v in metrics.scores(sp, op).items() if k in ("r", "RMSE", "PBIAS", "NSE", "KGE")}})
    if name in ("GLDAS", "SMAP L4"):
        adjust.append({"dataset": name, "offset_mm_d": float(np.mean(o - s)), "N": len(o),
                       "start": str(months.min()), "end": str(months.max())})

t1 = pd.DataFrame(table).set_index("dataset")
t1.to_csv(RESULTS / "table1_agreement.csv", float_format="%.4f")
pd.concat(classes).to_csv(RESULTS / "pbias_by_class.csv", index=False, float_format="%.4f")
pd.concat(clim).to_csv(RESULTS / "monthly_climatology.csv", index=False, float_format="%.4f")
pd.DataFrame(adjust).to_csv(RESULTS / "mean_adjustment.csv", index=False, float_format="%.5f")
pd.DataFrame(sens).to_csv(RESULTS / "weighting_sensitivity.csv", index=False, float_format="%.4f")
prod.to_csv(RESULTS / "product_series_monthly.csv", float_format="%.5f")
pd.DataFrame({"CWatM": sim, "CWatM_polygon": sim_poly}).to_csv(RESULTS / "cwatm_basin_monthly.csv", float_format="%.5f")

pd.set_option("display.width", 200)
print(t1[["start", "end", "N", "r", "RMSE", "PBIAS", "NSE", "KGE", "PBIAS_upper_decile", "N_upper_decile"]].round(2))
window = slice(str(start), str(end))
print(f"\nCWatM 2000-2019 mean: mask {sim[window].mean():.3f} mm/d, polygon weights {sim_poly[window].mean():.3f} mm/d")
print("mean adjustment:", {a["dataset"]: round(a["offset_mm_d"], 3) for a in adjust})
