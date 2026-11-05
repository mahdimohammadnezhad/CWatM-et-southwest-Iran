"""Projected changes in basin-mean ET relative to the 1995-2019 retrospective baseline.

Each GCM-driven monthly series is scaled to the retrospective run over 2015-2019 before the statistics.
Outputs (results/):
  changes_by_gcm.csv        % change of Mean/Max/Min/Std for each GCM, SSP, period and time unit
  table2_changes.csv        mean of the five GCM-specific changes with their range (Table 2)
  tableS1_overlap.csv       seasonal-cycle comparison over the 2015-2019 overlap (Table S1)
  tableS2_statistics.csv    ensemble-mean changes of all statistics (Table S2)
  monthly_signals.csv       absolute changes, sign agreement and seasonal shares (Fig. 8)
  annual_series.csv         annual ET of all runs and the products (Fig. 7)
  baseline_statistics.csv
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from etcompare import cwatm, projections  # noqa: E402
from etcompare.config import GCMS, MONTHS, RESULTS, SEASONS, SSPS  # noqa: E402

runs = cwatm.load_runs()
changes, base, scale = projections.change_statistics(runs)
summary = projections.summarize(changes)
changes.to_csv(RESULTS / "changes_by_gcm.csv", index=False, float_format="%.4f")
base.to_csv(RESULTS / "baseline_statistics.csv", index=False, float_format="%.5f")

# Table 2: annual and seasonal means, standard deviation of monthly values
t2 = summary[((summary["stat"] == "Mean") & summary["unit"].isin(["Annual", *SEASONS])) |
             ((summary["stat"] == "Std") & (summary["unit"] == "Monthly"))].copy()
t2["column"] = t2["unit"].where(t2["stat"] == "Mean", "Monthly Std")
t2["value"] = [f"{m:.1f} [{a:.1f}, {b:.1f}]" for m, a, b in zip(t2["mean"], t2["min"], t2["max"])]
t2.pivot(index=["period", "ssp"], columns="column", values="value")[["Annual", *SEASONS, "Monthly Std"]].to_csv(
    RESULTS / "table2_changes.csv")

# Table S2: ensemble-mean changes of all statistics
order = ["Monthly", "Seasonal", "Annual", *SEASONS, *MONTHS]
s2 = summary.pivot_table(index=["period", "unit"], columns=["ssp", "stat"], values="mean")
s2 = s2.reindex(pd.MultiIndex.from_product([sorted(summary["period"].unique()), order]))
s2.to_csv(RESULTS / "tableS2_statistics.csv", float_format="%.2f")

projections.overlap_shape(runs).to_csv(RESULTS / "tableS1_overlap.csv", index=False, float_format="%.3f")

signals = projections.monthly_signals(summary, base)
signals.to_csv(RESULTS / "monthly_signals.csv", index=False, float_format="%.4f")

# annual series for Fig. 7
hist = cwatm.run_series(runs, "historical")
prod = pd.read_csv(RESULTS / "product_series_monthly.csv", index_col=0)
prod.index = pd.PeriodIndex(prod.index, freq="M")
adjust = pd.read_csv(RESULTS / "mean_adjustment.csv", index_col="dataset")["offset_mm_d"]
annual = {"CWatM": projections.annual_means(hist)}
for name in ["GLDAS", "SMAP L4"]:
    p = prod[name]
    p = p[(p.index >= pd.Period("2000-01", "M")) & (p.index <= pd.Period("2019-12", "M"))]
    annual[name] = projections.annual_means(p)
    annual[f"{name} adjusted"] = projections.annual_means(p - adjust[name])
for s in SSPS:
    scaled = {}
    for g in GCMS:
        x, _ = projections.scale_to_overlap(hist, cwatm.run_series(runs, f"{g}_{s}"))
        scaled[g] = projections.annual_means(x)
        annual[f"{g}_{s}"] = scaled[g]
    annual[f"ensemble_{s}"] = pd.concat(scaled, axis=1).mean(axis=1)
pd.DataFrame(annual).rename_axis("year").to_csv(RESULTS / "annual_series.csv", float_format="%.5f")

pd.set_option("display.width", 200)
print("scale factors:", {k: round(v, 2) for k, v in scale.items()})
print(pd.read_csv(RESULTS / "table2_changes.csv").to_string(index=False))
