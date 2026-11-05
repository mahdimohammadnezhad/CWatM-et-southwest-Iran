"""Projected changes in basin-mean ET relative to the retrospective baseline."""
import numpy as np
import pandas as pd

from .config import BASELINE, GCMS, MONTHS, OVERLAP, PERIODS, SEASONS, SSPS
from .cwatm import run_series

STATS = ["Mean", "Max", "Min", "Std"]


def window(s, y0, y1):
    return s[(s.index.year >= y0) & (s.index.year <= y1)]


def scale_to_overlap(hist, gcm):
    """Multiply a GCM-driven series by mean(retrospective)/mean(GCM-driven) over the overlap years."""
    f = window(hist, *OVERLAP).mean() / window(gcm, *OVERLAP).mean()
    return gcm * f, f


def _st(x):
    x = np.asarray(x, float)
    return {"Mean": x.mean(), "Max": x.max(), "Min": x.min(), "Std": x.std(ddof=1)}


def period_stats(s):
    """Statistics of a monthly series: all months, seasonal and annual means, each season and month.

    Seasons are taken within each calendar year (December of the same year is part of winter).
    """
    df = pd.DataFrame({"v": s.to_numpy(), "y": s.index.year, "m": s.index.month})
    out = {"Monthly": _st(df["v"])}
    annual = df.groupby("y")["v"].mean()
    seasonal, per_season = [], {k: [] for k in SEASONS}
    for _, g in df.groupby("y"):
        for name, months in SEASONS.items():
            v = g.loc[g["m"].isin(months), "v"]
            if len(v) == 3:
                seasonal.append(v.mean())
                per_season[name].append(v.mean())
    out["Seasonal"] = _st(seasonal)
    out["Annual"] = _st(annual)
    for name in SEASONS:
        out[name] = _st(per_season[name])
    for m in range(1, 13):
        out[MONTHS[m - 1]] = _st(df.loc[df["m"] == m, "v"])
    return out


def pct_change(future, base):
    return {k: {st: 100 * (future[k][st] - base[k][st]) / base[k][st] for st in STATS} for k in base}


def change_statistics(runs):
    """Percentage change of each statistic for every GCM, SSP and period (long table)."""
    hist = run_series(runs, "historical")
    base = period_stats(window(hist, *BASELINE))
    rows, scale = [], {}
    for s in SSPS:
        for g in GCMS:
            x, f = scale_to_overlap(hist, run_series(runs, f"{g}_{s}"))
            scale[f"{g}_{s}"] = f
            for p, (y0, y1) in PERIODS.items():
                ch = pct_change(period_stats(window(x, y0, y1)), base)
                for unit, d in ch.items():
                    for st, v in d.items():
                        rows.append({"period": p, "ssp": SSPS[s], "gcm": g, "unit": unit, "stat": st, "change_pct": v})
    base_df = pd.DataFrame([{"unit": u, "stat": st, "baseline": v} for u, d in base.items() for st, v in d.items()])
    return pd.DataFrame(rows), base_df, scale


def summarize(changes):
    """Mean of the GCM-specific changes with their range."""
    g = changes.groupby(["period", "ssp", "unit", "stat"])["change_pct"]
    return g.agg(mean="mean", min="min", max="max", n_increase=lambda v: int((v > 0).sum())).reset_index()


def overlap_shape(runs):
    """Seasonal-cycle difference between each scaled GCM-driven run and the retrospective run (overlap years)."""
    def clim(x):
        w = window(x, *OVERLAP)
        return w.groupby(w.index.month).mean().to_numpy()

    hist = clim(run_series(runs, "historical"))
    seasons = {"DJF": (12, 1, 2), "MAM": (3, 4, 5), "JJA": (6, 7, 8), "SON": (9, 10, 11)}

    def compare(c):
        f = hist.mean() / c.mean()
        adj = c * f
        absdiff = 100 * np.abs(adj - hist) / hist.mean()
        k = int(np.argmax(absdiff))
        row = {"scale_factor": f, "MAD_pct": absdiff.mean(), "MAX_pct": absdiff[k], "MAX_month": MONTHS[k],
               "r": np.corrcoef(hist, c)[0, 1]}
        for name, ms in seasons.items():
            idx = [m - 1 for m in ms]
            row[f"{name}_pct"] = 100 * (adj[idx].mean() - hist[idx].mean()) / hist[idx].mean()
        return row, adj

    rows, scaled = [], []
    for g in GCMS:
        for s in SSPS:
            row, adj = compare(clim(run_series(runs, f"{g}_{s}")))
            rows.append({"GCM": g, "SSP": SSPS[s], **row})
            scaled.append(adj)
    ens, _ = compare(np.mean(scaled, axis=0))
    rows.append({"GCM": "Ensemble mean", "SSP": "all", **ens})
    return pd.DataFrame(rows)


def monthly_signals(summary, base):
    """Absolute changes, sign agreement and seasonal shares from the ensemble-mean changes of the means."""
    m = summary[summary["stat"] == "Mean"].merge(base[base["stat"] == "Mean"][["unit", "baseline"]], on="unit")
    m["change_mm_d"] = m["mean"] / 100 * m["baseline"]
    m["n_agree"] = np.where(m["mean"] > 0, m["n_increase"], len(GCMS) - m["n_increase"])
    seas = m[m["unit"].isin(SEASONS)].copy()
    seas["future"] = seas["baseline"] * (1 + seas["mean"] / 100)
    tot_b = seas.groupby(["period", "ssp"])["baseline"].transform("sum")
    tot_f = seas.groupby(["period", "ssp"])["future"].transform("sum")
    seas["share_base"] = 100 * seas["baseline"] / tot_b
    seas["share_future"] = 100 * seas["future"] / tot_f
    out = m.merge(seas[["period", "ssp", "unit", "share_base", "share_future"]], on=["period", "ssp", "unit"], how="left")
    return out.drop(columns="stat")


def annual_means(s):
    """Mean of the monthly values of each complete calendar year."""
    s = s.dropna()
    g = s.groupby(s.index.year)
    return g.mean()[g.size() == 12]
