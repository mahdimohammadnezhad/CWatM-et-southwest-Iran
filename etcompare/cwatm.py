"""Reading CWatM NetCDF output and computing basin-mean evapotranspiration."""
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio

from .config import DATA, ET_COMPONENTS

DAY0 = pd.Timestamp("1901-01-01")


def read_monthly(path):
    """Monthly CWatM map stack as (values [time, lat, lon] in mm/d, months, lat, lon).

    CWatM writes depths in m/day; values above 1e10 are fill values.
    """
    with rasterio.open(Path(path)) as src:
        a = src.read().astype(np.float64)
        days = [float(src.tags(b)["NETCDF_DIM_time"]) for b in range(1, src.count + 1)]
        t = src.transform
        lat = t.f + t.e * (np.arange(src.height) + 0.5)
        lon = t.c + t.a * (np.arange(src.width) + 0.5)
    a[~np.isfinite(a) | (np.abs(a) > 1e10)] = np.nan
    months = (DAY0 + pd.to_timedelta(days, unit="D")).to_period("M")
    return a * 1000.0, months, lat, lon


def read_run(run_dir):
    """ET components of one run: dict name -> (values, months, lat, lon)."""
    return {c: read_monthly(Path(run_dir) / f"{c}_monthavg.nc") for c in ET_COMPONENTS}


def cell_series(run_dir, cells):
    """Monthly ET components (mm/d) for the given cells (DataFrame with lat, lon)."""
    parts = read_run(run_dir)
    _, months, lat, lon = parts["totalET"]
    rows = []
    for la, lo in zip(cells["lat"], cells["lon"]):
        i, j = int(np.argmin(np.abs(lat - la))), int(np.argmin(np.abs(lon - lo)))
        rec = {"ym": months.astype(str), "lat": round(la, 2), "lon": round(lo, 2)}
        for c, key in zip(ET_COMPONENTS, ["totalET", "waterbody", "channel"]):
            rec[key] = parts[c][0][:, i, j]
        rows.append(pd.DataFrame(rec))
    out = pd.concat(rows, ignore_index=True)
    out = out[out["totalET"].notna()]
    return out.fillna({"waterbody": 0.0, "channel": 0.0})


def basin_series(run_dir, weights):
    """Area-weighted basin-mean monthly ET (totalET + lake/reservoir + channel evaporation), mm/d.

    ``weights`` has columns lat, lon, w; cells with w = 0 are ignored.
    """
    cells = weights[weights["w"] > 0]
    c = cell_series(run_dir, cells).merge(cells[["lat", "lon", "w"]], on=["lat", "lon"])
    c["ET"] = c["totalET"] + c["waterbody"] + c["channel"]
    s = (c["ET"] * c["w"]).groupby(c["ym"]).sum() / c["w"].groupby(c["ym"]).sum()
    s.index = pd.PeriodIndex(s.index, freq="M")
    return s.rename("ET")


def load_weights(column="w"):
    """Cell weights; 'w' = 61-cell study-basin mask (cell area), 'w_poly' = area inside the basin polygon."""
    w = pd.read_csv(DATA / "cwatm" / "basin_cells.csv")
    return w.assign(w=w[column])


def load_cells_historical():
    """Monthly ET components of the retrospective run for all cells touching the basin."""
    c = pd.read_csv(DATA / "cwatm" / "cells_historical_monthly.csv")
    c["ET"] = c["totalET"] + c["waterbody"] + c["channel"]
    return c


def historical_basin_mean(column="w"):
    """Basin-mean monthly ET of the retrospective run from the stored cell series."""
    w = load_weights(column)
    c = load_cells_historical().merge(w[["lat", "lon", "w"]], on=["lat", "lon"])
    c = c[c["w"] > 0]
    s = (c["ET"] * c["w"]).groupby(c["ym"]).sum() / c["w"].groupby(c["ym"]).sum()
    s.index = pd.PeriodIndex(s.index, freq="M")
    return s.rename("CWatM")


def load_runs():
    """Basin-mean monthly ET of all runs (long table: run, ym, ET)."""
    r = pd.read_csv(DATA / "cwatm" / "basin_monthly_et.csv")
    r["ym"] = pd.PeriodIndex(r["ym"], freq="M")
    return r


def run_series(runs, name):
    return runs[runs["run"] == name].set_index("ym")["ET"]
