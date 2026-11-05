"""Inter-annual association of annual ET with annual precipitation and mean air temperature.

Basin means are taken over all grid cells touching the basin polygon (simple mean). Annual ET and P are
sums of daily values (mm/yr), T is the annual mean (degC). ET is CWatM land-surface ET (totalET_daily.nc)
without the level scaling used for the projections.
Retrospective: GSWP3-W5E5 forcing, 1981-2019. Projections: ISIMIP3b forcing of each GCM and SSP, 2015-2100.
Needs CWATM_OUTPUT, FORCING_HIST and FORCING_FUTURE; with --tables-only, Table S6 is recomputed from
results/annual_climate.csv and no raw data are read.
Outputs:
  data/climate/basin_monthly_climate.csv  monthly ET and P (mm/month) and T (degC) of every run
  results/annual_climate.csv              annual ET, P and T of every run (Fig. 9 uses the five-GCM mean)
  results/tableS6_relationships.csv       Pearson, Spearman and linear regression of ET on P and on T (Table S6)
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.features import rasterize
from rasterio.windows import from_bounds
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from etcompare.config import CWATM_OUTPUT, DATA, FORCING_FUTURE, FORCING_HIST, GCMS, RESULTS, SSPS  # noqa: E402
from etcompare.grid import basin_geometry  # noqa: E402

UNIT_SECONDS = {"days": 86400.0, "hours": 3600.0, "minutes": 60.0, "seconds": 1.0}


def basin_daily(path, start, end):
    """Daily basin-mean values of the single variable in a NetCDF file, in its native units."""
    geom = basin_geometry()
    with rasterio.open(path) as src:
        m = re.match(r"(\w+) since (.+)", src.tags().get("time#units", "days since 1901-01-01"), re.I)
        unit, origin = m.group(1).lower(), pd.Timestamp(m.group(2).strip()[:19])
        offsets = [float(src.tags(b)["NETCDF_DIM_time"]) for b in range(1, src.count + 1)]
        time = origin + pd.to_timedelta(np.array(offsets) * UNIT_SECONDS[unit], unit="s")
        keep = np.where((time >= pd.Timestamp(start)) & (time <= pd.Timestamp(end) + pd.Timedelta(hours=23)))[0]
        if len(keep) == 0:
            return pd.Series(dtype=float)
        x0, y0, x1, y1 = geom.bounds
        win = from_bounds(x0 - 1, y0 - 1, x1 + 1, y1 + 1, src.transform).round_offsets().round_lengths()
        a = src.read([int(b) + 1 for b in keep], window=win).astype(float)
        nodata = src.nodata
        tr = src.window_transform(win)
        units = src.tags(1).get("units", "")
    mask = rasterize([(geom, 1)], out_shape=a.shape[1:], transform=tr, fill=0, all_touched=True, dtype="uint8") == 1
    if nodata is not None:
        a[a == nodata] = np.nan
    a[np.abs(a) > 1e10] = np.nan
    s = pd.Series(np.nanmean(a[:, mask], axis=1), index=time[keep].normalize())
    s.attrs["units"] = units
    return s


def to_mm_day(s):
    u = s.attrs.get("units", "").replace(" ", "").lower()
    if "kgm-2s-1" in u or "mm/s" in u or "mms-1" in u:
        return s * 86400.0
    if u == "m":
        return s * 1000.0
    return s


def to_celsius(s):
    return s - 273.15 if s.attrs.get("units", "").strip().lower() in ("k", "kelvin") else s


def read_files(files, start, end, convert):
    parts = [convert(basin_daily(f, start, end)) for f in files]
    s = pd.concat([p for p in parts if len(p)]).sort_index()
    return s[~s.index.duplicated()]


def aggregate(et, pr, tas, freq):
    """Sums of ET and P and mean of T over complete calendar months ('M') or years ('Y').

    A year needs at least 365 days and February at least 28, so that series without leap days are kept.
    """
    df = pd.concat({"ET": et, "P": pr, "T": tas}, axis=1, join="inner")
    g = df.groupby(df.index.to_period(freq))
    out = pd.DataFrame({"ET_mm": g["ET"].sum(), "P_mm": g["P"].sum(), "T_C": g["T"].mean()})
    need = [365 if freq == "Y" else min(p.days_in_month, 28 if p.month == 2 else 31) for p in out.index]
    return out[g.size().to_numpy() >= need]


def annual_table(et, pr, tas):
    out = aggregate(et, pr, tas, "Y")
    out.index = out.index.year
    return out.rename_axis("year")


def monthly_table(et, pr, tas):
    out = aggregate(et, pr, tas, "M")
    out.index = out.index.astype(str)
    return out.rename_axis("ym")


def relationship(x, y):
    r, p = stats.pearsonr(x, y)
    rho, p_rho = stats.spearmanr(x, y)
    fit = stats.linregress(x, y)
    return {"N": len(x), "pearson_r": r, "pearson_p": p, "spearman_rho": rho, "spearman_p": p_rho,
            "slope": fit.slope, "intercept": fit.intercept, "slope_se": fit.stderr, "R2": fit.rvalue ** 2}


def forcing_files(folder, variable, pattern):
    files = sorted(Path(folder).glob(pattern.format(v=variable)))
    return [f for f in files if f"_{variable}_" in f.name]  # 'pr' must not match 'prsn'


def extract():
    """Read the raw data; write the monthly and annual tables and return the annual one."""
    annual, monthly = [], []

    def add_run(name, et, pr, tas):
        annual.append(annual_table(et, pr, tas).assign(run=name))
        monthly.append(monthly_table(et, pr, tas).assign(run=name))
        print(name, "years", len(annual[-1]), "months", len(monthly[-1]), flush=True)

    hist = ("1981-01-01", "2019-12-31")
    et = read_files([CWATM_OUTPUT / "historical" / "totalET_daily.nc"], *hist, to_mm_day)
    pr = read_files(forcing_files(FORCING_HIST, "pr", "gswp3-w5e5_obsclim_{v}_global_daily_*.nc"), *hist, to_mm_day)
    tas = read_files(forcing_files(FORCING_HIST, "tas", "gswp3-w5e5_obsclim_{v}_global_daily_*.nc"), *hist, to_celsius)
    add_run("historical", et, pr, tas)
    future = ("2015-01-01", "2100-12-31")
    pattern = "*_{v}_lon36.0to53.0lat25.0to41.0_daily_*.nc"
    for g in GCMS:
        for s in SSPS:
            folder = FORCING_FUTURE / g / s
            et = read_files([CWATM_OUTPUT / g / s / "totalET_daily.nc"], *future, to_mm_day)
            pr = read_files(forcing_files(folder, "pr", pattern), *future, to_mm_day)
            tas = read_files(forcing_files(folder, "tas", pattern), *future, to_celsius)
            add_run(f"{g}_{s}", et, pr, tas)

    (DATA / "climate").mkdir(exist_ok=True)
    monthly = pd.concat(monthly).reset_index()[["run", "ym", "ET_mm", "P_mm", "T_C"]]
    monthly.to_csv(DATA / "climate" / "basin_monthly_climate.csv", index=False, float_format="%.6f")
    annual = pd.concat(annual).reset_index()
    annual.to_csv(RESULTS / "annual_climate.csv", index=False, float_format="%.6f")
    return annual


# --tables-only: recompute Table S6 from results/annual_climate.csv (no raw data needed)
annual = pd.read_csv(RESULTS / "annual_climate.csv") if "--tables-only" in sys.argv else extract()

rows = []
groups = [("historical", annual[annual["run"] == "historical"])]
groups += [(f"{g}_{s}", annual[annual["run"] == f"{g}_{s}"]) for g in GCMS for s in SSPS]
for s in SSPS:
    ens = annual[annual["run"].isin([f"{g}_{s}" for g in GCMS])].groupby("year")[["ET_mm", "P_mm", "T_C"]].mean()
    groups.append((f"ensemble_{s}", ens))
for name, d in groups:
    for var in ["P_mm", "T_C"]:
        rows.append({"run": name, "driver": var, **relationship(d[var], d["ET_mm"])})
table = pd.DataFrame(rows)
table.to_csv(RESULTS / "tableS6_relationships.csv", index=False, float_format="%.8g")
print(table[table["run"].str.startswith(("historical", "ensemble"))].round(3).to_string(index=False))
