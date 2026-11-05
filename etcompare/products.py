"""Reference ET products: basin-mean monthly series and mean-ET maps exported from Google Earth Engine."""
import numpy as np
import pandas as pd
import rasterio

from .config import DATA

PRODUCTS = ["GLDAS", "TerraClimate", "MODIS 1 km", "MODIS 500 m", "PML-V2", "SMAP L4"]
LABELS = {"GLDAS": "GLDAS-Noah"}
SERIES_COLUMNS = {"GLDAS": "GLDAS", "TerraClimate": "TerraClimate", "MODIS_1km": "MODIS 1 km",
                  "MODIS_500m": "MODIS 500 m", "PML_V2": "PML-V2", "SMAP_L4": "SMAP L4"}
MAP_FILES = {"GLDAS": "GLDAS", "TerraClimate": "TerraClimate", "MODIS_1km": "MODIS 1 km",
             "MODIS_500m": "MODIS 500 m", "PML": "PML-V2", "SMAP": "SMAP L4"}


def load_series():
    """Monthly basin-mean ET (mm/d) of each product, indexed by calendar month.

    Values are the GEE chart exports (one value per monthly window). Windows start at the first
    image of each collection and advance by one month, so the PML-V2 windows start on the 26th and the
    SMAP L4 windows on the last day of a month. Each value is assigned to the month that holds most of
    its window (export date in Iran local time plus 15 days). TerraClimate was exported as monthly
    total / 30 and is rescaled to the number of days in each month.
    """
    raw = pd.read_csv(DATA / "products" / "et_products_basin_monthly.csv")
    t = pd.to_datetime(raw["date"], utc=True).dt.tz_convert("Asia/Tehran").dt.tz_localize(None)
    out = {}
    for col, name in SERIES_COLUMNS.items():
        s = pd.Series(raw[col].to_numpy(float), index=t).dropna()
        m = s.groupby((s.index + pd.Timedelta(days=15)).to_period("M")).mean()
        if name == "TerraClimate":
            m = m * 30 / m.index.days_in_month
        out[name] = m
    return pd.DataFrame(out).sort_index()


def load_maps(transform, shape):
    """Mean-ET maps (mm/d) placed on the analysis grid (files data/products/maps/ET_mean_<name>.tif)."""
    out = {}
    for key, name in MAP_FILES.items():
        f = DATA / "products" / "maps" / f"ET_mean_{key}.tif"
        if not f.exists():
            continue
        with rasterio.open(f) as src:
            if abs(src.transform.a - transform.a) > 1e-9 or abs(src.transform.e - transform.e) > 1e-9:
                raise ValueError(f"{f.name}: resolution differs from the analysis grid")
            a = src.read(1).astype(float)
            if src.nodata is not None:
                a[a == src.nodata] = np.nan
            r0 = round((src.transform.f - transform.f) / transform.e)
            c0 = round((src.transform.c - transform.c) / transform.a)
        a[~np.isfinite(a)] = np.nan
        grid = np.full(shape, np.nan)
        rs, cs = max(r0, 0), max(c0, 0)
        re, ce = min(r0 + a.shape[0], shape[0]), min(c0 + a.shape[1], shape[1])
        grid[rs:re, cs:ce] = a[rs - r0:re - r0, cs - c0:ce - c0]
        out[name] = grid
    return {n: out[n] for n in PRODUCTS if n in out}


def label(name):
    return LABELS.get(name, name)
