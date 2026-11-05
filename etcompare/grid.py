"""Common 0.01-degree analysis grid with elevation, land cover, Koppen-Geiger class and basin mask."""
import gzip
import json
import urllib.request

import geopandas as gpd
import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.features import rasterize
from rasterio.transform import from_origin
from rasterio.warp import reproject, transform_bounds
from rasterio.windows import from_bounds

from .config import BASIN_SHP, DATA, KOPPEN_SHP

# grid aligned with the 0.5-degree CWatM grid
LON0, LAT0, LON1, LAT1 = 44.5, 37.0, 52.5, 30.0

LULC = {1: "Water", 2: "Trees", 4: "Flooded vegetation", 5: "Crops", 7: "Built area", 8: "Bare ground",
        9: "Snow/ice", 10: "Clouds", 11: "Rangeland"}
LULC_URL = "https://ai4edataeuwest.blob.core.windows.net/io-lulc/io-annual-lulc-v02/{tile}_20170101-20180101.tif"
LULC_TOKEN = "https://planetarycomputer.microsoft.com/api/sas/v1/token/io-lulc-annual-v02"
LULC_TILES = ["38R", "38S", "39R", "39S"]
SRTM_URL = "https://s3.amazonaws.com/elevation-tiles-prod/skadi/{lat}/{name}.hgt.gz"
KOPPEN = {22: "BWh", 27: "BSh", 26: "BSk", 34: "Csa", 45: "Dsa", 46: "Dsb"}   # Kottek et al. (2006) codes


def make_grid(res):
    w, h = round((LON1 - LON0) / res), round((LAT0 - LAT1) / res)
    return from_origin(LON0, LAT0, res, res), (h, w)


def basin_geometry():
    return gpd.read_file(BASIN_SHP).to_crs(4326).union_all()


def koppen_polygons():
    # the shapefile has no .prj; coordinates are Web Mercator
    return gpd.read_file(KOPPEN_SHP).set_crs(3857, allow_override=True).to_crs(4326)


def srtm_tiles(tile_dir):
    """Download the SRTM 1-arcsecond tiles covering the basin (AWS Terrain Tiles, skadi format)."""
    tile_dir.mkdir(parents=True, exist_ok=True)
    x0, y0, x1, y1 = basin_geometry().bounds
    files = []
    for la in range(int(np.floor(y0)), int(np.floor(y1)) + 1):
        for lo in range(int(np.floor(x0)), int(np.floor(x1)) + 1):
            name = f"N{la:02d}E{lo:03d}"
            f = tile_dir / f"{name}.hgt.gz"
            if not f.exists():
                urllib.request.urlretrieve(SRTM_URL.format(lat=name[:3], name=name), f)
            files.append(f)
    return files


def elevation(transform, shape, tile_dir):
    out = np.full(shape, np.nan, np.float32)
    for f in srtm_tiles(tile_dir):
        la, lo = int(f.name[1:3]), int(f.name[4:7])
        a = np.frombuffer(gzip.decompress(f.read_bytes()), ">i2").reshape(3601, 3601).astype(np.float32)
        a[a == -32768] = np.nan
        tile = np.full(shape, np.nan, np.float32)
        src_tr = from_origin(lo - 0.5 / 3600, la + 1 + 0.5 / 3600, 1 / 3600, 1 / 3600)
        reproject(a, tile, src_transform=src_tr, src_crs="EPSG:4326", src_nodata=np.nan, dst_nodata=np.nan,
                  dst_transform=transform, dst_crs="EPSG:4326", resampling=Resampling.average)
        out = np.where(np.isnan(out), tile, out)
    return out


def land_cover(transform, shape, transform_map, shape_map, factor=8):
    """Class fractions on the analysis grid and the dominant class on the map grid.

    The 10 m Esri/Impact Observatory map is read at 80 m (``factor`` = 8) from the Planetary Computer.
    """
    token = json.load(urllib.request.urlopen(LULC_TOKEN, timeout=60))["token"]
    num = {c: np.zeros(shape, np.float32) for c in LULC}
    den = np.zeros(shape, np.float32)
    dom = np.zeros(shape_map, np.uint8)
    bounds = basin_geometry().buffer(0.05).bounds
    for tile in LULC_TILES:
        with rasterio.open(LULC_URL.format(tile=tile) + "?" + token) as src:
            l, b, r, t = transform_bounds("EPSG:4326", src.crs, *bounds)
            l, b = max(l, src.bounds.left), max(b, src.bounds.bottom)
            r, t = min(r, src.bounds.right), min(t, src.bounds.top)
            if l >= r or b >= t:
                continue
            win = from_bounds(l, b, r, t, src.transform).round_offsets().round_lengths()
            h, w = int(win.height // factor), int(win.width // factor)
            a = src.read(1, window=win, out_shape=(h, w), resampling=Resampling.nearest)
            tr = src.window_transform(win) * rasterio.Affine.scale(win.width / w, win.height / h)
            for code, target in [(None, den)] + [(c, num[c]) for c in LULC]:
                d = np.zeros(shape, np.float32)
                layer = (a > 0) if code is None else (a == code)
                reproject(layer.astype(np.float32), d, src_transform=tr, src_crs=src.crs, dst_transform=transform,
                          dst_crs="EPSG:4326", resampling=Resampling.average)
                target += d
            m = np.zeros(shape_map, np.uint8)
            reproject(a, m, src_transform=tr, src_crs=src.crs, src_nodata=0, dst_nodata=0, dst_transform=transform_map,
                      dst_crs="EPSG:4326", resampling=Resampling.mode)
            dom = np.where(dom == 0, m, dom)
    frac = {c: np.where(den > 0, num[c] / np.maximum(den, 1e-6), np.nan) for c in LULC}
    return frac, dom


def burn_basin_and_koppen(transform, shape):
    basin = rasterize([(basin_geometry(), 1)], out_shape=shape, transform=transform, fill=0, dtype="uint8")
    kop = koppen_polygons()
    k = rasterize(((g, int(c)) for g, c in zip(kop.geometry, kop["GRIDCODE"])), out_shape=shape,
                  transform=transform, fill=0, dtype="uint8")
    return basin, k


def load(res="001"):
    """Stored grid: 'env_grid_001.npz' (0.01 degree, analysis) or 'env_grid_00025.npz' (0.0025 degree, maps)."""
    env = np.load(DATA / "grid" / f"env_grid_{res}.npz")
    return env, rasterio.Affine(*env["transform"])


def pixel_centres(transform, shape):
    h, w = shape
    lon = transform.c + transform.a * (np.arange(w) + 0.5)
    lat = transform.f + transform.e * (np.arange(h) + 0.5)
    return np.meshgrid(lon, lat)


def area_weights(transform, shape):
    _, lat = pixel_centres(transform, shape)
    return np.cos(np.deg2rad(lat))


def cwatm_to_grid(cell_values, transform, shape):
    """Assign values of 0.5-degree cells (Series indexed by (lat, lon) of cell centres) to the grid pixels."""
    lon, lat = pixel_centres(transform, shape)
    key_lat = np.round(np.floor(lat / 0.5) * 0.5 + 0.25, 2)
    key_lon = np.round(np.floor(lon / 0.5) * 0.5 + 0.25, 2)
    look = {(round(a, 2), round(b, 2)): v for (a, b), v in cell_values.items()}
    vals = [look.get((a, b), np.nan) for a, b in zip(key_lat.ravel(), key_lon.ravel())]
    return np.array(vals).reshape(shape), key_lat, key_lon


def cwatm_mean_map(cells, start="2000-01", end="2019-12"):
    """Mean ET (mm/d) of each CWatM cell over [start, end] from the stored cell series."""
    c = cells[(cells["ym"] >= start) & (cells["ym"] <= end)]
    return c.groupby(["lat", "lon"])["ET"].mean()


def weighted_mean(values, weights):
    m = np.isfinite(values) & (weights > 0)
    return float(np.sum(values[m] * weights[m]) / np.sum(weights[m])) if m.any() else np.nan
