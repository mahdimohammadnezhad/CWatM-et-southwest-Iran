"""Shared figure style and map helpers."""
import geopandas as gpd
import matplotlib as mpl
import numpy as np
from rasterio.features import rasterize

from .config import BASIN_SHP, DATA, SUB_BASINS

INK, MUTED, GRID = "#1f1f1f", "#5f5f5f", "#e6e6e6"
WIDTH = 7.09  # double-column width (inches)
STYLE = {
    "CWatM": dict(color="#000000", marker="o", ls="-"),
    "GLDAS": dict(color="#0072B2", marker="s", ls="-"),
    "TerraClimate": dict(color="#E69F00", marker="^", ls="-"),
    "MODIS 1 km": dict(color="#009E73", marker="v", ls="--"),
    "MODIS 500 m": dict(color="#009E73", marker="D", ls="-"),
    "PML-V2": dict(color="#CC79A7", marker="P", ls="-"),
    "SMAP L4": dict(color="#D55E00", marker="X", ls="-"),
}
MAP_VIEW = [45.2, 52.2, 30.05, 35.95]


def setup():
    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8,
        "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
        "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "savefig.dpi": 400,
        "savefig.bbox": "tight", "savefig.pad_inches": 0.03, "legend.frameon": False,
    })


def panel_label(ax, s, x=-0.02, y=1.02):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="right", color=INK)


def main_basin():
    """Largest connected part of the basin polygon (the part drawn on maps)."""
    u = gpd.read_file(BASIN_SHP).to_crs(4326).union_all()
    parts = list(u.geoms) if u.geom_type == "MultiPolygon" else [u]
    return max(parts, key=lambda g: g.area)


def map_mask(transform, shape):
    return rasterize([(main_basin(), 1)], out_shape=shape, transform=transform, fill=0, dtype="uint8") == 1


def basin_outlines():
    main = main_basin().buffer(1e-6)
    out = {}
    for name, f in SUB_BASINS.items():
        g = gpd.read_file(DATA / "basin" / f).to_crs(4326)
        g["geometry"] = g.geometry.intersection(main)
        out[name] = g[~g.geometry.is_empty]
    return out


def draw_outlines(ax, outlines, lw=0.6, color=INK):
    for g in outlines.values():
        g.boundary.plot(ax=ax, color=color, linewidth=lw)


def map_axes(ax, extent=MAP_VIEW):
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_aspect(1 / np.cos(np.deg2rad(33.0)))
    ax.set_xticks([46, 48, 50, 52])
    ax.set_yticks([31, 33, 35])
    ax.set_xticklabels([f"{v}°E" for v in [46, 48, 50, 52]])
    ax.set_yticklabels([f"{v}°N" for v in [31, 33, 35]])
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_linewidth(0.5)
        s.set_color("#8c8c8c")
    ax.tick_params(length=2)
