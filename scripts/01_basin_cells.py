"""Cell weights for basin-mean CWatM series.

w      : area of the 61 cells of the study-basin mask (data/cwatm/basin_mask_cells.csv); used throughout
w_poly : area of each cell inside the basin polygon (sensitivity test)
Cell areas and the share of each cell inside the basin polygon are computed in an equal-area projection.
Output : data/cwatm/basin_cells.csv
"""
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from etcompare.config import BASIN_SHP, DATA  # noqa: E402

EQUAL_AREA = "EPSG:6933"

basin = gpd.read_file(BASIN_SHP).to_crs(4326).union_all()
basin_ea = gpd.GeoSeries([basin], crs=4326).to_crs(EQUAL_AREA).iloc[0]
x0, y0, x1, y1 = basin.bounds
lats = np.arange(np.floor(y0 * 2) / 2 + 0.25, y1 + 0.5, 0.5)
lons = np.arange(np.floor(x0 * 2) / 2 + 0.25, x1 + 0.5, 0.5)
mask = pd.read_csv(DATA / "cwatm" / "basin_mask_cells.csv").round(2)
cells = pd.DataFrame([(round(la, 2), round(lo, 2)) for la in lats for lo in lons], columns=["lat", "lon"])
cells = pd.concat([cells, mask]).drop_duplicates().reset_index(drop=True)
geoms = [box(lo - 0.25, la - 0.25, lo + 0.25, la + 0.25) for la, lo in zip(cells["lat"], cells["lon"])]
g = gpd.GeoDataFrame(cells, geometry=geoms, crs=4326).to_crs(EQUAL_AREA)
g["cell_km2"] = g.area / 1e6
g["frac"] = g.intersection(basin_ea).area / g.area

in_mask = set(zip(mask["lat"], mask["lon"]))
g["in_mask"] = [(a, b) in in_mask for a, b in zip(g["lat"], g["lon"])]
out = pd.DataFrame(g.drop(columns="geometry"))
out = out[(out["frac"] > 0) | out["in_mask"]].sort_values(["lat", "lon"], ascending=[False, True])
out["w"] = out["cell_km2"].where(out["in_mask"], 0.0)
out["w_poly"] = out["frac"] * out["cell_km2"]
out[["lat", "lon", "cell_km2", "frac", "in_mask", "w", "w_poly"]].to_csv(
    DATA / "cwatm" / "basin_cells.csv", index=False, float_format="%.6f")
print(f"basin area {basin_ea.area / 1e6:,.0f} km2; cells touching the basin {(out['frac'] > 0).sum()}; "
      f"mask cells {out['in_mask'].sum()}")
