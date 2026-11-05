"""Mean ET by elevation band, Koppen-Geiger class and land-cover class on the 0.01-degree grid.

All datasets are compared on the same pixels: those inside the basin and covered by CWatM cells.
CWatM cell means (2000-2019) are assigned to every pixel of the cell; land-cover means are weighted
by class fraction. Also checks the product maps against the mean of their basin series.
Outputs: results/gradients.json, results/cwatm_cells_environment.csv, results/product_map_check.csv
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from etcompare import cwatm, grid, products  # noqa: E402
from etcompare.config import RESULTS  # noqa: E402

BANDS = [(-1e9, 250, "<250"), (250, 500, "250–500"), (500, 1000, "500–1000"), (1000, 1500, "1000–1500"),
         (1500, 2000, "1500–2000"), (2000, 2500, "2000–2500"), (2500, 3000, "2500–3000"), (3000, 1e9, ">3000")]
LULC = {5: "Crops", 11: "Rangeland", 8: "Bare ground", 7: "Built area", 1: "Water", 2: "Trees", 4: "Flooded vegetation"}
MAP_PERIODS = {"GLDAS": ("2000-01", "2019-12"), "TerraClimate": ("2000-01", "2019-12"), "MODIS 1 km": ("2000-01", "2014-12"),
               "MODIS 500 m": ("2001-01", "2019-12"), "PML-V2": ("2000-03", "2019-12"), "SMAP L4": ("2015-04", "2019-12")}

env, tr = grid.load("001")
shape = env["basin"].shape
area = grid.area_weights(tr, shape)
cells = cwatm.load_cells_historical()
cw, key_lat, key_lon = grid.cwatm_to_grid(grid.cwatm_mean_map(cells), tr, shape)
basin = env["basin"] == 1
inside = basin & np.isfinite(cw)
print(f"basin area covered by CWatM cells: {area[inside].sum() / area[basin].sum() * 100:.1f}%")
datasets = {"CWatM": cw, **products.load_maps(tr, shape)}
elev, kop = env["elevation"], env["koppen"]

out = {"basin_mean": {}, "elevation": {}, "koppen": {}, "lulc": {}, "r_elevation_pixel": {}, "share": {}}
total = area[inside].sum()
out["share"]["elevation"] = {lab: float(area[inside & (elev >= a) & (elev < b)].sum() / total) for a, b, lab in BANDS}
out["share"]["koppen"] = {lab: float(area[inside & (kop == c)].sum() / total) for c, lab in grid.KOPPEN.items()}
out["share"]["lulc"] = {lab: float(np.nansum((env[f"lulc_{c}"] * area)[inside]) / total) for c, lab in LULC.items()}
for name, et in datasets.items():
    out["basin_mean"][name] = grid.weighted_mean(et[inside], area[inside])
    out["elevation"][name] = {}
    for a, b, lab in BANDS:
        m = inside & (elev >= a) & (elev < b)
        out["elevation"][name][lab] = grid.weighted_mean(et[m], area[m])
    out["koppen"][name] = {lab: grid.weighted_mean(et[inside & (kop == c)], area[inside & (kop == c)])
                           for c, lab in grid.KOPPEN.items()}
    out["lulc"][name] = {lab: grid.weighted_mean(et[inside], (env[f"lulc_{c}"] * area)[inside]) for c, lab in LULC.items()}
    m = inside & np.isfinite(et) & np.isfinite(elev)
    out["r_elevation_pixel"][name] = float(np.corrcoef(et[m], elev[m])[0, 1])

# environment of the 61 CWatM cells of the basin mask
weights = cwatm.load_weights("w")
cell_means = grid.cwatm_mean_map(cells)
rows = []
for _, c in weights[weights["in_mask"]].iterrows():
    m = inside & (key_lat == round(c["lat"], 2)) & (key_lon == round(c["lon"], 2))
    if not m.any():
        continue
    codes, counts = np.unique(kop[m][kop[m] > 0], return_counts=True)
    rows.append({"lat": c["lat"], "lon": c["lon"], "ET": cell_means.get((c["lat"], c["lon"]), np.nan),
                 "elevation": grid.weighted_mean(elev[m], area[m]),
                 "koppen": grid.KOPPEN.get(int(codes[np.argmax(counts)]), "") if len(codes) else "",
                 "f_crops": grid.weighted_mean(env["lulc_5"][m], area[m]),
                 "f_rangeland": grid.weighted_mean(env["lulc_11"][m], area[m]),
                 "f_bare": grid.weighted_mean(env["lulc_8"][m], area[m]),
                 "f_water": grid.weighted_mean(env["lulc_1"][m], area[m])})
cdf = pd.DataFrame(rows)
cdf.to_csv(RESULTS / "cwatm_cells_environment.csv", index=False, float_format="%.5f")
out["cells_r"] = {v: float(np.corrcoef(cdf["ET"], cdf[v])[0, 1]) for v in ["elevation", "f_crops", "f_rangeland", "f_bare", "f_water"]}
out["cells_n"] = len(cdf)
json.dump(out, open(RESULTS / "gradients.json", "w"), indent=1)

# product maps against their basin series over the same period
series = products.load_series()
check = []
for name, a in products.load_maps(tr, shape).items():
    v = a[basin]
    ok = np.isfinite(v)
    s = series[name].dropna()
    p0, p1 = (pd.Period(x, "M") for x in MAP_PERIODS[name])
    sm = float(s[(s.index >= p0) & (s.index <= p1)].mean())
    mm = float(np.sum(v[ok] * area[basin][ok]) / np.sum(area[basin][ok]))
    check.append({"dataset": name, "valid_share": ok.mean(), "map_mean": mm, "series_mean": sm, "ratio": mm / sm})
pd.DataFrame(check).to_csv(RESULTS / "product_map_check.csv", index=False, float_format="%.4f")

pd.set_option("display.width", 200)
for k in ["elevation", "koppen", "lulc"]:
    print(f"\n{k} (mm/d)")
    print(pd.DataFrame(out[k]).round(2).assign(share=pd.Series(out["share"][k]).round(3)))
print("\nbasin means", {k: round(v, 2) for k, v in out["basin_mean"].items()})
print("cell correlations", {k: round(v, 2) for k, v in out["cells_r"].items()})
print(pd.DataFrame(check).round(3))
