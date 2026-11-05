"""Build the environmental layers on the 0.01-degree analysis grid and a 0.0025-degree grid for maps.

elevation : mean SRTM 1-arcsecond elevation per pixel (tiles downloaded to data/grid/srtm)
lulc_<c>  : fraction of each Esri/Impact Observatory land-cover class in 2017
koppen    : Koppen-Geiger code at the pixel centre
basin     : pixel centre inside the basin polygon
Needs internet access. Outputs: data/grid/env_grid_001.npz, data/grid/env_grid_00025.npz
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from etcompare import grid  # noqa: E402
from etcompare.config import DATA  # noqa: E402

tiles = DATA / "grid" / "srtm"
t1, s1 = grid.make_grid(0.01)
t2, s2 = grid.make_grid(0.0025)
b1, k1 = grid.burn_basin_and_koppen(t1, s1)
b2, k2 = grid.burn_basin_and_koppen(t2, s2)
e1 = grid.elevation(t1, s1, tiles)
e2 = grid.elevation(t2, s2, tiles)
frac, dominant = grid.land_cover(t1, s1, t2, s2)

np.savez_compressed(DATA / "grid" / "env_grid_001.npz", transform=np.array(t1)[:6], basin=b1, koppen=k1, elevation=e1,
                    **{f"lulc_{c}": frac[c].astype(np.float32) for c in grid.LULC})
np.savez_compressed(DATA / "grid" / "env_grid_00025.npz", transform=np.array(t2)[:6], basin=b2, koppen=k2, elevation=e2,
                    lulc=dominant)
inside = b1 == 1
print(f"basin pixels {inside.sum()}; elevation {np.nanmin(e1[inside]):.0f}-{np.nanmax(e1[inside]):.0f} m")
for c, name in grid.LULC.items():
    print(f"  {name:20s} {np.nanmean(frac[c][inside]) * 100:5.1f}%")
