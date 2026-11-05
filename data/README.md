# Data

| Path | Content |
|---|---|
| `basin/3basins.*` | study-basin polygon (Karun, Karkheh, Marzi-e-Gharb), WGS 84 |
| `basin/basin1–3.*` | Karun, Karkheh, and Marzi-e-Gharb sub-basins |
| `basin/koppen_geiger_3basins.*` | Köppen–Geiger classes clipped to the basin (Kottek et al., 2006); no `.prj`, coordinates in EPSG:3857; `GRIDCODE` 22 BWh, 26 BSk, 27 BSh, 34 Csa, 45 Dsa, 46 Dsb |
| `cwatm/basin_mask_cells.csv` | centres of the 61 CWatM cells (0.5°) of the study-basin mask |
| `cwatm/basin_cells.csv` | cell area (`cell_km2`), share inside the polygon (`frac`), mask weights `w`, polygon weights `w_poly` |
| `cwatm/cells_historical_monthly.csv` | retrospective run, monthly mean ET components per cell (mm/d): `totalET`, `waterbody` (EvapWaterBodyM), `channel` (EvapoChannel) |
| `cwatm/basin_monthly_et.csv` | basin-mean monthly ET (mm/d) of the retrospective run and the 15 GCM-driven runs (`run` = `historical` or `<GCM>_<SSP>`) |
| `climate/basin_monthly_climate.csv` | monthly CWatM land-surface ET (`totalET`) and precipitation (mm/month) and mean air temperature (°C) of the 16 runs, simple mean over the cells touching the basin polygon (script 07) |
| `products/et_products_basin_monthly.csv` | basin-mean ET (mm/d) of the six reference datasets per monthly window (GEE export; `date` = window start); TerraClimate values are monthly totals / 30 |
| `products/maps/ET_mean_<dataset>.tif` | mean ET (mm/d) of each dataset on the 0.01° grid, over the periods of Table 1 |
| `grid/env_grid_001.npz` | 0.01° grid: `basin`, `koppen`, `elevation` (m), `lulc_<class>` fractions, `transform` |
| `grid/env_grid_00025.npz` | 0.0025° grid for maps: `basin`, `koppen`, `elevation`, dominant land-cover class `lulc` |

Periods of the reference datasets: GLDAS-Noah and TerraClimate 2000–2019, MODIS 1 km 2000–2014,
MODIS 500 m 2001–2019, PML-V2 March 2000–2019, SMAP L4 April 2015–2019.

Land-cover classes (Esri/Impact Observatory): 1 water, 2 trees, 4 flooded vegetation, 5 crops,
7 built area, 8 bare ground, 9 snow/ice, 10 clouds, 11 rangeland.
