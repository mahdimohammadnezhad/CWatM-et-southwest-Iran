# cwatm-et-southwest-iran

Code and derived data for the evaluation of actual evapotranspiration (ET) simulated by the
Community Water Model (CWatM) in the Karun, Karkheh, and Marzi-e-Gharb basins of southwest Iran,
and for its projected change under CMIP6 climate scenarios.

Related article: Mohammadnezhad, M., Rezaee, A., & Smilovic, M. Evaluating CWatM evapotranspiration
against six global products and its CMIP6 projections in southwest Iran's regulated basins for
agricultural water planning (submitted to *Agricultural Water Management*).

## What the code does

1. Basin-mean monthly CWatM ET (land-surface ET + lake/reservoir + channel evaporation) for the
   retrospective run (GSWP3-W5E5, 1981–2019) and 15 GCM-driven runs (5 GCMs × SSP1-2.6, SSP3-7.0,
   SSP5-8.5, 2015–2100).
2. Agreement with six reference datasets (GLDAS-Noah, TerraClimate, MOD16A2 on 1 km and 500 m grids,
   PML-V2, SMAP L4): r, RMSE, PBIAS, NSE, KGE, and PBIAS by flux class and calendar quarter.
3. Mean ET of all datasets by elevation band, Köppen–Geiger class, and land-cover class on a common
   0.01° grid.
4. Projected changes for 2025–2049, 2050–2074, and 2075–2099 relative to 1995–2019, after scaling each
   GCM-driven series to the retrospective run over 2015–2019.
5. Inter-annual association of annual ET with precipitation and air temperature.

## Repository layout

```
etcompare/        Python package (I/O, statistics, grids, projections, plotting)
scripts/          numbered analysis steps
gee/              Google Earth Engine script for the reference products
data/             inputs and derived data (see data/README.md)
results/          tables written by the scripts
figures/          figures written by the scripts
run_all.py
```

## Installation

Python 3.10 or later.

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

NetCDF files are read through the GDAL netCDF driver included in the `rasterio` wheels.

## Running

From the repository root:

```bash
python run_all.py            # steps 03, 05, 06, 08: uses only the files in data/
python run_all.py --full     # all steps; needs CWatM output, climate forcing and internet access
```

| Step | Script | Needs | Output |
|---|---|---|---|
| 01 | `01_basin_cells.py` | basin polygon | `data/cwatm/basin_cells.csv` |
| 02 | `02_cwatm_series.py` | CWatM output (`CWATM_OUTPUT`) | `data/cwatm/*.csv` |
| 03 | `03_evaluation.py` | data/ | Table 1, flux-class and quarterly PBIAS, mean adjustment |
| 04 | `04_environment_grid.py` | internet (SRTM tiles, Planetary Computer) | `data/grid/*.npz` |
| 05 | `05_spatial_gradients.py` | data/ | ET by elevation, climate, and land cover |
| 06 | `06_projections.py` | data/ | Table 2, Tables S1–S2, monthly signals, annual series |
| 07 | `07_climate_associations.py` | CWatM daily output, forcing (`FORCING_HIST`, `FORCING_FUTURE`) | Table S3 |
| 08 | `08_figures.py` | results/ | Figs. 2–8, S1–S2 |

Paths to the raw data are set with environment variables (defaults in `etcompare/config.py`):

```bash
export CWATM_OUTPUT=/path/to/CWatM_output/observed     # historical/ and <GCM>/<SSP>/ subfolders
export FORCING_HIST=/path/to/GSWP3-W5E5
export FORCING_FUTURE=/path/to/ISIMIP3b                # <GCM>/<SSP>/ subfolders
```

## Outputs and article items

| Article | File |
|---|---|
| Fig. 2 | `figures/fig2_environment.png` |
| Fig. 3, Table 1 | `figures/fig3_seasonal_cycle.png`, `results/table1_agreement.csv` |
| Fig. 4 | `figures/fig4_pbias.png`, `results/pbias_by_class.csv` |
| Fig. 5 | `figures/fig5_mean_et_maps.png` |
| Fig. 6 | `figures/fig6_gradients.png`, `results/gradients.json` |
| Fig. 7 | `figures/fig7_annual_et.png`, `results/annual_series.csv`, `results/mean_adjustment.csv` |
| Fig. 8 | `figures/fig8_monthly_change.png`, `results/monthly_signals.csv` |
| Table 2 | `results/table2_changes.csv` |
| Table S1 | `results/tableS1_overlap.csv` |
| Table S2 | `results/tableS2_statistics.csv` |
| Table S3, Fig. S2 | `results/tableS3_relationships.csv`, `figures/figS2_et_climate.png` |
| Fig. S1 | `figures/figS1_monthly_distributions.png` |

Fig. 1 (study area and infrastructure) is taken from Mohammadnezhad et al. (2026).

## Notes

- The basin-mean CWatM series use the 61 cells of the study-basin mask with cell-area weights.
  Weighting cells by their area inside the basin polygon (`w_poly`) changes the 2000–2019 mean by 0.4%.
- The additive mean adjustment of GLDAS-Noah and SMAP L4 is used only to compare inter-annual
  variability in Fig. 7; all statistics use the original series.
- Table 2 values reproduce the article to within 0.1 percentage point; the article used spherical
  cell areas and this code uses ellipsoidal (EPSG:6933) cell areas.
- Maps show the contiguous basin complex; all statistics use the full basin polygon.

## Data sources

CWatM: Burek et al. (2019), https://github.com/iiasa/CWatM. Climate forcing: ISIMIP3a GSWP3-W5E5 and
ISIMIP3b (https://data.isimip.org). Reference ET: Google Earth Engine data catalogue. Elevation: SRTM
1-arcsecond (AWS Terrain Tiles). Land cover: Esri/Impact Observatory 10 m annual land use/land cover
2017 (Microsoft Planetary Computer). Climate classes: Kottek et al. (2006).

## License

Code: MIT (see `LICENSE`). Derived data in `data/` are provided under CC BY 4.0; the original
datasets keep their own licenses.

## Citation

See `CITATION.cff`.
