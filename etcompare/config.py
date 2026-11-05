"""Paths and study settings.

Raw model output and climate forcing are large and are not part of the repository. Point the
environment variables below to local copies before running the extraction steps
(scripts 02 and 07); all other steps use the files in ``data/``.
"""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"

# CWatM output: <CWATM_OUTPUT>/historical and <CWATM_OUTPUT>/<GCM>/<SSP>, monthly NetCDF files
CWATM_OUTPUT = Path(os.environ.get("CWATM_OUTPUT", "E:/CWatM_output/observed"))
# climate forcing used by CWatM (GSWP3-W5E5 and ISIMIP3b), daily NetCDF files with pr and tas
FORCING_HIST = Path(os.environ.get("FORCING_HIST", "E:/GSWP3-W5E5"))
FORCING_FUTURE = Path(os.environ.get("FORCING_FUTURE", "E:/CWatM_input/Future"))

BASIN_SHP = DATA / "basin" / "3basins.shp"
SUB_BASINS = {"Karun": "basin1.shp", "Karkheh": "basin2.shp", "Marzi-e-Gharb": "basin3.shp"}
KOPPEN_SHP = DATA / "basin" / "koppen_geiger_3basins.shp"

GCMS = ["GFDL", "IPSL", "MPI", "MRI", "UKESM"]
SSPS = {"126": "SSP1-2.6", "370": "SSP3-7.0", "585": "SSP5-8.5"}
PERIODS = {"2025-2049": (2025, 2049), "2050-2074": (2050, 2074), "2075-2099": (2075, 2099)}
BASELINE = (1995, 2019)
OVERLAP = (2015, 2019)          # retrospective and GCM-driven runs overlap
EVALUATION = ("2000-01", "2019-12")
SEASONS = {"Winter": (12, 1, 2), "Spring": (3, 4, 5), "Summer": (6, 7, 8), "Autumn": (9, 10, 11)}
MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

# ET components written by CWatM (monthly means of daily values, m/day)
ET_COMPONENTS = ["totalET", "EvapWaterBodyM", "EvapoChannel"]


def runs():
    """(name, directory) of the retrospective run and the 15 GCM-driven runs."""
    out = [("historical", CWATM_OUTPUT / "historical")]
    out += [(f"{g}_{s}", CWATM_OUTPUT / g / s) for g in GCMS for s in SSPS]
    return out
