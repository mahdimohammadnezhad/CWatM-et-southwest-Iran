"""Basin-mean monthly ET of the retrospective and GCM-driven CWatM runs.

ET = totalET + EvapWaterBodyM + EvapoChannel (monthly means of daily values, mm/d).
Needs the CWatM output folder (environment variable CWATM_OUTPUT).
Outputs: data/cwatm/cells_historical_monthly.csv, data/cwatm/basin_monthly_et.csv
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from etcompare import cwatm  # noqa: E402
from etcompare.config import DATA, runs  # noqa: E402

weights = cwatm.load_weights("w")
cells = pd.read_csv(DATA / "cwatm" / "basin_cells.csv")

hist_dir = dict(runs())["historical"]
c = cwatm.cell_series(hist_dir, cells)
c.to_csv(DATA / "cwatm" / "cells_historical_monthly.csv", index=False, float_format="%.6g")
print("retrospective cells:", c[["lat", "lon"]].drop_duplicates().shape[0])

out = []
for name, folder in runs():
    if not folder.exists():
        print("missing", folder)
        continue
    s = cwatm.basin_series(folder, weights)
    out.append(pd.DataFrame({"run": name, "ym": s.index.astype(str), "ET": s.to_numpy()}))
    print(name, len(s), "months")
pd.concat(out).to_csv(DATA / "cwatm" / "basin_monthly_et.csv", index=False, float_format="%.6f")
