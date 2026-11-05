"""Run the analysis from the files in data/ and write results/ and figures/.

    python run_all.py              # steps 03, 05, 06 and 08 (data in the repository only)
    python run_all.py --full       # also 01, 02, 04 and 07 (CWatM output, climate forcing, internet)
"""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
FULL = ["01_basin_cells.py", "02_cwatm_series.py", "03_evaluation.py", "04_environment_grid.py",
        "05_spatial_gradients.py", "06_projections.py", "07_climate_associations.py", "08_figures.py"]
LOCAL = ["03_evaluation.py", "05_spatial_gradients.py", "06_projections.py", "08_figures.py"]

for script in FULL if "--full" in sys.argv else LOCAL:
    print(f"\n=== {script}")
    subprocess.run([sys.executable, str(HERE / "scripts" / script)], check=True)
