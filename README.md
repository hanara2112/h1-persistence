# Ripser + Persim evaluation

A small study of how Gaussian noise affects persistent homology of a circle,
using a filled disk as a control. It also checks sampling density and runtime.

- [Report](REPORT.md): question, methods, results, and limitations.
- [Analysis](src/analysis.py): experiments and plotting; settings are at the top.
- [Results](results/): figures, CSV tables, and the saved environment.
- [Hole-detection investigation](DETECTION_REPORT.md): independent null calibration
  and held-out detection rates, false positives, and ROC-AUC.

## Setup and run

Use Python 3.11. From the repository root:

```sh
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python src/analysis.py
```

On Windows, activate with `.venv\Scripts\activate` instead.
After setup, rerun with `python src/analysis.py`.
Run the separate detection study with `python src/hole_detection.py`; its outputs
go to `results/hole_detection/`.

The script regenerates five figures, six CSVs, and `results/environment.json`.
No data download or Jupyter installation is needed. Seeds are fixed; runtime
varies with the machine and system load. Dependencies are pinned in
[requirements.txt](requirements.txt), and the tested environment is recorded in
[environment.json](results/environment.json).
