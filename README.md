# Phiddle: GUI for phase identification for laser experiments
## Installation
1. Install [uv](https://docs.astral.sh/uv/getting-started/installation/)
```console
curl -LsSf https://astral.sh/uv/install.sh | sh
```

2. Install Phiddle
```console
git clone git@github.com:MingChiangChang/phiddle.git
cd phiddle
uv sync
```
This creates a `.venv` with the exact dependency versions in `uv.lock` (uv will also fetch Python 3.11 if needed).

Without uv, `pip install -r requirements.txt` (the same locked versions) or `pip install .` (latest compatible versions) also work.

3. Install CrystalShiftAPI

Install [julia](https://julialang.org/downloads/), then install CrystalShift and CrystalTree using the Julia package manager (from the repository root):
```console
julia --project=CrystalShiftAPI -e 'using Pkg; Pkg.instantiate()'
```

## Usage
Start Phiddle from the repository root:
```console
./phiddle/start.sh
```
The script starts the phase labeling backend, waits for it to be ready, then launches the GUI (through `uv run` if uv is installed, otherwise with the active `python`). It can be run from any directory, and double clicking it also works. Optional arguments preload files:
```console
./phiddle/start.sh --h5 data.h5 --csv phases.csv --json progress.json
```

To run the two parts separately:
```console
julia --project=CrystalShiftAPI --threads=6 CrystalShiftAPI/src/CrystalShiftAPI.jl
cd phiddle && uv run python phiddle.py
```
