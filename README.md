# Phiddle: GUI for phase identification for laser experiments
## Installation
1. Install Phiddle with [uv](https://docs.astral.sh/uv/)
```console
git clone git@github.com:MingChiangChang/phiddle.git
cd phiddle
uv sync
```
Without uv, `pip install .` (latest compatible versions) or `pip install -r requirements.txt` (the exact locked versions) also work.
2. Install CrystalShiftAPI
Install [julia](https://julialang.org/downloads/).
Then install CrystalShift and CrystalTree using julia package manager
```julia
julia --project=../CrystalShiftAPI -e 'using Pkg; Pkg.instantiate()'
```

## Usage
1. Start phiddle
``` console
./start.sh
```
The shell script starts both the phase labeling backend and GUI frontend.
Double clicking the `start.sh` shell script should also work for Window systems.
