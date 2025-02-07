# Phiddle: GUI for phase identification for laser experiments
## Installation
1. Install CrystalShift 
Install [julia](https://julialang.org/downloads/).
Then install CrystalShift and CrystalTree using julia package manager
```julia
using Pkg
Pkg.add(url="https://github.com/MingChiangChang/CrystalShift.jl")
Pkg.add(url="https://github.com/MingChiangChang/CrystalTree.jl")
```
2. Install Phiddle
```console
git clone git@github.com:MingChiangChang/phiddle.git
cd phiddle
pip install .
```

## Usage
1. Start phiddle
``` console
./start.sh
```
The shell script starts both the phase labeling backend and GUI frontend.
