#!/bin/bash
julia --project=../CrystalShiftAPI --threads=6 ../CrystalShiftAPI/src/CrystalShiftAPI.jl &
pid=$!

sleep 1

python phiddle.py $*
kill $pid
