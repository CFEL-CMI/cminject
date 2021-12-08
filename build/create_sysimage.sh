#!/bin/bash
set -e

tmpout="tmpout.h5"
rm -f "$tmpout"

julia -t 2 --project=@. --startup-file=no --trace-compile=precompile_statements.jl ../cminject.jl \
    -f ../data/nitrogen_1.8mbar_extended_nan.h5 -n 100\
    --x G[0,1e-3] --z -0.1285 --vx G[0,0.1] --vz 10 --r 13.5e-9\
    --rho 19320 -d -0.001 0.0 0.001\
    -o "$tmpout" -t 0 0.05 --plot

julia --project=@. --startup-file=no -e '
using PackageCompiler;
create_sysimage([:ArgParse, :Plots, :DifferentialEquations, :Random, :GR, :PlotlyBase, :CMInject];
    sysimage_path="CMInject.so",
    precompile_statements_file=["precompile_statements.jl"])
'

rm "$tmpout"
