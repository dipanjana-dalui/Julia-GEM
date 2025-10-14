"""
Running this file will install all the Julia packages 
necessary for running the scripts.

Run this only once in the beginning. 
"""

using Pkg

Pkg.add("DifferentialEquations")
Pkg.add("Plots")
Pkg.add("Distributions")
Pkg.add("DataFrames")
Pkg.add("Random")
Pkg.add("Statistics")
Pkg.add("StaticArrays")
Pkg.add("StatsBase")
Pkg.add("AlgebraOfGraphics")
Pkg.add("CSV")
Pkg.add("Makie")
Pkg.add("CairoMakie")
Pkg.add("Tidier")
Pkg.add("UnPack")
