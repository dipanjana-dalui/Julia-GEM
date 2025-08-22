"""
The following packages are already added into the env if you activate
an environment inside /JGEM. Check Projects.toml for more details, and 
Manifest.toml for current versions of each dependency.

This startup file is aimed at being loaded when you first activate the env
so all the packages will be ready to use.

Note: instead of loading each package in it's entirety, you can also specify 
exactly what you are using from each package.
Example syntax: 
using StaticArrays: SVector
"""

using DifferentialEquations
using Plots
using Distributions
using DataFrames
using Random
using Statistics
using StaticArrays
using StatsBase
using AlgebraOfGraphics
using CSV
using Makie
using CairoMakie 
using Tidier
using Threads
#using Dates