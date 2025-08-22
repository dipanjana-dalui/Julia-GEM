using Pkg

Pkg.activate(".") #bring active env to current dir
Pkg.activate("./JGEM") #move it to subdir
Pkg.status()
#Status `~/Documents/Work/UNL Postdoc Research /GEM_pkg/JGEM/Project.toml` (empty project)

#add what needs to be added 
Pkg.add([
"Plots",
"Distributions",
"DataFrames",
"Random",
"Statistics",
"StaticArrays",
"BenchmarkTools",
"StatsBase",
"AlgebraOfGraphics",
"Makie",
"CairoMakie"])

Pkg.status() #will list packages. Also in project.toml

#Check path
Base.LOAD_PATH

#not added in this pkg env but in the main env
using BenchmarkTools

"""
Aug21:

"""