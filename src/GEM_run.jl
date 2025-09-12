""" 
This file provides the space to run your GEM 
"""

#include("Module_JGEM.jl")

include("functions/Packages.jl")

include("bdLM_model_definition.jl") 
include("bdLM_model_config.jl")

include("functions/AuxiliaryFunc.jl")



#include("functions/GEM_main.jl")

#check what functions are exported and are ready for use
#names(JGEM)

run_sim = GEM_sim(N0, model_par_vect, design_choices, mappings, 
                sim_params,sim_output)
#Tuple{Array{Float64, 4}, Array{Float64, 5}, Array{Float64, 5}}

# dataframe for population time series
pop_dat = run_sim.pop_df

# 2 dataframes: mean and variance
trait_dat = run_sim.trait_df

# trait mean dataframe:
trait_dat.mean

# trait variance dataframe
trait_dat.var

# Pop_Plot(pop data, stateID)
Pop_Plot(run_sim.pop_df, 1)

# Trait_Plot(mean, var, stateID, "trait name")
Trait_Plot(trait_dat.mean, trait_dat.var, 1, "b_max")
