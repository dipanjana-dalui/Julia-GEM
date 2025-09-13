""" 
This file provides the space to run your GEM 
"""

# run only once the first time
# include("install_pkgs.jl")

# load packages 
include("functions/Packages.jl")

# load model definition  
include("bdLM_model_definition.jl") 

# load model configuration
include("bdLM_model_config.jl")

# load all functions 
include("functions/GEM_Functions.jl")

# run the GEM simulation
run_sim = GEM_sim(N0, model_par_vect, design_choices, mappings, 
                sim_params,sim_output)
# output: 
# Tuple{Array{Float64, 4}, Array{Float64, 5}, Array{Float64, 5}}

# dataframe for population time series
pop_dat = run_sim.pop_df

# 2 dataframes: mean and variance
trait_dat = run_sim.trait_df

# trait mean dataframe:
trait_dat.median 

# trait variance dataframe
trait_dat.var

# Pop_Plot(pop data, stateID)
Pop_Plot(run_sim.pop_df, 1)

# Trait_Plot(mean, var, stateID, "trait name")
Trait_Plot(trait_dat.median, trait_dat.var, 1, "b_s")
