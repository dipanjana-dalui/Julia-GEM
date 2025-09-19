""" 
This file provides the space to run your GEM model 
"""
# run only once the first time
# include("install_pkgs.jl")

# load all required packages 
include("functions/Packages.jl")

# load model definition  
include("bdLM_model_definition.jl") 
#include("2spp_model_definition.jl")

# load model configuration
include("bdLM_model_config.jl")
#include("2spp_model_config.jl")

# load all functions 
include("functions/GEM_Functions.jl")

# run the GEM simulation

using Dates
current_time = now()
run_sim = GEM_sim(
                  N0, # initial state
                  model_par_vect, # model parameters in a vector
                  design_choices, # evolution decisions
                  mappings, # parameter-to-state mapping
                  sim_params,# simulation parameters
                  sim_output, # output containers
                  verbose=false # show time on console
                  ) #


current_time = now()
# output: 
# Tuple{Array{Float64, 4}, Array{Float64, 5}, Array{Float64, 5}}

# dataframe for population time series
pop_dat = run_sim.pop_df

# 2 dataframes: mean and variance
trait_dat = run_sim.trait_df

# accessing the two trait dataframes
# trait mean dataframe:
trait_dat.median 

# trait variance dataframe
trait_dat.var

# Pop_Plot(pop data, stateID)
Pop_Plot(run_sim.pop_df, 1)

# Trait_Plot(mean, var, stateID, "trait name")
Trait_Plot(trait_dat.median, trait_dat.var, 1, "d_min")

# Geno_Plot(mean, stateID, "trait name")
Geno_Freq_Plot(trait_dat.median, 1, "g_3")