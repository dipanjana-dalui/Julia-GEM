""" 
This is the entru script.
This file provides the space to run your GEM model.
Before you start, please go through README.txt
"""
# run only once the first time
# include("Install_Pkgs.jl")

# load all required packages 
include("functions/Packages.jl")

# load model definition file for corresponding model 
include("bdLM_model_definition.jl") 

include("2spp_model_definition.jl") 

# load model configuration for corresponding model
include("bdLM_model_config.jl")

include("2spp_model_config.jl")

# load all functions 
include("functions/GEM_functions.jl")


# run the GEM simulation

start_t = now() # run this line to save staring time
run_sim = GEM_sim(
                  N0, # initial state
                  model_par_vect, # model parameters in a vector
                  gem_const_vect,
                  design_choices, # evolution decisions
                  mappings, # parameter-to-state mapping
                  sim_params,# simulation parameters
                  sim_output;
                  verbose=false # show time on console
                  ) #
end_t = now() # run this line to save ending time



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

# pop_plot(pop data, stateID)
p = pop_plot(pop_data=run_sim.pop_df, stateID= 1, add_mean=true)
save("population_plot.svg", p)

# trait_plot(mean, var, stateID, "trait name")
t = trait_plot(mediandf=trait_dat.median, vardf=trait_dat.var, 
stateID=1, trait_to_plot="d_min", add_mean=true)
save("d_min_trait_mean_plot.svg", t)

# geno_plot(mean, stateID, "trait name")
g = geno_freq_plot(freqdf=trait_dat.median, stateID= 1, geno_names= "g_1")
save("geno_freq_plot.svg", g)
