module JGEM

# load packages
include("functions/Packages.jl")

# some dummy exports
#export ModelParameters, PopulationState, event_terms, GEM_sim

# include files
# 1. Load all dependent functions
include("functions/AuxiliaryFunc.jl")

# 2. Setup
include("setup_struct.jl") # definations of all of your setup params
include("setup_config.jl") # instantiating the parameters 

include("bdLM_setup_struct.jl") # definations of all of your setup params
include("bdLM_setup_config.jl")


run = run_replicate(N0, model_params, design_choices, mappings, 
                sim_params,1)

run.pop_time_series

run.trait_mom1
run.trait_mom2

end # module JGEM
