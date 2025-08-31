module JGEM

# load packages
include("functions/Packages.jl")

# some dummy exports
#export ModelParameters, PopulationState, event_terms, GEM_sim

# include files
# 1. Load all dependent functions
include("functions/AuxiliaryFunc.jl")

# 2. Setup
# 2.1 Example 2 species model
include("setup_struct.jl") # definations of all of your setup params
include("setup_config.jl") # instantiating the parameters 

# 2.2 Example birth-death Logistic Model
include("bdLM_setup_struct.jl") # definations of all of your setup params
include("bdLM_setup_config.jl")

# Include the main GEM file 
include("GEM_main.jl")

run_sim = GEM_sim(N0, model_par_vect, design_choices, mappings, 
                sim_params,sim_output)

pop_df = make_pop_df_long(sim_output, sim_params, design_choices)
trait_df = make_trait_df_long(sim_output, sim_params, design_choices, mappings)

trait_df[1]

typeof(trait_df)

Pop_Plot(pop_df, 1)
Trait_Plot()

"""
to check how many cores are available
versioninfo()

julia --threads=4
or 
julia --threads=auto
check:
using Base.Threads
Threads.nthreads()
To actually use these threads in your code, you'll need to use macros like @threads for parallelizing loops 
    or @spawn for creating tasks that run on different threads.

"""


end # module JGEM
