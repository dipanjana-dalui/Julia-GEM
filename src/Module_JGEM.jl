module JGEM

# load packages
include("functions/Packages.jl")

# functions you wish to be available to use from this module
export run_replicate, GEM_sim, make_trait_df_long, make_pop_df_long,
Pop_Plot, Trait_Plot

# include files
# 1. Load all dependent functions
include("functions/AuxiliaryFunc.jl")
include("bdLM_setup_struct.jl") 

include("functions/GEM_main.jl")


# 2. Setup
# 2.1 Example birth-death Logistic Model
#include("bdLM_setup_struct.jl") # definations of all of your setup params


# 2.1 Example 2 species model
#include("setup_struct.jl") # definations of all of your setup params
#include("setup_config.jl") # instantiating the parameters 

# 2.3. Include any other example model you want to load 

# 3. Load the core GEM function 



end

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


