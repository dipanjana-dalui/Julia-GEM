module JGEM

# load packages
include("Packages.jl")

# some dummy exports
export ModelParameters, PopulationState, event_terms, GEM_sim

# include files

# 1. Setup
include("setup_struct.jl")
include("prelim_setup.jl")


# 2. Load all dependent functions
include("AuxiliaryFunc.jl")



include("GEM_main_struct.jl")

end # module JGEM
