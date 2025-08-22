# This is the complete, refactored Julia script.
# It uses a module to contain all the functions and types,
# making the code organized and reusable.

module Julia_GEM

# load packages needed 
include("packages.jl")

# some dummy exports
export ModelParameters, PopulationState, event_terms, GEM_sim

#include files

end # end of module

