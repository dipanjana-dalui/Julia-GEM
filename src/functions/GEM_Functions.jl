# ==================================
# the core GEM simulation function. 
# ==================================
include("GEM_main.jl")

# ==================================
#       auxiliary functions 
# ==================================
# function to calculate mean, median, varaince 
include("CalcAvgFreq.jl")

# helper functions 
include("ConveniencFunc.jl")

# function that draws offsping traits 
include("DrawNewTrait.jl")

# function to set up inital population-trait matrix 
include("InitiatePop.jl")

# function to pick an event (the "wheel-of-fortune")
include("PickEvent.jl")

# function to pick individual event happens to 
include("PickIndiv.jl")

# plot function
include("PlotFunc.jl")

# function to find the next individual
include("WhoIsNext.jl")






  