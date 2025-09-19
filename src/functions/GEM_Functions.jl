"""
A list of all functions
"""
# the core GEM simulation function. 
include("GEM_main.jl")

# function to calculate mean, median, varaince 
include("CalcAvgFreq.jl")

# function to set up inital population-trait matrix 
include("InitiatePop.jl")

# function to pick individual event happens to 
include("PickIndiv.jl")

# function to find the next individual
include("WhoIsNext.jl")

# function to pick an event (the "wheel-of-fortune")
include("PickEvent.jl")

# function that draws offsping traits 
include("DrawNewTrait.jl")

# plot function
include("PlotFunc.jl")

# other functions 
include("ConveniencFunc.jl")



  