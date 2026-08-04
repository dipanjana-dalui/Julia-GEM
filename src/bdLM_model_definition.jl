"""
This is the model definition file. There are 3 parts to this document.
1. Your model definition.
2. birth-death terms, and events function.  
3. struct definitons for the model, parameters, and simulation.
"""

# ======================================================================
#           1.   MODEL 
# ======================================================================
"""
This part does not run when you call the code. This space is provided for 
the user to write down their model and refer to it when needed. 

Define your model here:
dNdt := the derivative value
N := integer state value
p := parameters: (b_max, d_min, b_s, d_s)
t := time

ODE for birth-death logistic: 
dNdt = r_max * N * (1 - N / K) 

constants:
r_max = b_max-d_min
K = floor(vec((b_max - d_min)/(b_s + d_s))[1])

list of traits (these will be unique to your system. One may 
also have parameters that are not traits under selection.)
b_max: maximum birth rate
d_min: minimum death rate
b_s: density dependence of birth
d_s: density dependence of death
"""

# ======================================================================
#          2. BIRTH-DEATH FUNCTIONS
#
# Note: Symmetric ODE systems can be generalized using linear algebraic 
# matrix notation, but care must be taken regarding dimensionality of the 
# matrix operations. 
# If your birth or death function is outputting an array corresponding to multiple states, 
# then be extra careful while using the output downstream
# ======================================================================
"""
Define the birth and death terms for each state below. 
If you have more than one interacting state, make sure to
count interaction terms in the corresponding state's birth and
death values. 

"""

function birth(b_max::Float64, b_s::Float64, R::Vector{Int})
    b_new = max((b_max - b_s*R[1])*R[1],0)
end

function death(d_min::Float64, d_s::Float64, R::Vector{Int})
    d_new = (d_min + d_s*R[1])*R[1]
end

"""
This is a VERY CRITICAL function. This function creates the 
wheel of fortune (the events), and the sampling at each run 
picks an event. 
"""

function event_terms(param_next::Matrix{Float64}, const_vect::Any,R::Vector{Int})
    b_max = param_next[1,1] # max birth
    d_min = param_next[1,2] # min death
    b_s = param_next[1,3] # density dependence of birth
    d_s = param_next[1,4] # density dependence of death 

    birth_N =  birth(b_max, b_s, R)
    death_N =  death(d_min, d_s, R)
    return birth_N, death_N
end

# ======================================================================
# ======================================================================
#          3.  STRUCTURES
# only make changes to the following block if you are addint new struct
# ======================================================================
# ======================================================================

"""
    1. InitState
Initial state struct
"""
mutable struct InitState
    N::Vector{Int} # initial state abundances vector
end

"""
    2. ModelParVector
"""
struct ModelParVector{T}
    param_init::Vector{T} # initial parameter vector
end

"""
    3. DesignChoices{T}
Eco-evo choices
"""

struct DesignChoices
    h2::Array{Float64} # narrow sense heritability
    cv::Array{Float64} # coefficient of variation
    GEM_ver::Vector{String} # GEM version
end

"""
    4. SimulationMapping
"""
struct SimulationMaps
    state_par_match::Matrix{Int} # match parameters to state 
    state_geno_match::Matrix{Int} # match genotype to state
    par_names::Vector{String} # parameter names vector 
    geno_names::Vector{String} # genotype names vector 
end

"""
    5. SimulationParameters
"""
struct SimulationParameter
    no_state::Int # total number of states 
    no_param::Int # total number of parameter 
    no_columns::Int # total number of columns 
    num_time_steps::Int # standardized number of time steps 
    num_rep::Int # number of replication per GEM version 
    t_max::Float64 # total time to run simulation for 
    min_time_step_to_store::Float64 # time points to store data at 
end

"""
    6. GEMOutput: storage for all replicates and versions of GEM
"""
struct GEMOutput
    pop_stand_out_all::Array{Float64, 4} # population data storage 
    x_stand_out_all::Array{Float64, 5} # trait median & genotype frequency data storage 
    x_var_stand_out_all::Array{Float64, 5} # trait variance data storage
end

"""
7. CONSTANTS
"""
struct GEMConstant{T} 
    const_vect::Vector{T}
end

