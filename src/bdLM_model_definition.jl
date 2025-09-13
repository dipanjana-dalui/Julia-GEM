"""
This is the model definition file. There are 3 parts to this document.
1. Your model definition.
2. Birth-death terms, and events function.  
3. struct definitons for the model, parameters, and simulation.
"""

# ======================================================================
#           1.   MODEL 
# ======================================================================
"""
Define your model here:
dNdt := the derivative value
N := integer state value
p := parameters: (b_max, d_min, b_s, d_s)
t := time

ODE: 
dNdt = r_max * N * (1 - N / K) 

# calculate constants:
r_max = b_max-d_min
K = floor(vec((b_max - d_min)/(b_s + d_s))[1])

"""

# ======================================================================
#          2. BIRTH-DEATH FUNCTIONS
# ======================================================================
"""
Define the birth and death terms for each state below. 
"""
function Birth(b_max, b_s, R)
    b_new = max((b_max - b_s*R[1])*R[1],0)
end

function Death(d_min, d_s, R)
    d_new = (d_min + d_s*R[1])*R[1]
end


function Event_Terms(param_next::Matrix{Float64}, R::Vector{Int})
    b_max = param_next[1,1] # max birth
    d_min = param_next[1,2] # min death
    b_s = param_next[1,3] # density dependence of birth
    d_s = param_next[1,4] # density dependence of death 

    birth =  Birth(b_max, b_s, R)
    death =  Death(d_min, d_s, R)
    return birth, death
end

# ======================================================================
#          3.  STRUCTURES
# ======================================================================
"""
    1. InitState
Initial state struct
"""
mutable struct InitState
    N::Vector{Int}
end

"""
    2. ModelParVector
"""
struct ModelParVector{T}
    param_init::Vector{T}
end

"""
    3. DesignChoices{T}
Eco-evo choices
"""
struct DesignChoice
    h2::Matrix{Float64}
    cv::Matrix{Float64}
    GEM_ver::Vector{String}
end

struct DesignChoices
    h2::Array{Float64}
    cv::Array{Float64}
    GEM_ver::Vector{String}
end

"""
    4. SimulationMapping
"""
struct SimulationMap
    state_par_match::Matrix{Int}
    state_geno_match::Matrix{Int}
    geno_par_match::Matrix{Int}
    which_par_quant::Matrix{Int}
    par_names::Vector{String}
    geno_names::Vector{String}
end

"""
    5. SimulationParameters
Our constants for the simulation
"""
struct SimulationParameters
    no_species::Int
    no_param::Int
    no_columns::Int 
    num_time_steps::Int
    num_rep::Int
    t_max::Float64
    min_time_step_to_store::Float64 
end

"""
    6. GEMOutput
"""
struct GEMOutput
    pop_stand_out_all::Array{Float64, 4}
    x_stand_out_all::Array{Float64, 5}
    x_var_stand_out_all::Array{Float64, 5}
end


