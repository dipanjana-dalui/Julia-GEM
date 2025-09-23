
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
p := parameters: (b_max, d_min, scr, f, m)
t := time

ODE: 
dN1dt = r*N1 - s*N1*N2
dn2dt = f*s*N1*N2 - m*N2

# calculate constants:
r_max = b_max-d_min


"""
# ======================================================================
#                       BIRTH-DEATH FUNCTIONS
# ======================================================================

function birth_prey(b::Float64, N::Vector{Int})
    birth = max((b*N[1]),0)
end

function death_prey(d::Float64, N::Vector{Int}, s::Float64)
    death =  d .* N[1]
    pred_death = s .* N[1] .* N[2]
    pred_death = pred_death[1]
    death = death + pred_death
end

function birth_pred(s::Float64, N::Vector{Int}, f::Float64)
    birth = f .* s .* N[1] .* N[2]
end

#function death_pred(p::ModelParameters, s::InitState)
#    return p.m * s.N[2]
#end
function death_pred(m::Float64, N::Vector{Int})
    death = m*N[2]
end

function Event_Terms(param_next::Matrix{Float64}, N::Vector{Int})
    b = param_next[1,1] # max birth
    d = param_next[1,2] # min death
    #b_s = param_next[1,3] # density dependence of birth
    #d_s = param_next[1,4]
    s = param_next[2,3]
    f = param_next[2,4]
    m = param_next[2,5]

    birth_H =  birth_prey(b, N)
    death_H =  death_prey(d, N, s)
    birth_P =  birth_pred(s, N, f)
    death_P =  death_pred(m, N) 
    return birth_H, death_H, birth_P, death_P
end



# ======================================================================
#                         SETUP STRUCTURES
# ======================================================================

"""
    1. InitState
Initial state struct
"""
mutable struct InitState
    N::Vector{Int}
end

mutable struct InitStates
    N::Vector{Number}
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
struct SimulationMaps
    state_par_match::Matrix{Int}
    state_geno_match::Matrix{Int}
    par_names::Vector{String}
    geno_names::Vector{String}
end

"""
    5. SimulationParameters
Our constants for the simulation
"""
struct SimulationParameter
    no_state::Int
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

