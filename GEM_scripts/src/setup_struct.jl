""" 
This is the struct file for setup
"""
# ======================================================================
# SETUP PARAMETER STRUCTURES
# ======================================================================

"""
    ModelParameters{T}
This struct holds the model initial parameters. These are not mutable 
"""
struct ModelParameters{Float64}
    b_max::Float64  # Maximum birth rate of prey
    d_min::Float64  # Minimum death rate of prey
    scr::Float64    # Spsce clearace rate
    fec::Float64    # Fecundity of predators/pathogens
    m::Float64      # Mortality rate of predators/pathogens
end

"""
    InitState
Initial state struct
"""
mutable struct InitState
    N::Vector{Int}
end

"""
    SimulationParameters
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
    DesignChoices{T}
Eco-evo choices
"""
struct DesignChoices{T}
    h2_vect::Matrix{T}
    cv_vect::Matrix{T}
end

struct SimulationOutput
    pop_stand_out_all::Array{Float64, 4}
    x_stand_out_all::Array{Float64, 5}
    x_var_stand_out_all::Array{Float64, 5}
end

struct SimulationMapping
    state_par_match::Matrix{Int64}
    state_geno_match::Matrix{Int64}
    geno_par_match::Matrix{Int64}
    par_names::Vector{String}
    geno_names::Vector{String}
end

# ======================================================================
# BIRTH-DEATH FUNCTIONS
# ======================================================================

function birth_prey(p::ModelParameters, s::InitState)
    return max(p.b_max * s.N[1], 0.0)
end

function death_prey(p::ModelParameters, s::InitState)
    death = p.d_min * s.N[1]
    pred_death = p.scr * s.N[1] * s.N[2]
    return death + pred_death
end

function birth_pred(p::ModelParameters, s::InitState)
    return p.fec * p.scr * s.N[1] * s.N[2]
end

function death_pred(p::ModelParameters, s::InitState)
    return p.m * s.N[2]
end

"""
struct ParNext
    


end
"""
function event_terms(p::ParNext, s::InitState)
    birth_H = birth_prey(p, s)
    death_H = death_prey(p, s)
    birth_P = birth_pred(p, s)
    death_P = death_pred(p, s)
    return birth_H, death_H, birth_P, death_P
end

