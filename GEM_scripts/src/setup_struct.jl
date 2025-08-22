""" 
This is the struct file for setup
"""
# ======================================================================
# SETUP PARAMETER STRUCTURES
# ======================================================================

"""
    ModelParameters{T}
This struct holds the model initial parameters. These are not mutavle 
"""
struct ModelParameters{Float64}
    b_max::T  # Maximum birth rate of prey
    d_min::T  # Minimum death rate of prey
    scr::T    # Scavenging/interaction rate
    fec::T    # Fecundity of predators/pathogens
    m::T      # Mortality rate of predators/pathogens
end

"""
    PopulationState
Initial population struct
"""
mutable struct PopulationState
    N::Vector{Int64}
end

"""
    StateCount
State count and param count struct
"""
struct Count{Int64}
    no_species
    no_param
end


"""
    SimulationParameters

Holds all the constant settings for the simulation run.
"""
struct SimulationParameters
    num_rep::Int64
    t_max::Float64
    min_time_step_to_store::Float64
end

"""
    DesignChoices{T}

Holds the matrices related to design choices like heritability and CV.
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

function birth_prey(p::ModelParameters, s::PopulationState)
    return max(p.b_max * s.N[1], 0.0)
end

function death_prey(p::ModelParameters, s::PopulationState)
    death = p.d_min * s.N[1]
    pred_death = p.scr * s.N[1] * s.N[2]
    return death + pred_death
end

function birth_pred(p::ModelParameters, s::PopulationState)
    return p.fec * p.scr * s.N[1] * s.N[2]
end

function death_pred(p::ModelParameters, s::PopulationState)
    return p.m * s.N[2]
end

function event_terms(p::ModelParameters, s::PopulationState)
    birth_H = birth_prey(p, s)
    death_H = death_prey(p, s)
    birth_P = birth_pred(p, s)
    death_P = death_pred(p, s)
    return birth_H, death_H, birth_P, death_P
end

