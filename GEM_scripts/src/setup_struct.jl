""" 
This is the struct file for setup
"""
# ======================================================================
#           SETUP STRUCTURES
# ======================================================================
"""
    InitState
Initial state struct
"""
mutable struct InitState
    N::Vector{Int}
end

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
struct DesignChoices
    h2_vect::Matrix{Float64}
    cv_vect::Matrix{Float64}
    GEM_ver::Vector{String}
end


struct SimulationMapping
    state_par_match::Matrix{Int64}
    state_geno_match::Matrix{Int64}
    geno_par_match::Matrix{Int64}
    par_names::Vector{String}
    geno_names::Vector{String}
end

struct GEMOutput
    pop_stand::Array{Float64, 4}
    trait_moment1::Array{Float64, 5}
    trait_moment2::Array{Float64, 5}
end


# ======================================================================
# BIRTH-DEATH FUNCTIONS
# ======================================================================
"""
Issues with defining the functions with struct ModelParamsters:
 - the functions below will go to struct ModelParameters 
 and look for b_max. But when being used to calculate the
event_terms, we don't have a struct. It will pull the 
b_max value from the param_next structured created evolve_system_in_time_with_extinction_check
the GEM loop. For now, stick to the former definitiond 
without going through the struct for this set of funcs. 
"""

#function birth_prey(p::ModelParameters, s::InitState)
#    return max(p.b_max * s.N[1], 0.0)
#end
function birth_prey(b_max::Float64, N::Vector{Int64})
    birth = max((b_max*N[1]),0)
end

#function death_prey(p::ModelParameters, s::InitState)
#    death = p.d_min * s.N[1]
#    pred_death = p.scr * s.N[1] * s.N[2]
#    return death + pred_death
#end
function death_prey(d_min::Float64, N::Vector{Int64}, scr::Float64)
    death =  d_min .* N[1]
    pred_death = scr .* N[1] .* N[2]
    pred_death = pred_death[1]
    death = death + pred_death
end


#function birth_pred(p::ModelParameters, s::InitState)
#    return p.fec * p.scr * s.N[1] * s.N[2]
#end
function birth_pred(scr::Float64, N::Vector{Int64}, fec::Float64)
    birth = fec .* scr .* N[1] .* N[2]
end

#function death_pred(p::ModelParameters, s::InitState)
#    return p.m * s.N[2]
#end
function death_pred(m::Float64, N::Vector{Int64})
    death = m*N[2]
end

function event_terms(params_next::Matrix{Float64}, s::InitState)
    b_max = params_next[1,1] # max birth
    d_min = params_next[1,2] # min death
    #b_s = params_next[1,3] # density dependence of birth
    #d_s = params_next[1,4]
    scr = params_next[2,3]
    fec = params_next[2,4]
    m = params_next[2,5]

    birth_H =  birth_prey(b_max, N)
    death_H =  death_prey(d_min, N, scr)
    birth_P =  birth_pred(scr, N, fec)
    death_P =  death_pred(m, N) 
    return birth_H, death_H, birth_P, death_P
end

