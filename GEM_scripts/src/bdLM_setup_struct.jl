""" 
This is the struct file for setup
"""
# ======================================================================
#           SETUP STRUCTURES
# ======================================================================
"""
    1. InitState
Initial state struct
"""
mutable struct InitState
    N::Vector{Int}
end

"""
    2. ModelParameters{T}
This struct holds the model initial parameters. These are not mutable 
"""
struct ModelParameters{Float64}
    b_max::Float64  # Maximum birth rate of prey
    d_min::Float64  # Minimum death rate of prey
    b_s::Float64    # Density dependence of birth
    d_s::Float64    # Density dependence of death
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
function Birth(b_max, b_s, R)
    b_new = max((b_max - b_s*R[1])*R[1],0)
end

function Death(d_min, d_s, R)
    d_new = (d_min + d_s*R[1])*R[1]
end


function event_terms(param_next::Matrix{Float64}, R::Vector{Int})
    b_max = param_next[1,1] # max birth
    d_min = param_next[1,2] # min death
    b_s = param_next[1,3] # density dependence of birth
    d_s = param_next[1,4]

    birth =  Birth(b_max, b_s, R)
    death =  Death(d_min, d_s, R)
    return birth, death
end


