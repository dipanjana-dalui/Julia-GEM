
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
    scr::Float64    # Spsce clearace rate
    fec::Float64    # Fecundity of predators/pathogens
    m::Float64      # Mortality rate of predators/pathogens
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
struct GEMSimOutput
    pop_stand_out_all::Array{Float64, 4}
    x_stand_out_all::Array{Float64, 5}
    x_var_stand_out_all::Array{Float64, 5}
end

"""
    7. ModelParVector
"""
struct ModelParVector{T}
    param_init::Vector{T}
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
function birth_prey(b::Float64, N::Vector{Int64})
    birth = max((b*N[1]),0)
end

#function death_prey(p::ModelParameters, s::InitState)
#    death = p.d_min * s.N[1]
#    pred_death = p.scr * s.N[1] * s.N[2]
#    return death + pred_death
#end
function death_prey(d::Float64, N::Vector{Int64}, s::Float64)
    death =  d .* N[1]
    pred_death = s .* N[1] .* N[2]
    pred_death = pred_death[1]
    death = death + pred_death
end


#function birth_pred(p::ModelParameters, s::InitState)
#    return p.fec * p.scr * s.N[1] * s.N[2]
#end
function birth_pred(s::Float64, N::Vector{Int64}, f::Float64)
    birth = f .* s .* N[1] .* N[2]
end

#function death_pred(p::ModelParameters, s::InitState)
#    return p.m * s.N[2]
#end
function death_pred(m::Float64, N::Vector{Int64})
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


