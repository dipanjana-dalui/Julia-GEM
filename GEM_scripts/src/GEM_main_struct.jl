



# Include the functions for your simulation (run_replicate, GEM_sim, etc.)
include("simulation_functions.jl")


# Define a struct for simulation results to make the code cleaner
struct SimulationResults
    pop_df::DataFrame
    par_mean_df::DataFrame
    par_var_df::DataFrame
end

# Define a struct for the simulation parameters
struct SimulationParameters
    GEM_ver::Vector{String}
    t_max::Float64
    no_species::Int
    num_rep::Int
    no_columns::Int
    no_param::Int
    N0::Vector{Int}
    which_par_quant::Matrix{Int}
    state_geno_match::Matrix{Int}
    state_par_match::Matrix{Int}
    param_init::Vector{Float64}
    cv_vect::Matrix{Float64}
    h2_vect::Matrix{Float64}
    par_names::Vector{String}
    geno_names::Vector{String}
    stand_time::Vector{Float64}
end

