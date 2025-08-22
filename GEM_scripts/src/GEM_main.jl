"""

"""

# Include the files with the struct definitions and parameter values
include("GEM_main_structs.jl")
include("Setup.jl")

include("AuxiliaryFunc.jl")


# Function to run a single replicate
function run_replicate(params::SimulationParameters, j::Int, i::Int)
    
    # Unpack parameters
    @unpack GEM_ver, t_max, no_species, num_rep, no_columns, no_param, N0, which_par_quant, state_geno_match, state_par_match, param_init, cv_vect, h2_vect, stand_time = params

    # Re-declare and initialize all local variables within the function
    t = 0.0
    N = copy(N0)
    
    # These arrays now hold the data for a single replicate
    pop_slice = zeros(Int, no_species, length(stand_time)) 
    x_slice = fill(NaN, no_columns - 1, length(stand_time), no_species)
    x_var_slice = fill(NaN, no_columns - 1, length(stand_time), no_species)
    Conveniencfunc.jl
    pop_slice[:, 1] .= N0
    
    # Using a named tuple for return values from helper functions would be cleaner
    # For example: (x_dist_init, N_updated) = InitiatePop(...)
    x_dist_init = InitiatePop(N0, which_par_quant, state_geno_match, state_par_match, params, cv_vect, j)
    
    # Store initial state details
    for ii = 1:no_species
        x_slice[:, 1, ii] .= CalcMedian(ii, no_columns, no_param, x_dist_init)
        x_var_slice[1:no_param, 1, ii] .= CalcVar(ii, no_param, x_dist_init)
    end
    
    x_dist = x_dist_init
    time_step_index = 2
    
    while t < t_max && sum(N) > 0
        FindWhoNext = WhoIsNext(x_dist, no_species, no_columns, no_param, N0, state_par_match, state_geno_match)
        params_next = FindWhoNext.params_next # Using named tuple for cleaner access
        genotypes_next = FindWhoNext.genotypes_next
        whosnext = FindWhoNext.whosnext
        
        # We can use a more Julia-idiomatic approach here, e.g., a dictionary for event types
        rates = event_terms(params_next, N)
        
        # PickEvent should return a named tuple
        picked_event = PickEvent(rates, no_species)
        
        # Update state based on the picked event
        if picked_event.event_type == "birth"
            parent_traits = @view x_dist[Int(whosnext[picked_event.species_id]), 2:end]
            new_trait = DrawNewTraits(x_dist, parent_traits, h2_vect, no_param, no_columns, picked_event.species_id, j)
            new_trait_row = vcat(picked_event.species_id, new_trait...)
            x_dist = vcat(x_dist, new_trait_row')
        elseif picked_event.event_type == "death"
            row_to_delete = Int(whosnext[picked_event.species_id])
            x_dist = x_dist[1:size(x_dist, 1) .!= row_to_delete, :]
        end
        
        # Update abundances
        N .= [count(==(species), x_dist[:, 1]) for species in 1:no_species]
        
        # Store data at standardized time steps
        while time_step_index <= length(stand_time) && t > stand_time[time_step_index]
            pop_slice[:, time_step_index] .= N
            for ii in 1:no_species
                x_slice[:, time_step_index, ii] .= CalcMedian(ii, no_columns, no_param, x_dist)
                x_var_slice[1:no_param, time_step_index, ii] .= CalcVar(ii, no_param, x_dist)
            end
            time_step_index += 1
        end
        
        # Advance time
        rate_sum = sum(rates)
        if rate_sum > 0
            t += -log(rand()) / rate_sum
        else
            break
        end
    end
    
    # Return the results for this replicate
    return pop_slice, x_slice, x_var_slice
end

# Main simulation function
function GEM_sim(params::SimulationParameters)
    
    # Unpack parameters
    @unpack GEM_ver, t_max, no_species, num_rep, no_columns, no_param, N0, which_par_quant, state_geno_match, state_par_match, param_init, cv_vect, h2_vect, par_names, geno_names, stand_time = params
    
    # Pre-allocate output arrays
    pop_stand_out_all = zeros(no_species, length(stand_time), num_rep, length(GEM_ver))
    x_stand_out_all = fill(NaN, no_columns - 1, length(stand_time), no_species, num_rep, length(GEM_ver))
    x_var_stand_out_all = fill(NaN, no_columns - 1, length(stand_time), no_species, num_rep, length(GEM_ver))
    
    # Use `@threads` for multi-threading
    Threads.@threads for j = 1:length(GEM_ver)
        pop_stand = zeros(no_species, length(stand_time), num_rep)
        x_stand = fill(NaN, no_columns - 1, length(stand_time), no_species, num_rep)
        x_var_stand = fill(NaN, no_columns - 1, length(stand_time), no_species, num_rep)
        
        Threads.@threads for i = 1:num_rep
            pop_slice, x_slice, x_var_slice = run_replicate(params, j, i)
            pop_stand[:, :, i] .= pop_slice
            x_stand[:, :, :, i] .= x_slice
            x_var_stand[:, :, :, i] .= x_var_slice
        end
        
        pop_stand_out_all[:, :, :, j] .= pop_stand
        x_stand_out_all[:, :, :, :, j] .= x_stand
        x_var_stand_out_all[:, :, :, :, j] .= x_var_stand
    end
    
    # Process results into DataFrames
    pop_df = make_pop_df(pop_stand_out_all, params)
    par_mean_df = make_trait_df(x_stand_out_all, params, par_names, geno_names)
    par_var_df = make_trait_df(x_var_stand_out_all, params, par_names, geno_names)
    
    return SimulationResults(pop_df, par_mean_df, par_var_df)
end

# Helper function to create the population DataFrame
function make_pop_df(data, params)
    @unpack no_species, num_rep, GEM_ver, stand_time = params
    
    long_data = DataFrame(
        time = repeat(stand_time, inner=no_species * num_rep * length(GEM_ver)),
        GEM_ver = repeat(1:length(GEM_ver), inner=no_species * num_rep, outer=length(stand_time)),
        rep = repeat(1:num_rep, inner=no_species, outer=length(stand_time) * length(GEM_ver)),
        species_id = repeat(1:no_species, outer=length(stand_time) * num_rep * length(GEM_ver)),
        population = vec(data)
    )
    
    return long_data
end

# Helper function to create the trait DataFrames
function make_trait_df(data, params, par_names, geno_names)
    @unpack no_columns, no_species, num_rep, GEM_ver, stand_time = params
    
    col_names = vcat(par_names, geno_names)
    
    long_data = DataFrame(
        time = repeat(stand_time, inner=no_columns - 1, outer=no_species * num_rep * length(GEM_ver)),
        GEM_ver = repeat(1:length(GEM_ver), inner=(no_columns - 1) * no_species * num_rep),
        rep = repeat(1:num_rep, inner=(no_columns - 1) * no_species, outer=length(GEM_ver)),
        species_id = repeat(1:no_species, inner=no_columns - 1, outer=num_rep * length(GEM_ver)),
        trait = repeat(col_names, outer=length(stand_time) * no_species * num_rep * length(GEM_ver)),
        value = vec(data)
    )
    
    return long_data
end