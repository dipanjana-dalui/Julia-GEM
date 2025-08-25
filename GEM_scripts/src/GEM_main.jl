"""

"""

# Include the files with the struct definitions and parameter values



# Function to run a single replicate
function run_replicate(params::SimulationParameters, j::Int, i::Int)
    
    # Unpack parameters
    @unpack GEM_ver, t_max, no_species, num_rep, no_columns, no_param, N0, which_par_quant, state_geno_match, state_par_match, param_init, cv_vect, h2_vect, stand_time = params

    #########
N0 = [5,1]
no_species = 2
no_columns=10
num_time_steps = 20

    # Re-declare and initialize all local variables within the function
    t = 0.0
    N0 = copy(N0)
    N = copy(N0)
    params = copy(param_init) 

    init_comm_mat =  Array{Float64}(fill(NaN, Int(sum(N0)), no_columns))
    pop_slice = Array{Int}(fill(0, no_species, num_time_steps))
    x_slice = fill(NaN, no_columns-1, num_time_steps, no_species)
    x_var_slice = fill(NaN, no_columns-1, num_time_steps, no_species)

    # store details of initial state 
    pop_slice[:,1] .= N0 #first col/time step gets the initial pop 
    
    """ instantiate the struct for Initiate Population func here """
    x_dist_init = InitiatePop(N0, which_par_quant, state_geno_match, state_par_match, params, cv_vect, j)
    
    # Store initial state details
    for ii = 1:no_species
        x_slice[:, 1, ii] .= CalcMedian(ii, no_columns, no_param, x_dist_init)
        #x_slice[:, 1, ii] .= CalcMean(ii, no_columns, no_param, x_dist_init)
        x_var_slice[1:no_param, 1, ii] .= CalcVar(ii, no_param, x_dist_init)
    end
    """Do we need this in a function for one rep?
                # count up each individual for all states after the first sampling
            for jj = 1:length(N)
                x = init_comm_mat[:,1] #extract first col
                N[jj] = count(.==(jj), x)
            end

    """
    x_dist = x_dist_init
    time_step_index = 2
    time_step = stand_time[time_step_index]
    
    while t < t_max && sum(N) > 0
        """ Instantiate for FindWhoNext """
        FindWhoNext = WhoIsNext(x_dist, no_species, no_columns, no_param, N0, state_par_match, state_geno_match)
        params_next = FindWhoNext.params_next # FindWhoNext[1]; Using named tuple for cleaner access
        genotypes_next = FindWhoNext.genotypes_next # FindWhoNext[2]
        whosnext = FindWhoNext.whosnext # FindWhoNext[3]
        
                """Instantiate"""
        terms = collect(event_terms(params_next, N))
        
        """Instantiate"""
        picked_event = PickEvent(rates, no_species)
        c_sum = PickedEvent.c_sum # PickedEvent[1]
        row = PickedEvent.row # PickedEvent[2]
        col = PickedEvent.col # PickedEvent[3]
        #@show PickedEvent
        if row == 1 #&& col == 1
            parent_traits = x_dist[Int(whosnext[col]), 2:no_columns] 
            #@show parent_traits   
                 """Instantiate"""     
            new_trait = DrawNewTraits(x_dist,parent_traits,h2_vect,no_param,no_columns,col, j)
                    
            new_trait_row = hcat(col, new_trait)
            new_trait_row = hcat(new_trait_row[1],new_trait_row[2][1],new_trait_row[2][2])
            x_dist = vcat(x_dist, new_trait_row) 

        elseif row == 2 #&& col ==1 # death 
            # delete the individual by making a new copy of the matrix
            # without the row 
            x_dist = x_dist[1:size(x_dist, 1) .!= Int(whosnext[col]), :]
                nd

        # Update abundances
        ## UPDATE ABUNDANCES 
                for jj in 1:no_species
                    N[jj] = sum(x_dist[:,1].== jj)
                end
        
        
        while t > time_step
                    pop_slice[:,time_step_index] .= N  # assign current values to sliced standard times
                    for ii in 1:no_species
                        x_slice[:,time_step_index,ii] = CalcMedian(ii,no_columns,no_param,x_dist)
                        x_var_slice[1:no_param,time_step_index,ii] = CalcVar(ii, no_param, x_dist)
                    end
                    time_step_index +=  1 # advance to next standardized time
                    time_step = stand_time[time_step_index]
                end
        
        # Advance time
        if !isnan(time_advance) && time_advance > 0 ##
            t = t + time_advance
        else
            break
            println("Time advance error. Stopped at time:\nT $T")
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