"""
    1. Function to run one replication 
"""

# Function to run a single replicate
function run_replicate(init_state::InitState,
                        #mod_par::ModelParameters,
                        mod_par_vect::ModelParVector,
                        dc::DesignChoices,
                        sim_map::SimulationMap,
                        sim_par::SimulationParameters,
                        j::Int)
    @unpack N = init_state
    @unpack param_init = mod_par_vect
    @unpack h2, cv = dc
    @unpack state_geno_match, state_par_match, geno_par_match, which_par_quant = sim_map
    @unpack t_max, no_species, no_columns, no_param, num_time_steps,min_time_step_to_store  = sim_par    

    # make a copy 
    t = 0.0
    N0 = copy(N)
    N = copy(N)
    params = copy(param_init) 
    #@show params
    #@show N0
    
    # some internal container
    init_comm_mat =  Array{Float64}(fill(NaN, Int(sum(N0)), no_columns))
    pop_slice = Array{Int}(fill(0, no_species, num_time_steps))
    x_slice = fill(NaN, no_columns-1, num_time_steps, no_species)
    x_var_slice = fill(NaN, no_columns-1, num_time_steps, no_species)

    # store details of initial state 
    pop_slice[:,1] .= N0 #first col/time step gets the initial pop 
    
    """ Func Initiate Population """
    x_dist_init = InitiatePop(
        N0, 
        which_par_quant, 
        geno_par_match,
        state_geno_match, 
        state_par_match, 
        init_comm_mat,
        params, 
        cv, 
        j)
    # Store initial state details
    for ii = 1:no_species
        x_slice[:, 1, ii] = CalcMedian(ii, no_columns, no_param, x_dist_init)
        x_var_slice[1:no_param, 1, ii] = CalcVar(ii, no_param, x_dist_init)
    end

    # count up each individual for all states after the first sampling
    for jj = 1:length(N)
        x = init_comm_mat[:,1] #extract first col
        N[jj] = count(.==(jj), x)
    end
    
    x_dist = x_dist_init
    time_step_index = 2
    time_step = stand_time[time_step_index]
    
    while t < t_max && sum(N) > 0
        """ Func WhoIsNext """
        FindWhoNext = WhoIsNext(x_dist, no_species, no_columns, no_param, N0, state_par_match, state_geno_match)
        param_next = FindWhoNext.param_next # FindWhoNext[1]; Using named tuple for cleaner access
        genotype_next = FindWhoNext.genotype_next # FindWhoNext[2]
        whosnext = FindWhoNext.whosnext # FindWhoNext[3]
        
        """ Func Event Terms """
        terms = collect(Event_Terms(param_next, N))
        
        """ Func Pick Event """
        picked_event = PickEvent(terms, no_species)
        c_sum = picked_event.c_sum # PickedEvent[1]
        row = picked_event.row # PickedEvent[2]
        col = picked_event.col # PickedEvent[3]
        #@show PickedEvent
        if row == 1 #&& col == 1
            parent_traits = x_dist[Int(whosnext[col]), 2:no_columns] 
            #@show parent_traits   
                 """ Draw New Trait """     
            new_trait = DrawNewTraits(x_dist,parent_traits,h2,no_param,no_columns,col, j)
                    
            new_trait_row = hcat(col, new_trait)
            new_trait_row = hcat(new_trait_row[1],new_trait_row[2][1],new_trait_row[2][2])
            x_dist = vcat(x_dist, new_trait_row) 
        
        elseif row == 2 #&& col ==1 # death 
            # delete the individual by making a new copy of the matrix
            # without the row 
            x_dist = x_dist[1:size(x_dist, 1) .!= Int(whosnext[col]), :]
        end

        # Update abundances 
        for jj in 1:no_species
            N[jj] = sum(x_dist[:,1].== jj)
        end
        
        # while loop
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
        time_advance = exp(-1/c_sum[end])/(c_sum[end])
        if !isnan(time_advance) && time_advance > 0 
            t = t + time_advance
        else
            break
            println("Time advance error. Stopped at time:\nT $T")
        end

        # store the last value of the replicate
        pop_slice[1:no_species, time_step_index] = N
        for ii = 1:no_species
            x_slice[:,time_step_index,ii] = CalcMedian(ii,no_columns,no_param,x_dist);
            x_var_slice[1:no_param,time_step_index,ii] = CalcVar(ii, no_param,x_dist)
        end
        
    end
    
    # Return the results for this replicate
    return (pop_time_series=pop_slice, trait_mom1 = x_slice, trait_mom2 = x_var_slice)
end

# ==================================================================

"""
    2. parallel GEM simulation function
"""

# Main simulation function
function GEM_sim(init_state::InitState,
                        mod_par_vect::ModelParVector,
                        dc::DesignChoices,
                        sim_map::SimulationMap,
                        sim_par::SimulationParameters,
                        sim_op::GEMOutput)
    
    @unpack N = init_state
    @unpack param_init = mod_par_vect
    @unpack h2, cv, GEM_ver = dc
    @unpack state_geno_match, state_par_match, geno_par_match, which_par_quant = sim_map
    @unpack num_rep, t_max, no_species, no_columns, no_param, num_time_steps,min_time_step_to_store  = sim_par    
    @unpack pop_stand_out_all, x_stand_out_all, x_var_stand_out_all = sim_op
    
    # Use `@threads` for multi-threading
    for j = 1:length(GEM_ver) # loop through the GEM versions

        # some internal containers
        pop_stand = zeros(no_species, num_time_steps, num_rep)
        x_stand = fill(NaN, no_columns - 1, num_time_steps, no_species, num_rep)
        x_var_stand = fill(NaN, no_columns - 1, num_time_steps, no_species, num_rep)
        
        #Threads.@threads 
        for i = 1:num_rep # loop through the replicates
            @show j
            @show i
            @show Threads.threadid()
            pop_slice, x_slice, x_var_slice = run_replicate(init_state, mod_par_vect, design_choices, mappings, 
                sim_params,j)
            pop_stand[:, :, i] .= pop_slice
            x_stand[:, :, :, i] .= x_slice
            x_var_stand[:, :, :, i] .= x_var_slice
        end
        
        pop_stand_out_all[:, :, :, j] .= pop_stand
        x_stand_out_all[:, :, :, :, j] .= x_stand
        x_var_stand_out_all[:, :, :, :, j] .= x_var_stand
    end
    
    # turn the multidimensional output arrays into long dataframes
    pop_out = make_pop_df_long(sim_output, sim_params, design_choices) # population time series dataframe
    trait_out = make_trait_df_long(sim_output, sim_params, design_choices, mappings) # trait mean and var time seroes dataframes 
    # trait_out has two dataframes that can be accessed with named tuples: trait_out.mean and trait_out.var 

    return (pop_df = pop_out, trait_df = trait_out)
end

#===================================================================#

