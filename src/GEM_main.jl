"""
    1. Function to run one replication 
"""

# Function to run a single replicate
function run_replicate(init_state::InitState,
                        #mod_par::ModelParameters,
                        mod_par_vect::ModelParVector,
                        dc::DesignChoice,
                        sim_map::SimulationMap,
                        sim_par::SimulationParameters,
                        j::Int)
    @unpack N = init_state
    #@unpack b_max, d_min, scr, fec, m = mod_par
    #@unpack b_max, d_min, b_s, d_s = mod_par
    @unpack param_init = mod_par_vect
    @unpack h2, cv = dc
    @unpack state_geno_match, state_par_match, geno_par_match, which_par_quant = sim_map
    @unpack t_max, no_species, no_columns, no_param, num_time_steps,min_time_step_to_store  = sim_par    
    #j=1

    t = 0.0
    N0 = copy(N)
    N = copy(N)
    #params = [b_max, d_min, scr, fec, m] 
    #params = [b_max, d_min, b_s, d_s]
    params = copy(param_init) 
    @show params
    
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
        state_geno_match, 
        state_par_match, 
        init_comm_mat,
        params, 
        cv, 
        j)
    
    # Store initial state details
    for ii = 1:no_species
        x_slice[:, 1, ii] = CalcMedian(ii, no_columns, no_param, x_dist_init)
        #x_slice[:, 1, ii] .= CalcMean(ii, no_columns, no_param, x_dist_init)
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
                 """Instantiate"""     
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

"""parallel sim function"""

# Main simulation function
function GEM_sim(init_state::InitState,
                        #mod_par::ModelParameters,
                        mod_par_vect::ModelParVector,
                        dc::DesignChoice,
                        sim_map::SimulationMap,
                        sim_par::SimulationParameters,
                        sim_op::GEMOutput)
    
    @unpack N = init_state
    #@unpack b_max, d_min, scr, fec, m = mod_par
    #@unpack b_max, d_min, b_s, d_s = mod_par
    @unpack param_init = mod_par_vect
    @unpack h2, cv, GEM_ver = dc
    @unpack state_geno_match, state_par_match, geno_par_match, which_par_quant = sim_map
    @unpack t_max, no_species, no_columns, no_param, num_time_steps,min_time_step_to_store  = sim_par    
    @unpack pop_stand_out_all, x_stand_out_all, x_var_stand_out_all = sim_op
    
    # Use `@threads` for multi-threading
    Threads.@threads for j = 1:length(GEM_ver)
        pop_stand = zeros(no_species, num_time_steps, num_rep)
        x_stand = fill(NaN, no_columns - 1, num_time_steps, no_species, num_rep)
        x_var_stand = fill(NaN, no_columns - 1, num_time_steps, no_species, num_rep)
        
        Threads.@threads for i = 1:num_rep
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
    
   
    return (pop_stand_out_all, x_stand_out_all,x_var_stand_out_all)
end

#===================================================================#

function make_pop_df_long(sim_output::GEMOutput,
                        sim_par::SimulationParameters, 
                        dc::DesignChoice )
    @unpack pop_stand_out_all = sim_output
    @unpack t_max, no_species, no_columns, no_param, num_time_steps,min_time_step_to_store, num_rep  = sim_par    
    @unpack GEM_ver = dc


    pop_out = vec(pop_stand_out_all)
    pop_out_mat = reshape(pop_out, no_species, num_time_steps, num_rep, length(GEM_ver))

    pop_out_gem_v_store = Vector{DataFrame}(undef, length(GEM_ver)) 
    pop_out_spp_store = Vector{DataFrame}(undef, no_species)

    for k = 1:no_species
        for i = 1:length(GEM_ver)
            temp = DataFrame(hcat(stand_time, fill(i,num_time_steps),
                            fill(k,num_time_steps)),:auto)
            rename!(temp, :x1 => :time, :x2 => :GEM_ver, :x3 => :state_ID)
            pop_out_temp = DataFrame(pop_out_mat[k,:,:,i], :auto)
            pop_out_temp = hcat(temp, pop_out_temp, makeunique=true)
            pop_out_gem_v_store[i] = pop_out_temp
        end
        pop_out_spp_store[k] = vcat(pop_out_gem_v_store...)
    end
    return pop_df = vcat(pop_out_spp_store...)
end

# =============================================================

function make_trait_df_long(sim_output::GEMOutput,
                                sim_par::SimulationParameters, 
                                dc::DesignChoice,
                                sim_map::SimulationMap )
                                
    @unpack x_stand_out_all,x_var_stand_out_all = sim_output
    @unpack t_max, no_species, no_columns, no_param, num_time_steps,min_time_step_to_store, num_rep  = sim_par    
    @unpack GEM_ver = dc
    @unpack par_names, geno_names = sim_map

    x_out = vec(x_stand_out_all)
    x_out_var = vec(x_var_stand_out_all)

    x_out_mat = reshape(x_out, no_columns-1,num_time_steps, no_species,num_rep, length(GEM_ver))
    x_out_var_mat = reshape(x_out_var, no_columns-1,num_time_steps, no_species,num_rep, length(GEM_ver))
 
    col_names = vcat(par_names, geno_names)
    
    x_out_spp_store = Vector{DataFrame}(undef, no_species)
    x_out_gem_v_store = Vector{DataFrame}(undef, length(GEM_ver)) 
    x_out_rep_store = Vector{DataFrame}(undef, num_rep)

    x_out_var_spp_store = Vector{DataFrame}(undef, no_species)
    x_out_var_gem_v_store = Vector{DataFrame}(undef, length(GEM_ver)) 
    x_out_var_rep_store = Vector{DataFrame}(undef, num_rep)


    for k = 1:no_species
        #k = 2
        for i = 1:length(GEM_ver)
            #i = 1
            for j = 1:num_rep
                col_df = DataFrame(hcat(stand_time,
                        fill(j,num_time_steps),
                        fill(i,num_time_steps),
                        fill(k,num_time_steps)),:auto)
                rename!(col_df, :x1 => :time, :x2 => :rep, :x3 => :GEM_ver, :x4 => :state_ID)
                x_temp = x_out_mat[:,:,k,j,i]
                x_temp = DataFrame(transpose(x_temp), :auto)
                rename!(x_temp, col_names)
                x_out_dat = hcat(col_df, x_temp)
                x_out_rep_store[j] = x_out_dat
            end
            x_out_gem_v_store[i] = vcat(x_out_rep_store...)
        end
        x_out_spp_store[k] = vcat(x_out_gem_v_store...)
    end

    for k = 1:no_species
        for i = 1:length(GEM_ver)
            for j = 1:num_rep
                col_df = DataFrame(hcat(stand_time,
                        fill(j,num_time_steps),
                        fill(i,num_time_steps),
                        fill(k,num_time_steps)),:auto)
                rename!(col_df, :x1 => :time, :x2 => :rep, :x3 => :GEM_ver, :x4 => :state_ID)
                x_temp = x_out_var_mat[:,:,k,j,i]
                x_temp = DataFrame(transpose(x_temp), :auto)
                rename!(x_temp, col_names)
                x_out_var_dat = hcat(col_df, x_temp)
                x_out_var_rep_store[j] = x_out_var_dat
            end
            x_out_var_gem_v_store[i] = vcat(x_out_var_rep_store...)
        end
        x_out_var_spp_store[k] = vcat(x_out_var_gem_v_store...)
    end
    
    par_mean_df = vcat(x_out_spp_store...)
    par_var_df = vcat(x_out_var_spp_store...)

    return (trait_mean_df=deci_threshold(par_mean_df), trait_var_df=deci_threshold(par_var_df))

end

