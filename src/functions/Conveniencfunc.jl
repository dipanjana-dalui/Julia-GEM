## Function to convert 0s to NaN in matrix
function zero_to_nan(Bool_mat::Matrix{Int64})
	return [x==0 ? NaN : x for x in Bool_mat]
end

# =============================================================

## function to set a threshold for zero
function deci_threshold(df::DataFrame)

	df_mod = mapcols(col -> round.(col, digits=8), df)
	return df_mod
end

# =============================================================


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
    
    par_median_df = vcat(x_out_spp_store...)
    par_var_df = vcat(x_out_var_spp_store...)

    return (median=deci_threshold(par_median_df), var=deci_threshold(par_var_df))

end

