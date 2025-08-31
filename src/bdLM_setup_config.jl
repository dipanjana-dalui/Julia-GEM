# ======================================================================
# MAIN SCRIPT EXECUTION
# This is where we set up the simulation and run it.
# It's good practice to wrap this in a main function.
# ======================================================================

#function set_init_conditions()
    # Set seed for reproducibility
    Random.seed!(42)

    N_init = [10]

    # bd-logistic parameters distribution  
    b_max_mu = 4.0
    b_max_sigma = 0.0

    b_s_mu = 0.0012
    b_s_sigma = 0.0

    d_s_mu = 1e-5
    d_s_sigma = 0.0

    d_min_mu = 1
    d_min_sigma = 0.0

    b_max = rand(LogNormal(log(b_max_mu), b_max_sigma), 1)[1]  # max birth
    d_min = rand(LogNormal(log(d_min_mu), d_min_sigma), 1)[1] # min death
    b_s = rand(LogNormal(log(b_s_mu), b_s_sigma), 1)[1] # density dependence of birth
    d_s = rand(LogNormal(log(d_s_mu), d_s_sigma), 1)[1] # density dependence of death
    # calculate initial constant 
    r_max = b_max-d_min
    K = floor(((b_max - d_min)/(b_s + d_s)))
    
    param_vect = [b_max, d_min, b_s, d_s] 
    no_species = length(N_init) 
    no_param = fieldcount(ModelParameters) 
    
    GEM_ver = ["ver1", "ver2"]
    #h2 = [0.0 0.0; 0.1 0.1] ## rows: GEM versions, cols: state ID
    #cv = [0.0 0.0; 0.2 0.2] ## rows: GEM versions, cols: state ID
    h2_vect = [0.0;
            0.2] 
    h2 = reshape(h2_vect, length(GEM_ver), no_species)

    cv_vect = [0.0;
            0.2]
    cv = reshape(cv_vect, length(GEM_ver), no_species)

    # Define the mapping arrays
    state_par_match = Array{Int64}([1 1 1 1]) # matching parameters to state
    state_geno_match = Array{Int64}([0 0 0 0]) # matching genotype to state
    geno_par_match = Array{Int64}([0 0 0 0]) # connection b/w parameter and genotype
    which_par_quant = state_par_match - geno_par_match
    no_columns = no_param + 1 + size(state_geno_match, 2) 
    par_names = ["b_max", "d_min", "b_s", "d_s"]
    geno_names = ["g_1", "g_2", "g_3", "g_4"]

    num_rep = 3
    t_max = 10.0
    min_time_step_to_store = 0.5
    stand_time = range(0, t_max, step = min_time_step_to_store)
    stand_time = collect(stand_time)
    num_time_steps = length(stand_time)


    # 
    pop_stand_out_all = fill(NaN, no_species, num_time_steps, num_rep, length(GEM_ver))
    x_stand_out_all = fill(NaN, no_columns-1,num_time_steps, no_species,num_rep, length(GEM_ver))
    x_var_stand_out_all = fill(NaN, no_columns-1,num_time_steps, no_species,num_rep, length(GEM_ver))

    """ 1. Instantiate the initial population state """
    N0 = InitState(N_init)
    
    """ 2. Instantiate the model parameters struct with the generated values """
    model_params = ModelParameters(b_max, d_min, b_s, d_s)

    """ 3. Instantiate DesignChoices """
    design_choices = DesignChoice(
        h2, # h2_vect::Matrix{Float64}
        cv,  # cv_vect::Matrix{Float64}
        GEM_ver # GEM_ver::Vector{String}
    )

    """ 4. Instantiate SimulationMap """
    mappings = SimulationMap(
        state_par_match,
        state_geno_match,
        geno_par_match,
        which_par_quant,
        par_names,
        geno_names
    )

    """  5. Instantiate SimulationParameters """
    sim_params = SimulationParameters(
        no_species, # Int
        no_param, # Int
        no_columns,  # Int 
        num_time_steps, #Int
        num_rep, # Int
        t_max, # Float64
        min_time_step_to_store # Float64 
        )
    """ 6. Instantiate output container """
    sim_output = GEMOutput(
        pop_stand_out_all,
        x_stand_out_all,
        x_var_stand_out_all
        )
        
    """ 7. Instantiate ModelParVector """
    model_par_vect = ModelParVector(
        param_vect
    )

        