# ======================================================================
# MAIN SCRIPT EXECUTION
# This is where we set up the simulation and run it.
# It's good practice to wrap this in a main function.
# ======================================================================

function main()
    # Set seed for reproducibility
    Random.seed!(42)

    # Instantiate the initial population state
    N0 = InitState([5, 1])

    # bd-logistic parameters distribution  
    b_max_mu = 0.8
    b_max_sigma = 0.0
    d_min_mu = 0.4  
    d_min_sigma = 0.0
    #=
    b_s_mu = 1e-2
    b_s_sigma = 0.0
    d_s_mu = 1e-5
    d_s_sigma = 0.0
    =#
    scr_mu = 0.005
    scr_sigma = 0
    fec_mu = 0.05
    fec_sigma = 0
    m_mu = 0.01
    m_sigma = 0

    b_max = rand(LogNormal(log(b_max_mu), b_max_sigma), 1)[1]  # max birth
    d_min = rand(LogNormal(log(d_min_mu), d_min_sigma), 1)[1] # min death
    #b_s = rand(LogNormal(log(b_s_mu), b_s_sigma), 1)[1] # density dependence of birth
    #d_s = rand(LogNormal(log(d_s_mu), d_s_sigma), 1)[1] # density dependence of death
    scr = rand(LogNormal(log(scr_mu), scr_sigma), 1)[1]
    fec = rand(LogNormal(log(fec_mu), fec_sigma), 1)[1]
    m = rand(LogNormal(log(m_mu), m_sigma), 1)[1]
    # calculate initial constant 
    r_max = b_max-d_min
    #K = floor(vec((b_max - d_min)/(b_s + d_s))[1])

    # Instantiate the model parameters struct with the generated values
    model_params = ModelParameters(b_max, d_min, scr, fec, m)

    no_species = length(N0.N) ## also, no_species = size(state_par_match, 1) 
    no_params = fieldcount(ModelParameters)  ## also, size(state_par_match, 2)

    h2 = [0.0 0.0; 0.1 0.1] ## rows: GEM versions, cols: state ID
    cv = [0.0 0.0; 0.2 0.2] ## rows: GEM versions, cols: state ID
    GEM_ver = ["ver1", "ver2"]

    design_choices = DesignChoices(
        h2, # h2_vect
        cv,  # cv_vect
        GEM_ver
    )

    

    # Define the mapping arrays
    state_par_match = [1 1 0 0 0; 0 0 1 1 1] #no_col = param_init, no_row = state
    state_geno_match = [0 0 0 0; 0 0 0 0]
    geno_par_match = [0 0 0 0 0; 0 0 0 0 0]
    which_par_quant = state_par_match - geno_par_match
    no_columns = no_params + 1 + size(state_geno_match, 2) 
    par_names = ["b_max", "d_min", "scr", "fec", "m"]
    geno_names = ["g_1", "g_2", "g_3", "g_4"]

    # Instantiate the new struct
    mappings = SimulationMapping(
        state_par_match,
        state_geno_match,
        geno_par_match,
        par_names,
        geno_names
    )
    # Define time and standardized time steps 
    num_rep = 3
    t_max = 10.0
    min_time_step_to_store = 0.5
    stand_time = range(0, t_max, step = min_time_step_to_store);
    stand_time = collect(stand_time);
    num_time_steps = length(stand_time);


    # Instantiate the new structs with your parameters
    sim_params = SimulationParameters(
        no_species, # no_species::Int
        no_params, # no_param
        no_columns,  # no_columns::Int 
        num_time_steps,
        num_rep, # num_rep::Int
        t_max, # t_max::Float64
        min_time_step_to_store #min_time_step_to_store::Float64 
        )

    # 
    pop_stand_out_all = fill(NaN, no_species, num_time_steps, num_rep, length(GEM_ver));
    x_stand_out_all = fill(NaN, no_columns-1,num_time_steps, no_species,num_rep, length(GEM_ver));
    x_var_stand_out_all = fill(NaN, no_columns-1,num_time_steps, no_species,num_rep, length(GEM_ver));


    sim_output = GEMOutput(
        pop_stand_out_all,
        x_stand_out_all,
        x_var_stand_out_all
        )

#=
    # Instantiate the SimulationParameters struct
        const SIMULATION_PARAMETERS = SimulationParameters(
            GEM_VERSION,
            T_MAX,
            NUM_SPECIES,
            NUM_REPLICATES,
            NUM_COLUMNS,
            NUM_PARAMETERS,
            INITIAL_POPULATION,
            PARAMETER_QUANTITIES,
            GENOTYPE_MATCH,
            PARAMETER_MATCH,
            INITIAL_PARAMETERS,
            CV_VECT,
            H2_VECT,
            PARAMETER_NAMES,
            GENOTYPE_NAMES,
            STANDARDIZED_TIME
        )


        # Call the refactored GEM_sim function
        result = GEM_sim(
            sim_params,
            design_choices,
            model_params,
            initial_state,
            GEM_ver
        )

        println("\nSimulation result: ", result)
=#

end

# To run the script from the command line, you would call `main()`
# at the end of the file.
main()