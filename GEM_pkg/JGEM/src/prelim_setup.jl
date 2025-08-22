# ======================================================================
# MAIN SCRIPT EXECUTION
# This is where we set up the simulation and run it.
# It's good practice to wrap this in a main function.
# ======================================================================

function main()
    # Set seed for reproducibility
    Random.seed!(42)

    # Instantiate the new structs with your parameters
    sim_params = SimulationParameters(
        3, # num_rep
        10.0, # t_max
        0.5 # min_time_step_to_store
    )

    design_choices = DesignChoices(
        [0.0 0.0; 0.1 0.1], # h2_vect
        [0.0 0.0; 0.2 0.2]  # cv_vect
    )

    # Generate model parameters using the same logic as your original script
    b_max = rand(LogNormal(log(0.8), 0.0), 1)[1]
    d_min = rand(LogNormal(log(0.4), 0.0), 1)[1]
    scr = rand(LogNormal(log(0.005), 0.0), 1)[1]
    fec = rand(LogNormal(log(0.05), 0.0), 1)[1]
    m = rand(LogNormal(log(0.01), 0.0), 1)[1]

    # Instantiate the model parameters struct with the generated values
    model_params = ModelParameters(b_max, d_min, scr, fec, m)

    # Instantiate the initial population state
    initial_state = PopulationState([5, 1])

    # Other non-struct parameters
    GEM_ver = ["ver1", "ver2"]

    # Define the mapping arrays
    state_par_match = [1 1 0 0 0; 0 0 1 1 1]
    state_geno_match = [0 0 0 0; 0 0 0 0]
    geno_par_match = [0 0 0 0 0; 0 0 0 0 0]
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
end

# To run the script from the command line, you would call `main()`
# at the end of the file.
main()