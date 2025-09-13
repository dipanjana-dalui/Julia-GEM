"""
This file confures the model and the simulation.
"""

# ======================================================================
#                    INITIAL STATE AND PARAMETERS
# ======================================================================

# Set seed for reproducibility
Random.seed!(42)

N_init = [10]

# bd-logistic parameters distribution mu and sigma
# maximum birth
b_max_mu = 4.0
b_max_sigma = 0.1

# minimum death
d_min_mu = 1.0
d_min_sigma = 0.0

# density dependence of birth
b_s_mu = 0.0012
b_s_sigma = 0.0

# density dependence of death
d_s_mu = 1e-5
d_s_sigma = 0.0


# randomly draw one sample from the lognormal distrinution
# (see PickTrait.jl for MU and SIGMA transformations )
include("functions/PickTrait.jl") # load sampling function

b_max = PickTrait(b_max_mu, b_max_sigma)  # max birth
d_min = PickTrait(d_min_mu, d_min_sigma) # min death
b_s = PickTrait(b_s_mu, b_s_sigma) # density dependence of birth
d_s = PickTrait(d_s_mu, d_s_sigma) # density dependence of death

# calculate initial constant 
r_max = b_max-d_min
K = floor(((b_max - d_min)/(b_s + d_s)))

param_vect = [b_max, d_min, b_s, d_s] # parameter vector 
par_names = ["b_max", "d_min", "b_s", "d_s"] # parameter names
no_species = length(N_init) 
no_param = length(param_vect) 

# overall simulation design
GEM_ver = ["ver1", "ver2"] # number of GEM version
#h2 = [0.0 0.0; 0.1 0.1] ## rows: GEM versions, cols: state ID
#cv = [0.0 0.0; 0.2 0.2] ## rows: GEM versions, cols: state ID
h2_vect = [0.0;
        0.2] # narrow sense heritability
#h2 = reshape(h2_vect, length(GEM_ver), no_species)

#cv_vect = [0.0;
#        0.2] # coefficient of variation
#cv = reshape(cv_vect, length(GEM_ver), no_species)

cv_vect = cat([ 0.0 0.0 0.0 0.0;
                0.0 0.0 0.0 0.0],
              [ 0.3 0.1 0.0 0.0;
                0.0 0.0 0.1 0.0], dims=3)

#= 4×1×2 Array{Float64, 3}: row: state; col = param; stack = GEM ver
[:, :, 1] =
 0.0  0.0  0.0  0.0
[:, :, 2] =
 0.2  0.1  0.0  0.0
 =#

# mapping arrays; row = state; col = param
state_par_match = Array{Int64}([1 1 1 1]) # matching parameters to state
state_geno_match = Array{Int64}([1 0 0 0]) # matching genotype to state
geno_par_match = Array{Int64}([0 0 0 0]) # connection b/w parameter and genotype
which_par_quant = state_par_match - geno_par_match
no_columns = no_param + 1 + size(state_geno_match, 2) 
geno_names = ["g_1", "g_2", "g_3", "g_4"] # genotype name

# replicate and time
num_rep = 2 # number of replicates
t_max = 5.0 # maximum time 
min_time_step_to_store = 0.5 #
stand_time = range(0, t_max, step = min_time_step_to_store) 
stand_time = collect(stand_time)
num_time_steps = length(stand_time)


# storage containers
pop_stand_out_all = fill(NaN, no_species, num_time_steps, num_rep, length(GEM_ver))
x_stand_out_all = fill(NaN, no_columns-1,num_time_steps, no_species,num_rep, length(GEM_ver))
x_var_stand_out_all = fill(NaN, no_columns-1,num_time_steps, no_species,num_rep, length(GEM_ver))


# ======================================================================
#                             INSTANTIATE 
# Only make changes to this block if you make changes to the variable 
# names that are used to instantiate the struct. 
# ======================================================================
""" 1. Instantiate the initial population state """
N0 = InitState(N_init) # Vector{Int}

""" 2. Instantiate ModelParVector """
model_par_vect = ModelParVector(
    param_vect # Vector{T}
)

""" 3. Instantiate DesignChoices """
design_choices = DesignChoice(
    h2, # Matrix{Float64}
    cv,  # Matrix{Float64}
    GEM_ver # Vector{String}
)

""" 4. Instantiate SimulationMap """
mappings = SimulationMap(
    state_par_match, # Matrix{Int}
    state_geno_match, # Matrix{Int}
    geno_par_match, # Matrix{Int}
    which_par_quant, # Matrix{Int}
    par_names, # Vector{String}
    geno_names # Vector{String}
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
    pop_stand_out_all, # Array{Float64, 4}
    x_stand_out_all, # Array{Float64, 5}
    x_var_stand_out_all # Array{Float64, 5}
    )
    

    