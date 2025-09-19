
using Distributions
using Random
include("functions/PickTrait.jl") # load sampling function

# Set seed for reproducibility
Random.seed!(42)

N_init = [5, 1]


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
scr_sigma = 0.0
fec_mu = 0.05
fec_sigma = 0.0
m_mu = 0.01
m_sigma = 0.0

# randomly draw one sample from the lognormal distrinution
# (see PickTrait.jl for MU and SIGMA transformations )

b_max = PickTrait(b_max_mu, b_max_sigma)  # max birth
d_min = PickTrait(d_min_mu, d_min_sigma) # min death
#b_s = rand(LogNormal(log(b_s_mu), b_s_sigma), 1)[1] # density dependence of birth
#d_s = rand(LogNormal(log(d_s_mu), d_s_sigma), 1)[1] # density dependence of death
scr = PickTrait(scr_mu, scr_sigma)
fec = PickTrait(fec_mu, fec_sigma)
m = PickTrait(m_mu, m_sigma)
# calculate initial constant 
r_max = b_max-d_min
#K = floor(vec((b_max - d_min)/(b_s + d_s))[1])

param_vect = [b_max, d_min, scr, fec, m]
par_names = ["b_max", "d_min", "scr", "fec", "m"]
no_species = length(N_init) 
no_param = length(param_vect) 

# Define the mapping arrays
# nrow = state, ncol = param
state_par_match = [1 1 0 0 0; 0 0 1 1 1]
# nrow = state, ncol = genotype 
state_geno_match = [1 0 0 0; 0 0 1 0] # startng genotype

no_columns = no_param + 1 + size(state_geno_match, 2) 
geno_names = ["g_1", "g_2", "g_3", "g_4"]

GEM_ver = ["ver1", "ver2"]
# nrow: state ID, ncol: GEM versions
h2 = [ 0.0 0.2 ;
       0.0 0.2 ]# narrow sense heritability

# cv = nrow:state ID, ncol:length(param), stack:GEM ver}
"""
Note: The first stack is for GEM ver 1; typically reserved for "no-evolution". All elements are set to 0.0
In stack 2, set cv value for parameters corresponding to each state. 
You can mirror the dimensions of state_parameter_match matrix defined above. 
1 -> cv value
0 -> n/a for this state
"""

cv = cat([ 0.0 0.0 0.0 0.0 0.0 ; 0.0 0.0 0.0 0.0 0.0],
         [ 0.3 0.1 0.0 0.0 0.0 ; 0.0 0.0 0.1 0.1 0.05], dims=3)


num_rep = 2
t_max = 5.0 #upwards of 6 the events times get very very small
min_time_step_to_store = 0.5
stand_time = range(0, t_max, step = min_time_step_to_store)
stand_time = collect(stand_time)
num_time_steps = length(stand_time)


# 
pop_stand_out_all = fill(NaN, no_species, num_time_steps, num_rep, length(GEM_ver))
x_stand_out_all = fill(NaN, no_columns-1,num_time_steps, no_species,num_rep, length(GEM_ver))
x_var_stand_out_all = fill(NaN, no_columns-1,num_time_steps, no_species,num_rep, length(GEM_ver))


# ======================================================================
#                             INSTANTIATE 
# Only make changes to this block if you make changes to the variable 
# names that are used to instantiate the struct. 
# ======================================================================
""" 1. Instantiate the initial population state """
N0 = InitStates(N_init) # Vector{Int}

""" 2. Instantiate ModelParVector """
model_par_vect = ModelParVector(
    param_vect # Vector{T}
)

""" 3. Instantiate DesignChoices """
design_choices = DesignChoices(
    h2, # Matrix{Float64}
    cv,  # Matrix{Float64}
    GEM_ver # Vector{String}
)

""" 4. Instantiate SimulationMap """
mappings = SimulationMaps(
    state_par_match, # Matrix{Int}
    state_geno_match, # Matrix{Int}
    #geno_par_match, # Matrix{Int}
    #which_par_quant, # Matrix{Int}
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
    

    