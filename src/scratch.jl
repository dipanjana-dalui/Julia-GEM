#mapping arrays; row = state; col = param
state_par_match = Array{Int64}([1 1 1 1]) # matching parameters to state
state_geno_match = Array{Int64}([1 0 0]) # matching genotype to state
geno_par_match =   Array{Int64}([0 0 0 0]) # connection b/w parameter and genotype
which_par_quant = state_par_match - geno_par_match
no_columns = no_param + 1 + size(state_geno_match, 2) 
geno_names = ["g_1", "g_2", "g_3"] # genotype name

N = N_init
N0 = copy(N)

param_init = param_vect 
params = copy(param_init) 

h2
cv 
state_geno_match
state_par_match
geno_par_match
which_par_quant

t_max
no_species
no_columns
no_param
num_time_steps
min_time_step_to_store

j=1