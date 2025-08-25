# scractch file for all structs
"""
InitiatePop
"""
struct InitiatePopArg
N0::Vector{Int64}
which_par_quant::Matrix{Int}
state_geno_match::Matrix{Int}, 
state_par_match::Matrix{Int}, 
init_comm_mat::Matrix{Float64}, 
params::Vector{Float64}, 
cv_vect::Matrix{Float64}, 
j::Int
end

"""
PickIndiv
"""
# n/a

"""
CalcAvgFreq
"""
struct CalcStatArg
    ii::Int
    no_columns::Int
    no_param::Int
    x_dist::Matrix{Float64}
end

"""
WhiIsNext
"""
struct WINArg
    x_dist::Matrix{Float64}, 
    no_species::Int 
    no_columns::Int 
    no_param::Int
    N0::Vector{Int} 
    state_par_match::Matrix{Int}, 
    state_geno_match::Matrix{Int}
end

"""
PickEvent
"""
struct PEArg
    terms::Vector{Float64}
    no_species::Int
end

"""
DrawNewTraits
"""
struct DNTArgs
    x_dist::Matrix{Float64}, 
    parent_traits::Vector{Float64}, 
	h2_vect::Matrix{Float64}, 
    no_param::Int, 
    no_columns::Int, 
    col::Int, 
    j::Int
end

################################################
################################################

