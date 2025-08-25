##############################################
#		  FUNCTION DRAW NEW TRAITS           #
##############################################
#function DrawNewTraits(x_dist::Matrix{Float64}, parent_traits::Vector{Float64}, 
#						h2_vect::Matrix{Float64}, no_param::Int64, no_columns::Int64, col::Int64, j::Int64)
function DrawNewTraits(p::DNTArgs)

    @unpack x_dist, parent_traits, h2_vect, no_param, no_columns, col, j = p

	pop_size = size(x_dist[x_dist[:, 1] .== col, 2:no_param+1])
	# QUANTITATIVE TRAITS
	pop_mean = mapslices(mean ∘ skipmissing, x_dist[x_dist[:, 1] .== col, 2:no_param+1], dims=1)

	if pop_size[1] == 1
		pop_var = Matrix{Float64}(fill(0, pop_size[1], pop_size[2]))
	else
		pop_var = mapslices(var ∘ skipmissing, x_dist[x_dist[:, 1] .== col, 2:no_param+1], dims=1)
	end
	pop_var = round.(pop_var, digits=6)
	pop_stdev = round.(sqrt.(pop_var), digits=6)

	exp_offspring_traits = h2_vect[j,col] .* reshape(parent_traits[1:no_param], size(pop_mean)[1], size(pop_mean)[2]) .+ (1-h2_vect[j,col]) .* pop_mean		
	sigma_offspring_traits = sqrt(1-(h2_vect[j,col])^2)*pop_stdev
	# MU and SIGMA for lognormal
	MU = log.(exp_offspring_traits.^2 ./ sqrt.(sigma_offspring_traits.^2 .+ exp_offspring_traits.^2)) 	
	SIGMA = sqrt.(log.(sigma_offspring_traits.^2 ./ exp_offspring_traits.^2 .+ 1))

	offspring_traits = fill(NaN, 1, length(exp_offspring_traits))
	for i = 1:length(exp_offspring_traits)
		if !isnan(MU[i]) && !isnan(SIGMA[i])
			if !iszero(SIGMA[i])
				offspring_traits[1,i] = rand(LogNormal(MU[i], SIGMA[i])) # pull traits when sigma != 0
			else
				offspring_traits[1,i] = exp_offspring_traits[i]
			end 
		end
	end
	
	## JPD: edit this please
	offspring_genotypes = collect(transpose(parent_traits[no_param+1:no_columns-1]))
	##### what happens to genotype?

	return (offspring_traits=offspring_traits, offspring_genotype=offspring_genotypes)
end

"""
MoJuWo:
split into multiple functions
- get_pop_stats
	retuns MU and SIGMA
- DrawNewTraits
	- the evaluation inside is fairly similar,
	with only stylistic differences. 
	- best to not change output dimensions just yet.
	- DON'T FIX WHAT AIN'T BROKEN.
	- multi-threading could be useful (using Base.Threads).


# Use broadcasting and a more robust way to handle single-individual variance
function get_pop_stats(x_dist::Matrix{Float64}, species_col::Int64, params::SimParams)
    @unpack no_param = params
    
    # Filter data for the specific species
    species_data = x_dist[x_dist[:, 1] .== species_col, 2:no_param+1]
    
    # Handle the case of no individuals
    if isempty(species_data)
        return fill(NaN, 1, no_param), fill(NaN, 1, no_param)
    end

    # Calculate mean and variance using mapslices
    # A more type-stable approach would be to use a package like StatsBase.jl
    pop_mean = mapslices(x -> mean(skipmissing(x)), species_data, dims=1)
    
    # Variance is zero for a single individual
    if size(species_data, 1) == 1
        pop_var = zeros(1, no_param)
    else
        pop_var = mapslices(x -> var(skipmissing(x)), species_data, dims=1)
    end
    
    return pop_mean, pop_var
end

# A dedicated function for log-normal distribution parameters
function get_lognormal_params(mean::Number, stdev::Number)
    if iszero(stdev)
        return log(mean), 0.0
    end
    
    mu = log(mean^2 / sqrt(stdev^2 + mean^2))
    sigma = sqrt(log(stdev^2 / mean^2 + 1))
    
    return mu, sigma
end

using Distributions, Statistics, NaNMath
using Base.Threads
using UnPack

function DrawNewTraits(x_dist::Matrix{Float64}, parent_traits::Vector{Float64},
                       species_col::Int64, j::Int64, params::SimParams)
    @unpack h2_vect, no_param, no_columns = params
    
    # 1. Get Population Stats
    pop_mean, pop_var = get_pop_stats(x_dist, species_col, params)
    pop_stdev = sqrt.(pop_var)
    
    # 2. Calculate Offspring Trait Parameters (using broadcasting)
    h2 = h2_vect[j, species_col]
    
    # Use broadcasting for element-wise operations on arrays
    exp_offspring_traits = h2 .* parent_traits[1:no_param]' .+ (1 - h2) .* pop_mean
    sigma_offspring_traits = sqrt(1 - h2^2) .* pop_stdev
    
    # 3. Draw New Quantitative Traits
    offspring_traits = Vector{Float64}(undef, no_param)
    @threads for i in 1:no_param
        mu, sigma = get_lognormal_params(exp_offspring_traits[i], sigma_offspring_traits[i])
        if !isnan(mu) && !iszero(sigma)
            offspring_traits[i] = rand(LogNormal(mu, sigma))
        else
            offspring_traits[i] = exp_offspring_traits[i]
        end
    end

    # 4. Handle Offspring Genotypes
    offspring_genotypes = parent_traits[no_param+1:no_columns-1]

    return offspring_traits, offspring_genotypes
end

"""