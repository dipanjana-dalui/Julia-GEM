


##############################################
#		FUNCTION CALC AVG FREQUENCIES        #
##############################################
function CalcMedian(ii::Int64, no_columns::Int64, no_param::Int64, x_dist::Matrix{Float64})
	# Step 1: Quantitative trait medians (ignoring NaN values)
	qt_medians = mapslices(median ∘ skipmissing, x_dist[x_dist[:, 1] .== ii, 2:no_param+1], dims=1)
	# Step 2: Discrete trait frequencies
	# take the sum of the genotype cols 
	gt_freqs = sum(x_dist[x_dist[:, 1] .== ii,2+no_param:no_columns], dims=1) ./ sum(x_dist[:, 1] .== ii)
	# Step 3: Combine the results
	median_freqs = hcat(qt_medians, gt_freqs)
	return median_freqs
end

function CalcMean(ii::Int64, no_columns::Int64, no_param::Int64, x_dist::Matrix{Float64})
	# Step 1: Quantitative trait medians (ignoring NaN values)
	qt_means = mapslices(mean ∘ skipmissing, x_dist[x_dist[:, 1] .== ii, 2:no_param+1], dims=1)
	# Step 2: Discrete trait frequencies
	# take the sum of the genotype cols 
	gt_freqs = sum(x_dist[x_dist[:, 1] .== ii,2+no_param:no_columns], dims=1) ./ sum(x_dist[:, 1] .== ii)
	# Step 3: Combine the results
	mean_freqs = hcat(qt_means, gt_freqs)
	return mean_freqs
end

function CalcVar(ii::Int64, no_param::Int64, x_dist::Matrix{Float64})
	temp_size = size(x_dist[x_dist[:,1] .== ii,2:no_param+1])
	if temp_size[1] == 1 #has only one row = only one individual
		trait_var = Matrix{Float64}(fill(0, temp_size[1], length(2:no_param+1)))
	else
		trait_var = mapslices(var ∘ skipmissing, x_dist[x_dist[:,1] .== ii,2:no_param+1],dims=1)
	end
	return trait_var
end

"""
MoJuWo below
"""


# Combine shared parameters into a single struct
# (This would be in your params.jl file)
"""
struct SimParams
    no_columns::Int
    no_param::Int
end
"""

# --- Adaptations for CalcMedian and CalcMean ---

function CalcStatistics(
    stats_func,
    ii::Int,
    x_dist::AbstractMatrix{<:Number},
    params::SimParams
)
    # The `params` struct holds the shared dimensions
    @unpack no_columns, no_param = params

    # Step 1: Filter the matrix for the current species
    species_data = x_dist[x_dist[:, 1] .== ii, :]
    
    # If no individuals of this species exist, return a vector of NaNs
    if isempty(species_data)
        return fill(NaN, 1, no_columns - 1)
    end

    # Step 2: Calculate quantitative trait statistics
    qt_data = @view species_data[:, 2:no_param+1]
    qt_stats = mapslices(stats_func ∘ skipmissing, qt_data, dims=1)

    # Step 3: Calculate discrete trait frequencies
    gt_data = @view species_data[:, no_param+2:end]
    gt_freqs = sum(gt_data, dims=1) ./ size(gt_data, 1)

    # Step 4: Combine and return
    return hcat(qt_stats, gt_freqs)
end

# Now define your specific functions using the generic one
CalcMedian(ii::Int, x_dist::AbstractMatrix{<:Number}, params::SimParams) = CalcStatistics(median, ii, x_dist, params)
CalcMean(ii::Int, x_dist::AbstractMatrix{<:Number}, params::SimParams) = CalcStatistics(mean, ii, x_dist, params)


# --- Adaptation for CalcVar ---

function CalcVar(
    ii::Int,
    x_dist::AbstractMatrix{<:Number},
    params::SimParams
)
    @unpack no_param = params

    # Step 1: Filter for the species' individuals
    species_data = x_dist[x_dist[:, 1] .== ii, :]
    
    # Check for empty or single-individual cases
    if size(species_data, 1) <= 1
        return fill(0.0, 1, no_param)
    else
        # Step 2: Calculate variance for quantitative traits
        qt_data = @view species_data[:, 2:no_param+1]
        return mapslices(var ∘ skipmissing, qt_data, dims=1)
    end
end

"""
call from main GEM sim function will look something like ~
	# The ii argument is defined in this scope, e.g., within a loop
    for ii in 1:params.no_species
        # Call the modernized functions
        median_result = CalcMedian(ii, x_dist, params)
        mean_result = CalcMean(ii, x_dist, params)
        var_result = CalcVar(ii, x_dist, params)
        
        # ... do something with the results ...

"""