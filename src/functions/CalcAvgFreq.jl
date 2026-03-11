
# ========================================== #
#		FUNCTION CALC AVG FREQUENCIES        #
# ========================================== #

function CalcMedian(ii::Int, no_columns::Int, no_param::Int, x_dist::Matrix{Float64})
    if !isempty(x_dist) && count(x_dist[:, 1] .== ii) > 0 
		# get the rows for the state of interest
		qt_medians = mapslices(median ∘ skipmissing, x_dist[x_dist[:, 1] .== ii, 2:no_param+1], dims=1)
		# Step 2: Discrete trait frequencies
		# take the sum of the genotype cols 
		gt_freqs = sum(x_dist[x_dist[:, 1] .== ii,2+no_param:no_columns], dims=1) ./ sum(x_dist[:, 1] .== ii)	
		# Step 3: Combine the results
		median_freqs = hcat(qt_medians, gt_freqs)
	return median_freqs
	end
end

function CalcMean(ii::Int, no_columns::Int, no_param::Int, x_dist::Matrix{Float64})
	if !isempty(x_dist) && count(x_dist[:, 1] .== ii) > 0
		# Step 1: Quantitative trait medians (ignoring NaN values)
		qt_means = mapslices(mean ∘ skipmissing, x_dist[x_dist[:, 1] .== ii, 2:no_param+1], dims=1)
		# Step 2: Discrete trait frequencies
		# take the sum of the genotype cols 
		gt_freqs = sum(x_dist[x_dist[:, 1] .== ii,2+no_param:no_columns], dims=1) ./ sum(x_dist[:, 1] .== ii)
			# Step 3: Combine the results
		mean_freqs = hcat(qt_means, gt_freqs)
	
	return mean_freqs
	#else
	#	return nothing # IDEA? if there are no individuals in the state, return nothing
	end
end

function CalcVar(ii::Int, no_param::Int, x_dist::Matrix{Float64})

	temp_size = size(x_dist[x_dist[:,1] .== ii,2:no_param+1])
	if temp_size[1] == 1 #has only one row = only one individual
		trait_var = zeros(Float64, temp_size[1], length(2:no_param+1))
	else 
		# get variance
		trait_var = mapslices(var ∘ skipmissing, x_dist[x_dist[:,1] .== ii,2:no_param+1],dims=1)
	end
	return trait_var
end
