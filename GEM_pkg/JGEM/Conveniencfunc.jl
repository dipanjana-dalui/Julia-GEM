## Function to convert 0s to NaN in matrix
function zero_to_nan(Bool_mat::Matrix{Int64})
	return [x==0 ? NaN : x for x in Bool_mat]
end

## function to set a threshold for zero
function deci_threshold(df::DataFrame)
	#trait = trait_mean_df.b_max
	#round.(trait_mean_df.b_max, digits=4)

	df_mod = mapcols(col -> round.(col, digits=6), df)
	return df_mod
end

"""
MoJuWo adapted
"""
function zeros_to_nan(M::AbstractMatrix{<:Number})
    # Create a new matrix of the correct type to hold the results
    # The `convert` call ensures type stability for the output
    result = similar(M, Union{eltype(M), Float64})
    
    # Use broadcasting for an efficient element-wise operation
    result .= ifelse.(M .== 0, NaN, M)
    
    return result
end

# An even more concise version using a single broadcasted operation
zeros_to_nan_concise(M::AbstractMatrix{<:Number}) = ifelse.(M .== 0, NaN, Float64.(M))

function round_df(df::DataFrame)
    # The `!` indicates an in-place modification
    df_mod = copy(df) # Create a copy to avoid modifying the original DataFrame
    for col in eachcol(df_mod)
        col .= round.(col, digits=6)
    end
    return df_mod
end